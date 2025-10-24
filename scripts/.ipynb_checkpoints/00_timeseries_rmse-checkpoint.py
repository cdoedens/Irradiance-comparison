#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logger
import datetime
import xesmf as xe

from read_datasets import read_dataset

from dask.distributed import Client, wait
import os, sys

year_month = sys.argv[1]
LOG = logger.get_logger(__name__)
# qsub -I -q normal -P er8 -l walltime=2:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+gdata/su28

def get_diff(ds1, ds2):
    return ds2.ghi - ds1.ghi

def calc_rmse(diff):
    LOG.info('start calc rmse')
    rmse = np.sqrt((np.square(diff)).weighted(np.cos(np.deg2rad(diff.lat))).mean(['lat', 'lon']))
    LOG.info('rmse calculated')
    rmse.name = 'rmse'
    return rmse

def rmse_workflow(ds1, ds2):
    diff = get_diff(ds1, ds2)
    LOG.info('diff calculated')
    return calc_rmse(diff)
    
if __name__ == '__main__':
    client = Client(
        n_workers=6,
        threads_per_worker=4,
        processes=False,
        memory_limit = int(os.environ['PBS_VMEM']) / int(os.environ['PBS_NCPUS']),
    )
    # OPEN DATA
    ds1= read_dataset(
            dataset='himawari',
            resolution='daily',
            date=year_month
        )
    ds2= read_dataset(
            dataset='barra-r2',
            resolution='daily',
            date=year_month
        )
    LOG.info('datasets opened')
    
    # Adjsust times so they match
    ds1 = ds1.assign_coords(time=ds1.time + np.timedelta64(1, "D"))
    ds1 = ds1.assign_coords(time=ds1.time.dt.floor("D"))
    ds2 = ds2.assign_coords(time=ds2.time.dt.floor("D"))
    LOG.info('dates aligned')

    # regrid fine resolution to coarse resolution so datasets match
    regridder_file = '/g/data/er8/users/cd3022/regridder_weights/himawari_to_barrar2_weights.nc'
    regridder =  xe.Regridder(ds1, ds2, "bilinear",
                         filename=regridder_file,
                         reuse_weights=True)
    ds1 = regridder(ds1)
    LOG.info('lat/lon regridded')

    file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/error_timeseries/')
    os.makedirs(file_path, exist_ok=True)
    for day in range(1,32):
        
        date = f'{year_month}-{day:02d}'
        ds1_day = ds1.sel(time=date)
        ds2_day = ds2.sel(time=date)
        if len(ds1_day) == 0:
            break
        rmse = client.submit(rmse_workflow(date))
        LOG.info('rmse calculated')



    rmse.to_netcdf(f'{file_path}/himawari-barrar2_{date}.nc')
    LOG.info(f'{date} data saved')
