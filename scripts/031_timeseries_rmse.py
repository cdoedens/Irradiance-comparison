#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import datetime
import xesmf as xe

from dask.distributed import Client, wait
import os, sys

sys.path.append('/home/548/cd3022/repos/Irradiance-comparisons/Irradiance-comparisons')
import logger
from read_datasets import read_dataset

year = sys.argv[1]
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
        n_workers=12,
        threads_per_worker=1,
    )
    ds1_list = []
    ds2_list = []
    for month in range(1,13):
        date = f'{year}-{month:02d}'
        ds1 = read_dataset(
                dataset='himawari',
                resolution='daily',
                date=date
            )
        ds1_list.append(ds1)
        ds2 = read_dataset(
                dataset='barra-r2',
                resolution='daily',
                date=date
            )
        ds2_list.append(ds2)
    ds1 = xr.concat(ds1_list, dim='time')
    ds2 = xr.concat(ds2_list, dim='time')

    ds1 = ds1.assign_coords(time=ds1.time + np.timedelta64(1, "D"))
    ds1 = ds1.assign_coords(time=ds1.time.dt.floor("D"))
    ds2 = ds2.assign_coords(time=ds2.time.dt.floor("D"))

    regridder_file = '/g/data/er8/users/cd3022/regridder_weights/himawari_to_barrar2_weights.nc'
    regridder =  xe.Regridder(ds1, ds2, "bilinear",
                        filename=regridder_file,
                        reuse_weights=True)
    ds1 = regridder(ds1)
    rmse = rmse_workflow(ds1, ds2)

    file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/error_timeseries/')
    os.makedirs(file_path, exist_ok=True)
    rmse.to_netcdf(f'{file_path}/himawari-barrar2_{year}.nc')
