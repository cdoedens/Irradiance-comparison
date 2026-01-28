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

LOG = logger.get_logger(__name__)

year = sys.argv[1]

# qsub -I -q normal -P er8 -l walltime=2:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+gdata/su28
    
if __name__ == '__main__':
    client = Client(
        n_workers=12,
        threads_per_worker=1,
    )
    LOG.info('dask client started')

    # OPEN DATA
    him_list = []
    bar_list = []
    for month in range(1,13):
        date = f'{year}-{month:02d}'
        him = read_dataset(
                dataset='himawari',
                resolution='daily',
                date=date
            )
        him_list.append(him)
        bar = read_dataset(
                dataset='barra-r2',
                resolution='daily',
                date=date
            )
        bar_list.append(bar)
    him = xr.concat(him_list, dim='time')
    bar = xr.concat(bar_list, dim='time')
    LOG.info('Datasets opened')

    # PREPROCESS
    him = him.assign_coords(time=him.time + np.timedelta64(1, "D"))
    him = him.assign_coords(time=him.time.dt.floor("D"))
    bar = bar.assign_coords(time=bar.time.dt.floor("D"))

    him = him.interp(
        lat=bar.lat,
        lon=bar.lon,
        method='linear'
    )
    LOG.info('Data preprocessed')

    diff = bar.ghi - him.ghi

    rmse = np.sqrt((np.square(diff)).weighted(np.cos(np.deg2rad(diff.lat))).mean(['lat', 'lon']))
    mbe = diff.weighted(np.cos(np.deg2rad(diff.lat))).mean(['lat', 'lon'])
    LOG.info('Error metrics calculated')

    ds_errors = xr.Dataset(
        {
            'rmse': rmse,
            'mbe': mbe
        }
    )

    ds_errors.attrs['history'] = datetime.date.today().strftime('%D')
    ds_errors.attrs['source_script'] = 'data produced by the script "031_timeseries_errors.py"'
    LOG.info('Metadata added')

    LOG.info('Loading data to memory...')
    ds_errors.load()

    file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/error_timeseries/')
    os.makedirs(file_path, exist_ok=True)
    ds_errors.to_netcdf(f'{file_path}/himawari-barrar2_{year}.nc')
    LOG.info('Files saved, job complete!')