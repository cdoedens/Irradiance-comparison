#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import logger
import datetime

from read_datasets import read_dataset

from dask.distributed import Client, as_completed
import os, sys

LOG = logger.get_logger(__name__)
# qsub -I -q normal -P er8 -l walltime=1:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74
def get_data(date):
    ds1= read_dataset(
            dataset='himawari',
            resolution='daily',
            date=date
        )
    ds2= read_dataset(
            dataset='barra-r2',
            resolution='daily',
            date=date
        )
    LOG.info('datasets opened')

    # Adjsust times so they match
    # e.g. shift Dataset A forward by 1 day so its timestamps match Dataset B
    ds1 = ds1.assign_coords(time=ds1.time + np.timedelta64(1, "D"))
    ds1 = ds1.assign_coords(time=ds1.time.dt.floor("D"))
    ds2 = ds2.assign_coords(time=ds2.time.dt.floor("D"))

    # regrid fine resolution to coarse resolution so datasets match
    if ds1.ghi.shape[1] > ds2.ghi.shape[1]:
        ds1 = ds1.interp(
            lat=ds2.lat,
            lon=ds2.lon,
            method='nearest' # 'linear', 'nearest', 'cubic'
        )
    else:
        ds2 = ds2.interp(
            lat=ds1.lat,
            lon=ds1.lon,
            method='linear'
        )
    LOG.info('dataset times adjusted')
    return ds1.ghi, ds2.ghi


def calc_rmse_timeseries(diff):
    LOG.info('start calc rmse')
    rmse = np.sqrt((np.square(diff)).weighted(np.cos(np.deg2rad(diff.lat))).mean(['lat', 'lon']))
    rmse.name = 'rmse'
    return rmse

if __name__ == '__main__':
    client = Client(
        n_workers = 24,
        memory_limit = int(os.environ['PBS_VMEM']) / int(os.environ['PBS_NCPUS']),
    )
    futures = []
    for month in range(1,13):
        for year in range(2016, 2025):
            date = f'{year}-{month:02d}'
            ds1, ds2 = get_data(date)
            diff =  ds2 - ds1
            rmse = client.submit(calc_rmse_timeseries, diff)
            futures.append(rmse)
        
    results = []
    for future in as_completed(futures):
        res = future.result()
        LOG.info('rmse calculated')
        results.append(res)
    combined = xr.concat(results, dim='time')
    print(f"final data shape: {combined.shape}')

