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

year = sys.argv[1]

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
    return ds1.ghi - ds2.ghi


def get_ct(date):
    ds = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")
    ds = ds.sel(time=date)
    daily = ds.sel(time=ds.time.dt.floor("D"))

if __name__ == '__main__':
    client = Client(
        n_workers = 24,
        memory_limit = int(os.environ['PBS_VMEM']) / int(os.environ['PBS_NCPUS']),
    )

