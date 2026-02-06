#!/g/data/xp65/public/apps/med_conda_scripts/analysis3-25.07.d/bin/python3
import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import datetime
import xesmf as xe
import re

from dask.distributed import Client, wait
import os, sys
import psutil

sys.path.append('/home/548/cd3022/repos/Irradiance-comparisons/Irradiance-comparisons')
import logger
from read_datasets import read_dataset

LOG = logger.get_logger(__name__)

if __name__ == '__main__':

    client = Client(
        n_workers=48,
        threads_per_worker=1
    )


    # BARRA Cloud Types
    BARRA = Path('/g/data/ob53/BARRA2/output/reanalysis/')
    BARRA_R2 = BARRA / "AUS-11/BOM/ERA5/historical/hres/BARRA-R2/v1"
    BARRA_R2_DIR = BARRA_R2 / "1hr"

    vars_to_open = [
        'clh',
        'clm',
        'cll'
    ]
    def _preprocess(ds):
        return ds.sel(lat=slice(-44.5, -10), lon=slice(112, 156.26), time=slice('2016', '2024'))

    ds_vars = []
    for var in vars_to_open:

        files = sorted([f for f in BARRA_R2_DIR.glob(f'{var}/latest/*1hr_20*.nc')])
        ds_b = xr.open_mfdataset(files,
                            chunks={'time':24, 'lat':157, 'lon':31},
                            concat_dim='time',
                            combine='nested',
                            data_vars='minimal',
                            coords='minimal',
                            compat='override',
                            parallel=True,
                            preprocess=_preprocess
                            )
        ds_vars.append(ds_b)
    barra_r2 = xr.merge(ds_vars)
    LOG.info('BARRA-R2 opened')

    # NEED TO CONTINUE TO REPLICATE WORKFLOW FROM HIM_CT.py
    # ...AND THEN CONTINUE FOR ERA5 AND BARRA-C2