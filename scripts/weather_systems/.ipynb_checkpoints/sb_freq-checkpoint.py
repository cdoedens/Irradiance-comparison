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

import logger
from read_datasets import read_dataset

'''
GET THE SPATIAL FREQ OF SEA BREEZE OBJECTS
'''

if __name__ == '__main__':

    client = Client(
        n_workers=48,
        threads_per_worker=1
    )

    sb_years = []
    for year in range(2016, 2025):
        sb_path = Path('/g/data/ng72/ab4502/sea_breeze_detection/barra_c_smooth_s2/filters/')
        sb_files = [f for f in sb_path.glob(f'*_F_{year}*.zarr')]
        ds_sb = xr.open_mfdataset(
            sb_files,
            engine="zarr",
            combine="by_coords",
            compat='override',
            coords='minimal',
            parallel=True,
        )

        mask = ds_sb.mask.fillna(False).astype(bool)
        sb_extension = 2 # hours
        sb_extended = (
            mask
            | mask.shift(time=sb_extension).fillna(False)
            | mask.shift(time=-sb_extension).fillna(False)
        )

        counts = xr.where(sb_extended != 0, 1, 0).sum(dim='time')
        counts = counts.assign_coords(year=year)
        sb_years.append(counts)
    all_counts = xr.concat(sb_years, dim='year')
    total_freq = all_counts.sum(dim='year')
    total_freq.to_netcdf('/g/data/er8/users/cd3022/Irradiance-comparisons/weather-features/sea_breeze/TOTAL_FREQ.nc')
        