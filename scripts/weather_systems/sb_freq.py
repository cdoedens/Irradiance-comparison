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

LOG = logger.get_logger(__name__)

year = sys.argv[1]

month_list = []
freq = 0
for month in range(1, 13):
    month=f'{month:02d}'
    sb_path = Path('/g/data/ng72/ab4502/sea_breeze_detection/barra_c_smooth_s2/filters/')
    sb_files = [f for f in sb_path.glob(f'*_F_{year}{month}*.zarr')][0]
    ds_sb = xr.open_zarr(sb_files)
    
    mask = ds_sb.mask
    sb_extended = mask.rolling(time = 5, center=True).max()
    LOG.info('Mask applied and extended')
    
    counts = xr.where(sb_extended != 0, 1, 0).sum(dim='time')
    counts = counts.assign_coords(month=month)
    freq += ds_sb.sizes['time']
    month_list.append(counts)
all_mon = xr.concat(month_list, dim='month')
year_freq = all_mon.sum(dim='month')
year_freq = year_freq.assign_coords(year=year, freq=freq)

year_freq.to_netcdf(f'/g/data/er8/users/cd3022/Irradiance-comparisons/weather-features/sea_breeze/TOTAL_FREQ_{year}.nc')
LOG.info('Job complete!')