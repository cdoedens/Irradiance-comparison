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
GET THE SPATIAL FREQ OF CLOUD TYPES
'''

LOG = logger.get_logger(__name__)

cloud_aggs = {
    'clear_sky':[1,2,3,4],
    'low':[5,6],
    'med':[7], 
    'high_opaque':[8,9],
    'high_semitransparent': [10,11,12,13,14] 
}

if __name__ == '__main__':

    client = Client(
        n_workers=48,
        threads_per_worker=1
    )

    ct = sys.argv[1]

    ds = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")
    ds = ds.sel(
        lat=slice(-10, -45),
        lon=slice(112, 156.4),
        time=slice('2016', '2024')
    )
    LOG.info('Data opened')
    # ct_years = []
    
    # TEST ALL IN ONE GO
    # for year in range(2016, 2025):
        # LOG.info(f'Start year: {year}')

    ct_codes = cloud_aggs[ct]
    ct_mask = xr.where(ds.ct.isin(ct_codes), 1, np.nan)
    LOG.info('CT aggreagated into mask')

    counts = ct_mask.sum(dim='time')
    LOG.info('Counts summed')
        
        # counts = counts.assign_coords(year=year)
        # sb_years.append(counts)
    # all_counts = xr.concat(sb_years, dim='year', coords='minimal')
    # total_freq = all_counts.sum(dim='year')
    counts.to_netcdf(f'/g/data/er8/users/cd3022/Irradiance-comparisons/cloud_types/TOTAL_FREQ_{ct}.nc')
    LOG.info('Job complete!')

        