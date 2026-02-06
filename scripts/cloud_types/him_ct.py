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

year = sys.argv[1]

if __name__ == '__main__':

    client = Client(
        n_workers=48,
        threads_per_worker=1
    )

    ds = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")

    for month in range(1, 13):
        month = f'{month:02d}'

        ds_ghi = read_dataset(
            dataset='himawari',
            resolution='instant',
            date=f'{year}-{month}'
        )

        him_ct = ds.ct.sel(time=f'{year}-{month}')

        cloud_aggs = {
            'clear_sky':[1,2,3,4],
            'low':[5,6],
            'med':[7], 
            'high_opaque':[8,9],
            'high_semitransparent': [10,11,12,13,14] 
        }

        for i, ct_agg in enumerate(cloud_aggs):
            ct_codes = cloud_aggs[ct_agg]
            ct_ghi = xr.where(him_ct.isin(ct_codes), ds_ghi.ghi, np.nan)
            ct_ghi_monthly_mean = ct_ghi.mean(dim='time')

            # add metadata
            ct_ghi_monthly_mean = ct_ghi_monthly_mean.to_dataset(name=f'{ct_agg}_ghi')
            ct_ghi_monthly_mean = ct_ghi_monthly_mean.assign_coords({'date':f'{year}-{month}'})
            ct_ghi_monthly_mean.attrs['date_generated'] = datetime.date.today().strftime('%D')
            ct_ghi_monthly_mean.attrs['source_script'] = 'data produced by the script "him_ct.py"'

            # SAVE DATA
            file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/cloud_types/')
            os.makedirs(file_path, exist_ok=True)
            ct_ghi_monthly_mean.to_netcdf(f'{file_path}/{ct_agg}-{year}-{month}.nc')
            LOG.info(f'CT {ct_agg} data saved')

        LOG.info(f'Month {month} finished')
    LOG.info('Job complete!')


