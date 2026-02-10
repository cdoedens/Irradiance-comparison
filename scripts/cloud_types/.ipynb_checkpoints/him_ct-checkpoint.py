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

dataset = sys.argv[1]
ct = sys.argv[2]
year = sys.argv[3]

if dataset == 'himawari':
    resolution = 'instant'
else:
    resolution = 'hourly'

if __name__ == '__main__':

    client = Client(
        n_workers=48,
        threads_per_worker=1
    )

    ds = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")

    for i, month in enumerate(range(1, 13)):

        month = f'{month:02d}'

        ds_ghi = read_dataset(
            dataset=dataset,
            resolution=resolution,
            date=f'{year}-{month}'
        )

        if i == 0:
            ds = ds.sel(
                lat=slice(ds_ghi.lat.max(), ds_ghi.lat.min()),
                lon=slice(ds_ghi.lon.min(), ds_ghi.lon.max())
            )
        him_ct = ds.ct.sel(time=f'{year}-{month}')

        if dataset == 'himawari':
            # ALIGN DATASETS
            ds_ghi = ds_ghi.interp(
                lat=him_ct.lat,
                lon=him_ct.lon,
                method='linear'
            )
            start = ds_ghi.time.min().item()
            end = ds_ghi.time.max().item()
            # fill missing overnight time steps
            full_time = pd.date_range(
                start=start,
                end=end,
                freq="10min"
            )
            ds_ghi = ds_ghi.reindex(time=full_time, fill_value=0)
        elif dataset == 'barra-r2':
            # ALIGN DATASETS
            him_ct = him_ct.interp(
                lat=ds_ghi.lat,
                lon=ds_ghi.lon,
                method='nearest'
            )
        
        # # trim ends
        # min_time = him_ct.time.min()
        # max_time = ds_ghi.time.max()
        # t_range = slice(min_time, max_time)
        # him_ct = him_ct.sel(time=t_range)
        # ds_ghi = ds_ghi.sel(time=t_range)

        t1 = ds_ghi.time.values
        t2 = him_ct.time.values
        
        good_t = np.intersect1d(t1, t2)
        
        him_ct = him_ct.sel(time=good_t)
        ds_ghi = ds_ghi.sel(time=good_t)

        # AGGREGATE CT FOR MASK
        cloud_aggs = {
            'clear_sky':[1,2,3,4],
            # 'low':[5,6],
            # 'med':[7], 
            'low_med':[5,6,7],
            'high_opaque':[8,9],
            'high_semitransparent': [10,11,12,13,14] 
        }

        # for i, ct in enumerate(cloud_aggs):
        # REMOVED LOOP AND INSTEAD DOING ONE CT PER JOB
        ct_codes = cloud_aggs[ct]
        ct_ghi = xr.where(him_ct.isin(ct_codes), ds_ghi.ghi, np.nan)

        # CT DATA
        ct_ghi_monthly_mean = ct_ghi.mean(dim='time')
        # remask himawari region
        if dataset == 'himawari':
            ct_ghi_monthly_mean = xr.where(ds_ghi.isel(time=5).ghi.isnull(), np.nan, ct_ghi_monthly_mean)

        # add metadata
        ct_ghi_monthly_mean = ct_ghi_monthly_mean.to_dataset(name=f'{ct}_ghi')
        ct_ghi_monthly_mean = ct_ghi_monthly_mean.assign_coords({'date':f'{year}-{month}'})
        ct_ghi_monthly_mean.attrs['date_generated'] = datetime.date.today().strftime('%D')
        ct_ghi_monthly_mean.attrs['source_script'] = 'data produced by the script "him_ct.py"'

        # SAVE DATA
        file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/cloud_types/{dataset}/')
        os.makedirs(file_path, exist_ok=True)
        ct_ghi_monthly_mean.to_netcdf(f'{file_path}/{ct}-{year}-{month}.nc')
        LOG.info(f'CT {ct} data saved')

        LOG.info(f'Month {month} finished')
    LOG.info('Job complete!')


