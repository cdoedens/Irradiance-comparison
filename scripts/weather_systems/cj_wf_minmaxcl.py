import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import datetime
import os, sys

from read_datasets import read_dataset
import logger

LOG = logger.get_logger(__name__)

from dask.distributed import Client


if __name__ == '__main__':
    
    client = Client(
        n_workers=12,
        threads_per_worker=1
    )

    dataset = sys.argv[1]
    wf = sys.argv[2]
    year = sys.argv[3]

    LOG.info(f'Starting analysis for dataset {dataset}, weather feature {wf}, year {year}')


    
    for month in range(1, 13):
        month = f'{month:02d}'
        # LOAD DATASETS
        # GHI dataset
        ds_ghi = read_dataset(
                dataset=dataset,
                resolution='hourly',
                date=f'{year}-{month}'
            )
        LOG.info(f'{dataset} data opened')
    


        # Weather feature dataset
        mask_var_names = {
            'maxcl': 'FLAG',
            'mincl': 'INPUT',
        }
        mask_var = mask_var_names[wf]
        def preprocess(_ds):
            # assign lat/lon values if missing
            if wf == 'maxcl':
                lat_vals = np.linspace(-90, 90, 361)
                lon_vals = np.linspace(-180, 179.75, 720)
                
                _ds = _ds.assign_coords(
                    lat=lat_vals,
                    lon=lon_vals
                )
            return _ds.sel(
                lat=slice(ds_ghi.lat.min(), ds_ghi.lat.max()),
                lon=slice(ds_ghi.lon.min(), ds_ghi.lon.max())
            )
        wf_path = Path(f'/g/data/su28/weatherfeatures.era5/{wf}/')
        wf_files = [f for f in wf_path.rglob(f'*{year}_{month}*.nc')]
        ds_wf = xr.open_mfdataset(wf_files) #, preprocess=preprocess)
        if wf == 'mincl':
            ds_wf = ds_wf.isel(dimz_INPUT=0)

        # PREPROCESS DATASETS BEFORE MASKING
        # Interpolate grid cells
        ds_wf = ds_wf.interp(
            lat=ds_ghi.lat,
            lon=ds_ghi.lon,
            method='nearest'
        )

        # adjust hourly times on the 30min to match 3hr times (on the hour) from wf
        ds_ghi = ds_ghi.assign_coords(time=ds_ghi.time + pd.Timedelta('30min'))

        # Extra Himawari preprocessing steps
        if dataset == 'himawari':
            start = ds_ghi.time.min().item()
            end = ds_ghi.time.max().item()
            # fill missing overnight time steps
            full_time = pd.date_range(
                start=start,
                end=end,
                freq="60min"
            )
            ds_ghi = ds_ghi.reindex(time=full_time, fill_value=0)
            
            # remove timesteps from the end of wf_ds to match himawari,
            # missing timesteps are overnight so will not impact final data
            # ds_wf = ds_wf.isel(time=slice(0, -14))

            # line up datasets
            min_time = ds_wf.time.min()
            max_time = ds_ghi.time.max()
            t_range = slice(min_time, max_time)
            ds_wf = ds_wf.sel(time=t_range)
            ds_ghi = ds_ghi.sel(time=t_range)

        ds_ghi_times = ds_ghi.sel(time=ds_wf.time)

        LOG.info('preprocessing complete')

        # APPLY MASK
        mask = ds_wf[mask_var]
        masked_ghi = xr.where(mask != 0, ds_ghi_times.ghi, np.nan)

        # TIME MEAN
        masked_ghi_mean = masked_ghi.mean(dim='time')

        # remask himawari region
        if dataset == 'himawari':
            masked_ghi_mean = xr.where(ds_ghi.isel(time=5).ghi.isnull(), np.nan, masked_ghi_mean)


        # add metadata
        masked_ghi_mean = masked_ghi_mean.to_dataset(name='monthly_mean_ghi')
        masked_ghi_mean = masked_ghi_mean.assign_coords(date=f'{year}-{month}')
        masked_ghi_mean.attrs['date_generated'] = datetime.date.today().strftime('%D')
        masked_ghi_mean.attrs['source_script'] = 'data produced by the script "cj_wf_minmaxcl.py"'

        # SAVE DATA
        file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/weather-features-cj/{wf}')
        os.makedirs(file_path, exist_ok=True)
        masked_ghi_mean.to_netcdf(f'{file_path}/{dataset}-{year}{month}.nc')
        LOG.info('Month data saved')
    LOG.info('Job complete!')

    

    

    