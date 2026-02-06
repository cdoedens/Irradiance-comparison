import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import sys, os
import datetime

from read_datasets import read_dataset
import logger

LOG = logger.get_logger(__name__)

from dask.distributed import Client

if __name__ == '__main__':

    client = Client(
        n_workers=6,
        threads_per_worker=1
    )

    dataset = sys.argv[1]
    year = sys.argv[2]
    LOG.info(f'Starting analysis for dataset {dataset}, sea breezes, year {year}')

    for month in range(1, 13):
        month = f'{month:02d}'
        # LOAD DATASETS
        ds = read_dataset(
                dataset=dataset,
                resolution='hourly',
                date=f'{year}-{month}'
            )
        LOG.info(f'{dataset} data opened')

        sb_path = Path('/g/data/ng72/ab4502/sea_breeze_detection/barra_c_smooth_s2/filters/')
        sb_files = [f for f in sb_path.glob(f'*_F_{year}{month}*.zarr')]
        ds_sb = xr.open_mfdataset(
            sb_files,
            engine="zarr",
            combine="by_coords",
            compat='override',
            coords='minimal',
            parallel=True,
        )
        LOG.info('sea breeze data opened')

        # PREPROCESS DATASETS TO MATCH TIMES/SHAPES
        # interp wf mask from era5 to barra
        ds_sb = ds_sb.interp(
            lat=ds.lat,
            lon=ds.lon,
            method='nearest'
        )
        
        # adjust hourly times on the 30min to match 3hr times (on the hour) from wf
        ds_shifted = ds.assign_coords(time=ds.time + pd.Timedelta('30min'))
        
        if dataset == 'himawari':
            # fill missing overnight time steps
            full_time = pd.date_range(
                start=ds_shifted.time.min().item(),
                end=ds_shifted.time.max().item(),
                freq="60min"
            )
            ds_shifted = ds_shifted.reindex(time=full_time, fill_value=0)
            
            # line up datasets
            min_time = ds_sb.time.min()
            max_time = ds_shifted.time.max()
            t_range = slice(min_time, max_time)
            ds_sb = ds_sb.sel(time=t_range)
            ds = ds.sel(time=t_range)

        ds_times = ds_shifted.sel(time=ds_sb.time)
        LOG.info('preprocessing complete')

        # Extend mask to capture larger area round sea breeze
        mask = ds_sb.mask.fillna(False).astype(bool)
        sb_extension = 2 # hours
        sb_extended = (
            mask
            | mask.shift(time=sb_extension).fillna(False)
            | mask.shift(time=-sb_extension).fillna(False)
        )
        LOG.info(f'sb mask extended forward and back {sb_extension} hours')
        
        # APPLY MASK
        sb_ghi = xr.where(sb_extended != 0, ds_times.ghi, np.nan)
        LOG.info('mask applied')

        # TIME MEAN
        sb_ghi_mean = sb_ghi.mean(dim='time')
        LOG.info('annual mean calculated')

        # remask himawari region
        if dataset == 'himawari':
            sb_ghi_mean = xr.where(ds.isel(time=5).ghi.isnull(), np.nan, sb_ghi_mean)

        # add metadata
        sb_ghi_mean = sb_ghi_mean.to_dataset(name='annual_mean_ghi')
        sb_ghi_mean = sb_ghi_mean.assign_coords({'year':year})
        sb_ghi_mean.attrs['date_generated'] = datetime.date.today().strftime('%D')
        sb_ghi_mean.attrs['source_script'] = 'data produced by the script "sea_breezes.py"'

        # SAVE DATA
        file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/weather-features/sea_breaze')
        os.makedirs(file_path, exist_ok=True)
        sb_ghi_mean.to_netcdf(f'{file_path}/{dataset}-{year}-{month}.nc')
        LOG.info('Monthly data saved')
    LOG.info('Job complete!')

















# sb_extended = ds_sb.mask.rolling(time=3, center=True).max()

# mask = ds_sb.mask.fillna(False).astype(bool)
# sb_extended = (
#     mask
#     | mask.shift(time=2).fillna(False)
#     | mask.shift(time=-2).fillna(False)
# ).compute()
# data = sb_extended.sel(time='2018-01-01').sum(dim='time').astype(float)
# data = ds_sb.mask.sel(time='2018-01-01').sum(dim='time')

# import matplotlib.pyplot as plt
# import cartopy.crs as ccrs
# fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
# mesh = ax.pcolormesh(
#         ds_sb.lon,
#         ds_sb.lat,
#         data,
#         cmap='Reds',
#         shading='auto',
#         transform=ccrs.PlateCarree()
#     )
# ax.coastlines()
# plt.savefig('/home/548/cd3022/figures/TEST_SEABREEZE_SHIFT_0.png')
# plt.close()
