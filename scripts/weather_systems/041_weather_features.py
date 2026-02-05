import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import sys, os
import datetime

from read_datasets import read_dataset
import logger

LOG = logger.get_logger(__name__)
# qsub -I -q normal -P er8 -l walltime=4:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+scratch/er8

from dask.distributed import Client


if __name__ == '__main__':

    client = Client(
        n_workers=12,
        threads_per_worker=1
    )

    

    dataset = sys.argv[1]
    wf = sys.argv[2]
    year = sys.argv[3]
    
    if wf == 'frovo_id':
        lvl = 850
    else:
        lvl = 'sfc'

    LOG.info(f'Starting analysis for dataset {dataset}, weather feature {wf}, year {year}')
    # GHI dataset
    ds_list = []
    for month in range(1, 13):
        # LOAD DATASETS
        ds_month = read_dataset(
                dataset=dataset,
                resolution='hourly',
                date=f'{year}-{month:02d}'
            )
        ds_list.append(ds_month)
    ds = xr.concat(ds_list, dim='time')
    LOG.info(f'{dataset} data opened')
            
    # Weather Feature dataset
    def preprocess(_ds):
        return _ds.sel(
            latitude=slice(ds.lat.max(), ds.lat.min()),
            longitude=slice(ds.lon.min(), ds.lon.max())
        ).rename(
            {'latitude':'lat',
            'longitude':'lon'}
        )
    

    wf_path = Path(f'/scratch/er8/cd3022/weather-features/cyc-anticyc-front-cao/')
    wf_files = [f for f in wf_path.glob(f'*{year}??.{lvl}.nc')]
    wf_ds = xr.open_mfdataset(wf_files) #, preprocess=preprocess)
    LOG.info(f'{wf} data opened')

    # remove the level dimension if looking at frontal volume
    if wf == 'frovo_id':
        wf_ds = wf_ds.isel(lev=0)

    # PREPROCESS DATASETS TO MATCH TIMES/SHAPES
    # interp wf mask from era5 to barra
    wf_ds = wf_ds.interp(
        lat=ds.lat,
        lon=ds.lon,
        method='nearest'
    )

    if dataset == 'himawari':
        # fill missing overnight time steps
        full_time = pd.date_range(
            start=ds.time.min().item(),
            end=ds.time.max().item(),
            freq="60min"
        )
        ds = ds.reindex(time=full_time, fill_value=0)
        
        # remove timesteps from the end of wf_ds to match himawari,
        # missing timesteps are overnight so will not impact final data
        wf_ds = wf_ds.isel(time=slice(0, -4))

    # adjust hourly times on the 30min to match 3hr times (on the hour) from wf
    ds_shifted = ds.assign_coords(time=ds.time - pd.Timedelta('30min'))
    ds_times = ds_shifted.sel(time=wf_ds.time)
    LOG.info('preprocessing complete')


    # APPLY MASK
    da_wf = xr.where(wf_ds[wf] != 0, ds_times.ghi, np.nan)
    LOG.info('mask applied')

    # TIME MEAN
    da_wf_mean = da_wf.mean(dim='time')
    LOG.info('annual mean calculated')

    # remask himawari region
    if dataset == 'himawari':
        da_wf_mean = xr.where(ds.isel(time=5).ghi.isnull(), np.nan, da_wf_mean)

    # add metadata
    ds_wf_mean = da_wf_mean.to_dataset(name='annual_mean_ghi')
    ds_wf_mean = ds_wf_mean.assign_coords({'year':year})
    ds_wf_mean.attrs['date_generated'] = datetime.date.today().strftime('%D')
    ds_wf_mean.attrs['source_script'] = 'data produced by the script "041_weather_features.py"'

    # SAVE DATA
    file_path = Path(f'/g/data/er8/users/cd3022/Irradiance-comparisons/weather-features/{wf}')
    os.makedirs(file_path, exist_ok=True)
    ds_wf_mean.to_netcdf(f'{file_path}/{dataset}-{year}.nc')
    LOG.info('data saved, job complete!')
        
        

