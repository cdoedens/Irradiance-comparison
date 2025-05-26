import xarray as xr
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta

start = '1-1-2016'
months_per_batch = 1
num_batches = 12 * 9
start_dt = datetime.strptime(start, "%d-%m-%Y")
for x in range(num_batches):
    date_dt = start_dt + relativedelta(months = months_per_batch * x)
    year = date_dt.strftime("%Y")
    month = date_dt.strftime("%m")
    
    
    
    ##### Himawari Data
    if date_dt <= datetime.strptime('2019-03-31', '%Y-%m-%d'):
        version = 'v1.0'
    else:
        version = 'v1.1'
    directory=Path(f'/g/data/rv74/satellite-products/arc/der/himawari-ahi/solar/p1s/{version}/{year}/{month}')
    files = sorted(str(p) for p in directory.rglob("*.nc"))
    def _preprocess(ds):
        return ds.drop_vars(set(ds.data_vars) - {'surface_global_irradiance'})
    
    himawari = xr.open_mfdataset(
        files,
        combine='by_coords',
        preprocess=_preprocess,
        engine='h5netcdf',
    )
    

    ##### ERA5 Data
    era5_dir = [
        Path(f"/g/data/rt52/era5/single-levels/reanalysis/msdwswrf/{year}"),
    ]
    era5_file = [f for d in era5_dirs for f in d.glob(f"msdwswrf_era5_oper_sfc_{year}{month}*.nc")][0]
    era5 = xr.open_dataset(era5_file)
    
    # Restrict ERA5 domain to Himawari domain
    lat_max = himawari.latitude.max().item()
    lat_min = himawari.latitude.min().item()
    lon_max = himawari.longitude.max().item()
    lon_min = himawari.longitude.min().item()
    era5_aus = era5.sel(
        latitude=slice(lat_max, lat_min),
        longitude=slice(lon_min, lon_max)
    )
    
    # Get monthly mean
    himawari_clim = himawari.surface_global_irradiance.mean(dim='time')
    era5_clim = era5_aus.msdwswrf.mean(dim='time')
    
    # Regrid ERA5 for comparison
    era5_clim_interp = era5_clim.interp(
        latitude=himawari_clim.latitude,
        longitude=himawari_clim.longitude,
        method='linear'
    )
    
    # Find the difference
    diff = era5_clim_interp - himawari_clim
    
    # Save results
    file_name = f'msdwswrf-era5-himawari_{year}{month}.nc'
    file_path = f'/g/data/er8/users/cd3022/Irradiance-comparisons/era5-himawari/'
        
    os.makedirs(file_path, exist_ok=True)
    diff.to_netcdf(f'{file_path}/{file_name}')

