import os
import xarray as xr
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import logger
import xesmf as xe

from dask.distributed import as_completed
import climtas.nci

LOG = logger.get_logger(__name__)

def himawari_era5_diff(date):

    date_dt = datetime.strptime(date, "%m-%Y")
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
        # engine='netcdf4',
    )
    LOG.info('OPEN HIMAWARI')
    # himawari = himawari.chunk({'time':100})
    LOG.info('CHUNK HIMAWARI')
    LOG.info(f'HIMAWARI SHAPE: {himawari.surface_global_irradiance.shape}')

    ##### ERA5 Data
    era5_dir = [
        Path(f"/g/data/rt52/era5/single-levels/reanalysis/msdwswrf/{year}"),
    ]
    era5_file = [f for d in era5_dir for f in d.glob(f"msdwswrf_era5_oper_sfc_{year}{month}*.nc")][0]
    era5 = xr.open_dataset(era5_file, chunks={'time':100})
    LOG.info('OPEN ERA5')
    
    # Restrict ERA5 domain to Himawari domain
    era5_aus = era5.sel(
        latitude=slice(-10, -44.5),
        longitude=slice(112, 156.26)
    )
    LOG.info(f'ERA5 SLICED SHAPE: {era5_aus.msdwswrf.shape}')
    
    # Get monthly mean
    himawari_clim = himawari.surface_global_irradiance.mean(dim='time')
    era5_clim = era5_aus.msdwswrf.mean(dim='time')
    LOG.info(f'HIMARI_CLIM SHAPE: {himawari_clim.shape}')
    LOG.info(f'ERA5_CLIM SHAPE: {era5_clim.shape}')
    
    # # Regrid ERA5 for comparison
    # era5_clim_interp = era5_clim.interp(
    #     latitude=himawari_clim.latitude,
    #     longitude=himawari_clim.longitude,
    #     method='linear'
    # )
    
    # Regrid himawari to coarse ERA5
    regridder = xe.Regridder(himawari_clim, era5_clim, method="conservative")
    himawari_clim_regridded = regridder(himawari_clim)
    
    diff = era5_clim - himawari_clim_regridded

    # Save results
    file_name = f'msdwswrf-era5-himawari_{year}{month}.nc'
    file_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/era5-himawari/monthly/')
    file_path.mkdir(parents=True, exist_ok=True)
    diff.to_netcdf(file_path / file_name)
    return f"Finished processing {year}-{month}"

if __name__ == '__main__':
    climtas.nci.GadiClient()
    date = sys.argv[1]
    futures = {}
    future = client.submit(himawari_era5_diff, date)
    futures[future] = f"The job for {date}"
    for future in as_completed(futures):
        _ = future.result()