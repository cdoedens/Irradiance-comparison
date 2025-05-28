import os
import xarray as xr
import numpy as np
from pathlib import Path
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
import concurrent.futures
import glob

def himawari_era5_diff(date):

    date_dt = datetime.strptime(date, "%d-%m-%Y")
    year = date_dt.strftime("%Y")
    month = date_dt.strftime("%m")
    day = date_dt.strftime("%d")
    
    ##### Himawari Data
    if date_dt <= datetime.strptime('2019-03-31', '%Y-%m-%d'):
        version = 'v1.0'
    else:
        version = 'v1.1'
    directory=Path(f'/g/data/rv74/satellite-products/arc/der/himawari-ahi/solar/p1s/{version}/{year}/{month}/{day}')
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
    era5_file = [f for d in era5_dir for f in d.glob(f"msdwswrf_era5_oper_sfc_{year}{month}*.nc")][0]
    era5 = xr.open_dataset(era5_file)
    era5 = era5.isel(time=slice((int(day)-1)*24,int(day)*24)) # probably a better way to get the day
    
    # Restrict ERA5 domain to Himawari domain
    era5_aus = era5.sel(
        latitude=slice(-10, -44.5),
        longitude=slice(156.26, 112)
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
    file_name = f'msdwswrf-era5-himawari_{year}{month}{day}.nc'
    file_path = f'/g/data/er8/users/cd3022/Irradiance-comparisons/era5-himawari/daily/'
        
    os.makedirs(file_path, exist_ok=True)
    diff.to_netcdf(f'{file_path}/{file_name}')
    return

# if __name__ == '__main__':


#     n_procs = os.cpu_count()
#     worker_pool = concurrent.futures.ProcessPoolExecutor(max_workers=10)


#     start = '01-01-2020'
#     end = '31-01-2020'
#     start_dt = datetime.strptime(start, "%d-%m-%Y")
#     end_dt = datetime.strptime(end, "%d-%m-%Y")

#     # Generate a list of dates
#     date_range = [start_dt + timedelta(days=i) for i in range((end_dt - start_dt).days + 1)]
    
#     futures = {}
#     # Loop over the dates
#     for date in date_range: 
#         date_s = date.strftime("%d-%m-%Y")
#         future = worker_pool.submit(himawari_era5_diff, date_s) 
#         futures[future] = f"The job for {date_s}"




###################### NOT TESTED YET ############################

def process_month(start_date, end_date):
    # Generate date range
    date_range = [start_date + timedelta(days=i) for i in range((end_date - start_date).days + 1)]

    # Create worker pool
    with concurrent.futures.ProcessPoolExecutor(max_workers=10) as executor:
        futures = {
            executor.submit(himawari_era5_diff, date.strftime("%d-%m-%Y")): date
            for date in date_range
        }

        # Wait for all tasks to complete
        for future in concurrent.futures.as_completed(futures):
            date = futures[future]
            try:
                future.result()
                print(f"Finished processing {date.strftime('%d-%m-%Y')}")
            except Exception as e:
                print(f"Error processing {date.strftime('%d-%m-%Y')}: {e}")

    # Load daily files and compute monthly mean
    pattern = os.path.join(
        '/g/data/er8/users/cd3022/Irradiance-comparisons/era5-himawari/daily/',
        f"{start_date.strftime('%Y-%m')}-??.nc"
    )
    file_list = sorted(glob.glob(pattern))

    if file_list:
        print(f"Found {len(file_list)} files, computing monthly mean...")

        datasets = [xr.open_dataset(f) for f in file_list]
        combined = xr.concat(datasets, dim="time")
        monthly_mean = combined.mean(dim="time")
        monthly_mean_dir = f'/g/data/er8/users/cd3022/Irradiance-comparisons/era5-himawari/monthly/'
        monthly_mean_path = f'{monthly_mean_dir}{start_date.strftime('%Y-%m')}.nc'
        os.makedirs(monthly_dir, exist_ok=True)
        monthly_mean.to_netcdf(monthly_mean_path)
        print(f"Saved monthly mean to {monthly_mean_path}")

        # Close and delete daily files
        for ds, f in zip(datasets, file_list):
            ds.close()
            os.remove(f)
            print(f"Deleted {f}")
    else:
        print("No files found to average.")

if __name__ == '__main__':
    start = '01-01-2021'
    end = '31-01-2021'
    
    start_dt = datetime.strptime(start, "%d-%m-%Y")
    end_dt = datetime.strptime(end, "%d-%m-%Y")
    
    process_month(start_dt, end_dt)







