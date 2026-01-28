import xarray as xr
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

file_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/error_timeseries')
files = [f for f in file_path.glob('*.nc')]
ds = xr.open_mfdataset(files)

# remove outliers
ds['rmse'] = xr.where(ds.rmse > 100, np.nan, ds.rmse)
ds['mbe'] = xr.where(abs(ds.mbe) > 50, np.nan, ds.mbe)

ds_year = ds.groupby('time.dayofyear').mean('time')

# Plot
plt.figure(figsize=(16,10))
plt.plot(ds_year.dayofyear, ds_year.rmse)
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/rmse_timeseries.png')
plt.close()

plt.figure(figsize=(16,10))
plt.plot(ds_year.dayofyear, ds_year.mbe)
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/mbe_timeseries.png')
plt.close()




# Plot
plt.figure(figsize=(16,10))
plt.plot(ds.time, ds.rmse)
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/rmse_timeseries.png')
plt.close()

plt.figure(figsize=(16,10))
plt.plot(ds.time, ds.mbe)
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/mbe_timeseries.png')
plt.close()