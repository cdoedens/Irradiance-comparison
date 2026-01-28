import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

file_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/error_timeseries')
files = [f for f in file_path.glob('*.nc')]
ds = xr.open_mfdataset(files)

# remove outliers
ds['rmse'] = xr.where(ds.rmse > 100, np.nan, ds.rmse)
ds['mbe'] = xr.where(abs(ds.mbe) > 50, np.nan, ds.mbe)

ds_year = ds.groupby('time.dayofyear').mean('time')
ds_year = ds_year.assign_coords(
    date=("dayofyear", pd.date_range("2000-01-01", periods=ds_year.dayofyear.size))
)

# Plot
fig, ax = plt.subplots(ncols=1, nrows=2, figsize=(16, 10))

ax[0].plot(ds_year.date, ds_year.rmse)
ax[1].plot(ds_year.date, ds_year.mbe)

ax[0].set_title('RMSE')
ax[1].set_title('MBE')

for i in [0,1]:
    ax[i].xaxis.set_major_locator(mdates.MonthLocator())
    ax[i].xaxis.set_major_formatter(mdates.DateFormatter("%b"))
    ax[i].set_xlabel("")

plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/rmse-mbe.png')
plt.close()


