import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Error data
err_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/error_timeseries')
err_files = [f for f in err_path.glob('*.nc')]
ds_err = xr.open_mfdataset(err_files)
# remove outliers
ds_err['rmse'] = xr.where(ds_err.rmse > 100, np.nan, ds_err.rmse)
ds_err['mbe'] = xr.where(abs(ds_err.mbe) > 50, np.nan, ds_err.mbe)

# Synoptic Weather Type data
swt_path = Path('/home/548/cd3022/repos/Australian-synoptic-weather-types/')
ds_swt = xr.open_dataset(swt_path / 'SWT_fields/SWT_climatology.nc')
# line up times with other dataset
ds_swt = ds_swt.assign_coords(time=ds_swt.time.dt.floor("D"))
ds_swt = ds_swt.sel(time=ds_err.time)

# Combine into one dataset
ds = xr.Dataset(
    {
        'rmse': ds_err.rmse,
        'mbe': ds_err.mbe,
        'SWTs': ds_swt.SWTs,
        'assigned_SWT': ds_swt.assigned_SWT
    }
)

# Calculate errors for individual SWTs
rmse_data = {}
mbe_data = {}

hot_months = [11, 12, 1, 2, 3]
rmse_hot_data = {}
mbe_hot_data = {}
for swt in ds.SWTs:
    data = ds.where(ds.assigned_SWT == swt.item(), drop=True)
    hot_data = ds.sel(time=ds.time.dt.month.isin(hot_months)).where(ds.assigned_SWT == swt.item(), drop=True)

    if hot_data.time.size < 10:
        continue

    rmse_data[swt.item()] = data.rmse.mean().compute().item()
    mbe_data[swt.item()] = data.mbe.mean().compute().item()

    rmse_hot_data[swt.item()] = hot_data.rmse.mean().compute().item()
    mbe_hot_data[swt.item()] = hot_data.mbe.mean().compute().item()


# Plot
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8, 5))

ax[0].bar(rmse_data.keys(), rmse_data.values())
ax[1].bar(mbe_data.keys(), mbe_data.values())

ax[0].set_ylabel('RMSE (Wm-2)')
ax[1].set_ylabel('MBE (Wm-2)')
ax[1].set_xlabel('SWT')

fig.suptitle('Errors Associated with Synoptic Weather Types')

for i in [0,1]:
    ax[i].tick_params(axis='x', rotation=90)
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/SWT-errors.png')
plt.close()


# Whole year vs hot months
SWTs = list(rmse_data.keys())
mbe_vals = list(mbe_data.values())
mbe_hot_vals = list(mbe_hot_data.values())

x = np.arange(len(SWTs))   # positions on x-axis
width = 0.35               # width of each bar

fig, ax = plt.subplots(figsize=(8, 4))

ax.bar(x - width/2, mbe_vals, width, label="All Data")
ax.bar(x + width/2, mbe_hot_vals, width, label="Nov - Mar")

ax.set_xticks(x)
ax.set_xticklabels(SWTs)
ax.tick_params(axis='x', rotation=90)
ax.set_ylabel("MBE (Wm-2)")
ax.set_title('Synoptic Weather Types MBE')
ax.legend()

plt.tight_layout()
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/SWT-mbe-hot.png')
plt.close()