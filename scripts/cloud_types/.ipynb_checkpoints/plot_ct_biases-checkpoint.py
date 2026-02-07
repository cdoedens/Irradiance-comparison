import xarray as xr
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pathlib import Path

ct_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/cloud_types')
low_cloud_files = [f for f in ct_path.glob('low*.nc')]
ds_low = xr.open_mfdataset(low_cloud_files, combine='nested', concat_dim='date')

low_ct_mean = ds_low.low_ghi.mean(dim='date')

# PLOTTING
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
ax.pcolormesh(
    low_ct_mean.lon,
    low_ct_mean.lat,
    low_ct_mean,
)
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/cloud_types//TEST.png')