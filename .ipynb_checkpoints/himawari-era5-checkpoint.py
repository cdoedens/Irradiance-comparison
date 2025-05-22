import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cartopy.crs as ccrs
from pathlib import Path

##### ERA5 Data
# era5_dirs = [
#     Path("/g/data/rt52/era5/single-levels/reanalysis/msdwswrf"),
# ]
# era5_files = sorted([f for d in era5_dirs for f in d.glob("*.nc")])
# era5 = xr.open_mfdataset(
#     era5_files,
#     combine='by_coords',   # Use coordinate-based merging
#     data_vars='minimal',   # Avoid trying to align data_vars across files
#     coords='minimal',      # Avoid "coords='different'" error
#     compat='override',     # Avoid unnecessary compatibility checks
# )


# test with smaller dataset first
era5 = xr.open_dataset('/g/data/rt52/era5/single-levels/reanalysis/msdwswrf/2020/msdwswrf_era5_oper_sfc_20200301-20200331.nc')

himawari = xr.open_dataset('/g/data/rv74/satellite-products/arc/der/himawari-ahi/solar/p1s/v1.1/2020/01/01/IDE00326.201912312210.nc')

lat_max = himawari.latitude.max().item()
lat_min = himawari.latitude.min().item()
lon_max = himawari.longitude.max().item()
lon_min = himawari.longitude.min().item()

era5_aus = era5.sel(
    latitude=slice(lat_min, lat_max),
    longitude=slice(lon_min, lon_max)
)

himawari_clim = himawari.surface_global_irradiance.mean(dim='time')
era5_clim = era5_aus.msdwswrf.mean(dim='time')

era5_clim_interp = era5_clim.interp(
    latitude=himawari_clim.latitude,
    longitude=himawari_clim.longitude,
    method='linear'
)

diff = era5_clim_interp - himawari_clim

# Create a figure
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

mesh=ax.pcolormesh(himawari.longitude.values, himawari.latitude.values, diff, cmap='BrBG', transform=ccrs.PlateCarree())

ax.set_title('ERA5 - Himawari')
ax.coastlines()
cbar = plt.colorbar(mesh,ax=ax,shrink=0.5)
cbar.ax.tick_params(labelsize=5)  # Set the font size for the colorbar ticks
cbar.set_label('Irradiance bias (W/m2)', fontsize=5) 
 
plt.tight_layout()

plt.savefig('/home/548/cd3022/figures/himawari-era5-diff.png')


