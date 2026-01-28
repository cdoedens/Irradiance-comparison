import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

date = '202001'

sb_path = Path('/g/data/ng72/ab4502/sea_breeze_detection/barra_c_smooth_s2/filters/')
sb_files = [f for f in sb_path.glob(f'*_F_{date}*.zarr')]

ds_sb = xr.open_mfdataset(
    sb_files,
    engine="zarr",
    combine="by_coords",
    parallel=True,
)

data = ds_sb.mask.isel(time=slice(0,5)).sum(dim='time')

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
mesh = ax.pcolormesh(
        ds_sb.lon,
        ds_sb.lat,
        data,
        cmap='Reds',
        shading='auto',
        transform=ccrs.PlateCarree()
    )
ax.coastlines()
plt.savefig('/home/548/cd3022/figures/TEST_SEABREAZE.png')
plt.close()
