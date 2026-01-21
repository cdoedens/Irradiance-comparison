import xarray as xr
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

from read_datasets import read_dataset

# qsub -I -q normal -P er8 -l walltime=4:00:00,ncpus=24,mem=120GB,jobfs=100MB,storage=gdata/xp65+gdata/er8+gdata/ob53+gdata/rt52+gdata/rv74+scratch/er8

year=2020
month=1

# LOAD DATASETS
barra = read_dataset(
        dataset='barra-r2',
        resolution='hourly',
        date=f'{year}-{month:02d}'
    )

himawari = read_dataset(
    dataset='himawari',
    resolution='hourly',
    date=f'{year}-{month:02d}'
)

def preprocess(ds):
    return ds.sel(
        latitude=slice(barra.lat.max(), barra.lat.min()),
        longitude=slice(barra.lon.min(), barra.lon.max())
    ).rename(
        {'latitude':'lat',
        'longitude':'lon'}
    )

wf_path = Path(f'/scratch/er8/cd3022/weather-features/cyc-anticyc-front-cao/')
cyc_files = [f for f in wf_path.glob(f'*{year}{month:02d}.sfc.cycmask.nc')]
cyc_ds = xr.open_mfdataset(cyc_files, preprocess=preprocess)

# PREPROCESS DATASETS TO MATCH TIMES/SHAPES
# interp him to barra
himawari = himawari.interp(
    lat=barra.lat,
    lon=barra.lon,
    method='nearest'
)

# interp wf mask from era5 to barra
cyc_ds = cyc_ds.interp(
    lat=barra.lat,
    lon=barra.lon,
    method='nearest'
)

# fill missing time steps
full_time = pd.date_range(
    start=himawari.time.min().item(),
    end=himawari.time.max().item(),
    freq="60min"
)
himawari = himawari.reindex(time=full_time, fill_value=0)

#

# adjust hourly times on the 30min to match 3hr times (on the hour) from wf
barra_shifted = barra.assign_coords(time=barra.time - pd.Timedelta('30min'))
barra_times = barra_shifted.sel(time=cyc_ds.time)

himawari_shifted = himawari.assign_coords(time=himawari.time - pd.Timedelta('30min'))
himawari_times = himawari_shifted.sel(time=cyc_ds.isel(time=slice(0, -4)).time)


# APPLY MASK
barra_cyc_ghi = xr.where(cyc_ds.cycmask == 1, barra_times.ghi, np.nan)
himawari_cyc_ghi = xr.where(cyc_ds.isel(time=slice(0, -4)).cycmask == 1, himawari_times.ghi, np.nan)

# TIME MEAN
barra_cyc_mean = barra_cyc_ghi.mean(dim='time')
himawari_cyc_mean = himawari_cyc_ghi.mean(dim='time')

diff = barra_cyc_mean - himawari_cyc_mean
diff = xr.where(himawari.isel(time=4).ghi.isnull(), np.nan, diff)

# PLOT
fig, ax = plt.subplots(figsize=(16,10), subplot_kw={'projection': ccrs.PlateCarree()})
mesh = ax.pcolormesh(
    diff.lon,
    diff.lat,
    diff,
    cmap='RdBu_r',
    vmin=-250, vmax=250,
    shading='auto',
    transform=ccrs.PlateCarree()
)
ax.coastlines()
ax.set_title('Bias during Cyclones')
cbar = fig.colorbar(mesh, ax=ax, orientation='vertical', fraction=0.02, pad=0.04)
cbar.set_label("Surface Downward SW bias (W/m²)")
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/weather-features/cyc-2020-01.png')
plt.close()