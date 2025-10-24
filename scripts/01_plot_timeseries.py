import xarray as xr
from pathlib import Path
import matplotlib.pyplot as plt

file_path = Path('/g/data/er8/users/cd3022/Irradiance-comparisons/error_timeseries')
files = [f for f in file_path.glob('*.nc')]
ds = xr.open_mfdataset(files)

plt.figure(figsize=(16,10))
plt.plot(ds.time, ds.rmse)
plt.savefig('/home/548/cd3022/figures/Irradiance_comparisons/rmse_timeseries.png')
plt.show()