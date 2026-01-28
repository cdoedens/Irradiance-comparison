import xarray as xr
from pathlib import Path

date = '202001'

sb_path = Path('/g/data/ng72/ab4502/sea_breeze_detection/barra_c_smooth_s2/filters/')
sb_files = [f for f in sb_path.glob(f'*_F_{date}*.zarr')]
print(len(sb_files))

ds_sb = xr.open_mfdataset(
    sb_files,
    engine="zarr",
    combine="by_coords",
    parallel=True,
)