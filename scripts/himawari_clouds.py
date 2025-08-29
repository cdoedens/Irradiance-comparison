import xarray as xr
import numpy as np
import scipy.stats as stats

ds = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")

ds.sel(time='2020-01-01T00:00').ct.plot(cmap='magma')

def mode_func(x, axis):
    mode_val, _ = stats.mode(x, axis=axis, nan_policy="omit")
    return mode_val

ds.sel(time='2020-01-01')..resample(time='1D').reduce(mode_func)