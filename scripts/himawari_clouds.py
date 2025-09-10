import xarray as xr
import numpy as np
import scipy.stats as stats

ds = xr.open_zarr("/g/data/su28/himawari-ahi/cloud/ct/aus_regional_domain/S_NWC_CT_HIMA08_HIMA-N-NR.zarr/")

ds = ds.sel(time='2020-01')

ds = ds.chunk({"time": 103, "lat": 100, "lon": 100})

def fast_mode(arr, axis=-1):
    arr = np.moveaxis(arr, axis, -1)
    out_shape = arr.shape[:-1]
    arr_flat = arr.reshape(-1, arr.shape[-1])

    modes = np.full(arr_flat.shape[0], np.nan)
    for i, row in enumerate(arr_flat):
        row = row[~np.isnan(row)]
        if row.size > 0:
            counts = np.bincount(row.astype(int))
            modes[i] = np.argmax(counts)
    return modes.reshape(out_shape)



ds.sel(time='2020-01-01').resample(time='1D').reduce(fast_mode)