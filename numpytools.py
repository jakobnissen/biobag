
"""Bioinformatics toolbox

This is a loose collection of tools that require the numpy package.

Author: Jakob Nybo Nissen, DTU Bioinformatics
"""



import numpy as _np



def zscore(array, axis=None, inplace=False):
    "Calculates zscore for an array. A cheap copy of scipy.stats.zscore."
    
    if axis is not None and axis >= array.ndim:
        raise _np.AxisError('array only has {} axes'.format(array.ndim))
        
    if inplace and array.dtype not in (_np.float, _np.float16, _np.float32, _np.float64, _np.float128):
            raise TypeError('Cannot convert a non-float array to zscores')
        
    mean = array.mean(axis=axis)
    std = array.std(axis=axis)
    
    if axis is None:
        if std == 0:
            std = 1 # prevent divide by zero
            
    else:
        std[std == 0.0] = 1 # prevent divide by zero
        shape = tuple(dim if ax != axis else 1 for ax, dim in enumerate(array.shape))
        mean.shape, std.shape = shape, shape

    if inplace:
        array -= mean
        array /= std
        return None
    else:
        return (array - mean) / std

