import numpy as np

def merge_structured_arrays(*arrays):
    fields = [field for array in arrays for field in array.dtype.descr]
    rv     = np.empty(len(arrays[0]), dtype = fields)
    for array in arrays:
        for field in array.dtype.names:
            rv[field] = array[field]
    return rv

def count(data, sorted = False):
    if len(data) == 0:
        return 0
    if not sorted:
        data = np.sort(data)

    if data.ndim == 1:
        diffs = data[:-1] != data[1:]
    else:
        diffs = np.any(data[:-1] != data[1:], 1)

    return np.sum(diffs) + 1


    
def group_apply(groups, data, func, sorted = False, args = (), kwargs = {}):
    """Group rows of data based on values in groups and apply func to each group"""

    if not len(data):
        return []

    if not sorted:
        order  = np.argsort(groups)
        groups = groups[order]
        data   = data[order]
    
    if groups.ndim == 1:
        diffs = np.nonzero(groups[:-1] != groups[1:])[0] + 1
    else:
        diffs = np.nonzero(np.any(groups[:-1] != groups[1:], 1))[0] + 1

    rvs  = []
    i = 0
    for j in diffs:
        rvs.append(func(data[i:j], *args, **kwargs))
        i = j
    rvs.append(func(data[i:], *args, **kwargs))

    return rvs

def is_struct_array(array):
    return array.dtype.fields != None

def group_aggregate(groups, data, func, sorted = False, dtype = None, args = (), kwargs = {}):
    if not sorted:
        order  = np.argsort(groups)
        groups = groups[order]
        data   = data[order]
    
    keys   = np.array(group_apply(groups, groups, lambda x: x[0], sorted = True), dtype = groups.dtype)
    values = np.array(group_apply(groups, data, func, sorted = True, args = args, kwargs = kwargs), dtype = dtype)
    
    if not is_struct_array(keys) and not is_struct_array(values):
        return np.column_stack([keys, values])
    else:
        if not is_struct_array(keys):
            keys = keys.view(dtype = [("key", keys.dtype, keys.shape[1:])])
            
        if not is_struct_array(values):
            values = values.view(dtype = [("value", values.dtype, values.shape[1:])]).flatten()
        
        return merge_structured_arrays(keys, values)

def masked_corrcoef(data, mask):
    mask              = np.array(mask, dtype = np.bool)
    keep              = mask == False
    masked_data       = np.array(data)
    masked_data[mask] = 0

    # create a series of matrices in which element (i, j) is some function of data[i] masked by i and j
    keepf             = np.array(keep, dtype = np.float)
    masked_size       = np.tensordot(keepf[:, None, :], keepf[None, :, :], (2, 2))[:, 0, 0, :]
    
    masked_mean       = np.tensordot(masked_data[:, None, :], keep[None, :, :], (2, 2))[:, 0, 0, :]/masked_size
    masked_moment_xx  = np.tensordot(np.square(masked_data[:, None, :]), keep[None, :, :], (2, 2))[:, 0, 0, :]/(masked_size - 1)
    masked_moment_xy  = np.dot(masked_data, masked_data.T)/(masked_size - 1)

    masked_var        = masked_moment_xx - masked_size/(masked_size - 1)*np.square(masked_mean)
    masked_cov        = masked_moment_xy - masked_size/(masked_size - 1)*masked_mean*masked_mean.T

    masked_corr       = masked_cov/np.sqrt(masked_var*masked_var.T)
    return np.maximum(np.minimum(masked_corr, 1.), -1.)

def outlier_mask(arr, num_outliers = 1, axis = 0, random_tiebreak = True):
    if num_outliers == 0:
        return np.zeros_like(arr, dtype = np.bool)

    if axis == 1:
        arr = arr.T
    
    if random_tiebreak:
        permutation = np.array([np.random.permutation(arr.shape[1]) for i in range(arr.shape[0])])
        permuted    = np.array([x[p] for (x, p) in zip(arr, permutation)])
    else:
        permuted    = arr

    ap = np.argpartition(permuted, permuted.shape[1] - num_outliers)

    mask = np.zeros_like(arr, dtype = np.bool)
    for i in range(len(ap)):
        mask[i, ap[i, -num_outliers:]] = True

    if random_tiebreak:
        permuted_mask = mask
        mask          = np.empty_like(mask)
        for i in range(len(mask)):
            mask[i, permutation[i]] = permuted_mask[i]

    return mask if axis == 0 else mask.T

def invert_indices(ind):
    inv = np.empty_like(ind)
    inv[ind] = np.arange(len(ind))
    return inv

class RandomSeedContext(object):
    def __init__(self, seed):
        self._seed = seed
    def __enter__(self):
        self._state = np.random.get_state()
        np.random.seed(self._seed)
        return self
    def __exit__(self, type, value, traceback):
        np.random.set_state(self._state)
        

