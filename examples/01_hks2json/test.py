#!/usr/bin/env python
import json
import numpy as np
from numpy import testing
from scipy.io import loadmat as load_mat
import h5py
import sys
from functools import partial

complex_keys = "H", "S"

def load_json(name):
    with open(name, 'r') as f:
        data = json.load(f)
    for k in data:
        if isinstance(data[k], (list, tuple)):
            data[k] = np.array(data[k])
        if k in complex_keys:
            data[k] = data[k][..., 0] + data[k][..., 1] * 1.j
    return data

def load_h5(name):
    data = {}
    with h5py.File(name, 'r') as f:
        for k in f:
            data[k] = np.array(f[k])
            if k in complex_keys:
                data[k] = data[k][..., 0] + data[k][..., 1] * 1.j
    return data

reference = load_json("reference_output/lead-chain.json")
fname = sys.argv[1]
ext = fname.split(".")[-1]
load_mat = partial(load_mat, squeeze_me=True)
actual = {"json": load_json, "h5": load_h5, "mat": load_mat}[ext](fname)

for k in reference:
    testing.assert_equal(reference[k], actual[k])

