import pickle
from mesas.sas.model import Model
import matplotlib.pyplot as plt
import numpy as np




filename = '/Users/ciaran/Documents/Research/TVTTD/MESAS/paper/python/runs/yr10_to_yr13_1D-2pt_scan_5_all_components_13_4_model.pickle'

with open(filename, 'rb') as f:
    model = pickle.load(f)
#model.options = {'n_substeps':20}# Didn't work
model.run()

PQ = model.result['PQ'][:,:,0]
print(PQ.min())
print(PQ.max())
print(PQ[430,501])


