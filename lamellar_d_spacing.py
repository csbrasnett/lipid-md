#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 10:47:09 2020

@author: Chris
"""


import ovito as ov
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import optimize
import pickle
from itertools import product
from multiprocessing import get_context


if __name__ == '__main__':    
    
    files = glob.glob('*.pdb')
    
    dims = np.zeros(0)
    
    for i in files:
        
        pipeline = ov.io.import_file(i)
            
        data = pipeline.compute()
        
        dims = np.append(dims, data.cell.matrix[2,2])
        
    
    print('d = %.3f ± %.3f' %(dims.mean(), dims.std()/np.sqrt(np.shape(dims))))
        
