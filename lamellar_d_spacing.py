#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec  4 10:47:09 2020

@author: Chris
"""


import ovito as ov
import numpy as np
import pickle
import glob

if __name__ == '__main__':    
    
    files = glob.glob('*.pdb')
    
    dims = np.zeros(0)
    
    for i in files:
        
        pipeline = ov.io.import_file(i)
            
        data = pipeline.compute()
        
        dims = np.append(dims, data.cell.matrix[2,2])
        
    
    d = {'lp': dims.mean(),
         'lp error': dims.std()/np.sqrt(np.shape(dims))}
    
    
    pickle.dump(d, open('lamellar lattice parameter.p', 'wb'))
    