#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""

import numpy as np
from .QIIDderivative import derivative

def nominator(F_x, F_y, F_z, F_xx, F_xy, F_yy, F_yz, F_zz, F_xz):
    m = np.array([[F_xx, F_xy, F_xz, F_x],
                 [F_xy, F_yy, F_yz, F_y],
                 [F_xz, F_yz, F_zz, F_z],
                 [F_x, F_y, F_z, 0]])
        
    d = np.linalg.det(m)
    
    return d

def denominator(F_x,F_y, F_z):
    
    g = np.array([F_x,F_y,F_z])
        
    mag_g = np.linalg.norm(g)
      
    return mag_g**4

def main(x, y, z, lamb):
    vals = derivative(x, y, z, lamb)
    
    n = nominator(vals[0],vals[1],vals[2],vals[3],vals[4],vals[5],vals[6],vals[7],vals[8])
    d = denominator(vals[0],vals[1],vals[2])
    K = -(n/d)
    
    return K