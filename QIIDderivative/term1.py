# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""
import numpy as np
import math
pi=math.pi

def term1(x,y,z,m):
    #t1_x
    cg = (-3 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
          3 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
          3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) + 
          3 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y))
    
    #t1_y
    cg0 = (-3 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) + 
           3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) - 
           3 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
           3 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z))
    
    #t1_z
    cg1 = (-3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) + 
           3 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) + 
           3 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           3 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z))
    
    #t1_xx
    cg2 = (-3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           3 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y))
    
    #t1_xy
    cg3 = (3 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) - 
           3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
           3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z))
    
    #t1_yy
    cg4 = (-3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           3 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y))
    
    #t1_yz
    cg5 = (3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
           3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           3 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) - 
           3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x))
    
    #t1_zz
    cg6 = (-3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           3 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y))
    
    #t1_xz
    cg7 = (3 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           3 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) + 
           3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           3 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z))

    return np.array([cg,cg0,cg1,cg2,cg3,cg4,cg5,cg6,cg7])