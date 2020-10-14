# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""
import numpy as np
import math
pi=math.pi

def term6(x,y,z,m):
    #t6_x
    cg = (-18 * pi * m * math.sin(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z) -
          18 * pi * m * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z) + 
          18 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi * m * math.cos(3 * pi * m * x) + 
          18 * math.cos(3 * pi * m * z) * pi * m * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y))
    
    #t6_y
    cg0 = (-18 * math.cos(3 * pi * m * z) * pi * m * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y) + 
           18 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi * m * math.cos(3 * pi * m * x) - 
           18 * pi * m * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z) + 
           18 * pi * m * math.sin(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z))
    
    #t6_z
    cg1 = (-18 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi * m * math.cos(3 * pi * m * x) + 
           18 * math.cos(3 * pi * m * z) * pi * m * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y) + 
           18 * pi * m * math.sin(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z) - 
           18 * pi * m * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z))
    
    #t6_xx
    cg2 = (-54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z) - 
           54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z) - 
           54 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) - 
           54 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y))
    
    #t6_xy
    cg3 = (54 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y) - 
           54 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) - 
           54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z) + 
           54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z))
    
    #t6_yy
    cg4 = (-54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z) - 
           54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z) - 
           54 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) - 
           54 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y))
    
    #t6_yz
    cg5 = (54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z) + 
           54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z) - 
           54 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y) - 
           54 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x))
    
    #t6_zz
    cg6 = (-54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z) - 
           54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z) - 
           54 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) - 
           54 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y))

    #t6_xz
    cg7 = (54 * math.cos(3 * pi * m * y) * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) - 
           54 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(3 * pi * m * y) + 
           54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(3 * pi * m * y) * math.cos(3 * pi * m * z) - 
           54 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(3 * pi * m * y) * math.sin(3 * pi * m * z))

    return np.array([cg,cg0,cg1,cg2,cg3,cg4,cg5,cg6,cg7])