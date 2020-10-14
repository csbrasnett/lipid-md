# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""
import numpy as np
import math
pi=math.pi

def term2(x,y,z,m):
    #t2_x
    cg = (-6 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
          6 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
          6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) + 
          6 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y))
    
    #t2_y
    cg0 = (-6 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) + 
           6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) - 
           6 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
           6 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z))
    
    #t2_z
    cg1 = (-6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) + 
           6 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) + 
           6 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           6 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z))
    
    #t2_xx
    cg2 = (-6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           6 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y))
    
    #t2_xy
    cg3 = (6 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) - 
           6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
           6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z))
    
    #t2_yy
    cg4 = (-6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           6 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y))
    
    #t2_yz
    cg5 = (6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
           6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           6 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) - 
           6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x))
    
    #t2_zz
    cg6 = (-6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           6 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y))

    #t2_xz
    cg7 = (6 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           6 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) + 
           6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           6 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z))

    return np.array([cg,cg0,cg1,cg2,cg3,cg4,cg5,cg6,cg7])