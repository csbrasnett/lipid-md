# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""
import numpy as np
import math
pi=math.pi

def term3(x,y,z,m):
    #t3_x
    cg = (-6 * pi * m * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
          2 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(pi * m * x) - 
          2 * math.cos(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) - 
          6 * pi * m * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
          2 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) + 
          2 * math.cos(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) - 
          6 * pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(pi * m * z) - 
          2 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) + 
          2 * math.sin(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y))
    
    #t3_y
    cg0 = (-2 * math.cos(3 * pi * m * x) * pi * m * math.sin(pi * m * y) * math.cos(pi * m * z) - 
           6 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) - 
           2 * math.cos(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) + 
           2 * pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(pi * m * z) - 
           6 * pi * m * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * math.sin(pi * m * x) + 
           2 * math.cos(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) + 
           2 * pi * m * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           6 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(pi * m * x) - 
           2 * math.sin(3 * pi * m * z) * math.cos(pi * m * x) * pi * m * math.cos(pi * m * y))
    
    #t3_z
    cg1 = (-2 * pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(pi * m * z) - 
           2 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) - 
           6 * math.sin(3 * pi * m * z) * math.cos(pi * m * x) * pi * m * math.cos(pi * m * y) + 
           2 * math.cos(3 * pi * m * x) * pi * m * math.sin(pi * m * y) * math.cos(pi * m * z) + 
           2 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(pi * m * x) - 
           6 * math.sin(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) - 
           2 * pi * m * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) + 
           2 * pi * m * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * math.sin(pi * m * x) - 
           6 * math.cos(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y))
    
    #t3_xx
    cg2 = (-18 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           2 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) - 
           2 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) - 
           18 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           2 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           2 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) + 
           18 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(pi * m * z) + 
           2 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) + 
           2 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y))
    
    #t3_xy
    cg3 = (6 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(pi * m * z) + 
           6 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) + 
           2 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) - 
           6 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(pi * m * z) - 
           6 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * math.cos(pi * m * x) + 
           2 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) + 
           6 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           6 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) + 
           2 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y))
    
    #t3_yy
    cg4 = (-2 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           18 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) - 
           2 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           18 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           2 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) + 
           2 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(pi * m * z) + 
           18 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) + 
           2 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y))
    
    #t3_yz
    cg5 = (2 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) + 
           6 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * math.cos(pi * m * x) + 
           6 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) + 
           2 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           6 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           6 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y) + 
           2 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(pi * m * z) + 
           6 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           6 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y))
    
    #t3_zz
    cg6 = (-2 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) - 
           2 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) - 
           18 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(pi * m * z) - 
           2 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) - 
           18 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) + 
           2 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(pi * m * z) + 
           2 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) + 
           18 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y))

    #t3_xz
    cg7 = (6 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(pi * m * z) + 
           2 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) + 
           6 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y) - 
           6 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(pi * m * z) + 
           2 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) - 
           6 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) - 
           6 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(pi * m * z) + 
           2 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * math.cos(pi * m * x) + 
           6 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y))

    return np.array([cg,cg0,cg1,cg2,cg3,cg4,cg5,cg6,cg7])