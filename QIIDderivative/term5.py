# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""
import numpy as np
import math
pi=math.pi

def term5(x,y,z,m):
    #t5_x
    cg = (-2 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
          10 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(5 * pi * m * x) - 
          2 * math.cos(pi * m * z) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * y) - 
          2 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) + 
          10 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(5 * pi * m * x) + 
          2 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * y) + 
          2 * pi * m * math.cos(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) + 
          10 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi * m * math.cos(5 * pi * m * x) - 
          2 * math.sin(pi * m * z) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * y) + 
          2 * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) - 
          10 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi * m * math.sin(5 * pi * m * x) + 
          2 * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * y))
    
    #t5_y
    cg0 = (-2 * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) - 
           2 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi * m * math.cos(5 * pi * m * x) - 
           10 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * y) + 
           2 * pi * m * math.cos(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
           2 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi * m * math.sin(5 * pi * m * x) + 
           10 * math.cos(pi * m * z) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * y) - 
           2 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) + 
           2 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(5 * pi * m * x) + 
           10 * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * y) + 
           2 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) + 
           2 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(5 * pi * m * x) - 
           10 * math.sin(pi * m * z) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * y))
    
    #t5_z
    cg1 = (-10 * pi * m * math.cos(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
           2 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(5 * pi * m * x) - 
           2 * math.sin(pi * m * z) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * y) + 
           10 * pi * m * math.cos(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) + 
           2 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(5 * pi * m * x) - 
           2 * math.sin(pi * m * z) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * y) + 
           10 * pi * m * math.sin(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           2 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi * m * math.sin(5 * pi * m * x) + 
           2 * math.cos(pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * y) - 
           10 * pi * m * math.sin(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) + 
           2 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi * m * math.cos(5 * pi * m * x) + 
           2 * math.cos(pi * m * z) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * y))
    
    #t5_xx
    cg2 = (-2 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           50 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           2 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) - 
           50 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           2 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
           50 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           2 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) - 
           50 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           2 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y))
    
    #t5_xy
    cg3 = (2 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) + 
           10 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           10 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
           10 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           10 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) + 
           10 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           10 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) + 
           2 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           10 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           10 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y))
    
    #t5_yy
    cg4 = (-2 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           2 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           50 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) - 
           2 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           50 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
           2 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           50 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) - 
           2 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) - 
           2 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           50 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y))
    
    #t5_yz
    cg5 = (10 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) + 
           2 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           10 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) + 
           10 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           2 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           10 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) - 
           10 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) - 
           2 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           10 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
           10 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) + 
           2 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           10 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y))
    
    #t5_zz
    cg6 = (-50 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           2 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           2 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
           50 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) - 
           2 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           2 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) - 
           50 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
           2 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           2 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) - 
           50 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) - 
           2 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           2 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y))

    #t5_xz
    cg7 = (10 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) + 
           10 * math.cos(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           2 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) - 
           10 * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) + 
           10 * math.cos(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           2 * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) + 
           10 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           10 * math.sin(pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           2 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) - 
           10 * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) - 
           10 * math.sin(pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           2 * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y))

    return np.array([cg,cg0,cg1,cg2,cg3,cg4,cg5,cg6,cg7])