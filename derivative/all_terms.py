# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""

from .term1 import term1
from .term2 import term2
from .term3 import term3
from .term4 import term4
from .term5 import term5
from .term6 import term6
from .term7 import term7

import numpy as np

def derivative(x,y,z,lamb):
    t1= term1(x,y,z,lamb)
    t2= 0.1407*term2(x,y,z,lamb)
    t3= 0.020*term3(x,y,z,lamb)
    t4= -0.0138*term4(x,y,z,lamb)
    t5= -0.0028*term5(x,y,z,lamb)
    t6= -0.0021*term6(x,y,z,lamb)
    t7= 0.001*term7(x,y,z,lamb)

    F_x = np.array([t1[0],t2[0],t3[0],t4[0],t5[0],t6[0],t7[0]])
    F_y = np.array([t1[1],t2[1],t3[1],t4[1],t5[1],t6[1],t7[1]])
    F_z = np.array([t1[2],t2[2],t3[2],t4[2],t5[2],t6[2],t7[2]])
    F_xx = np.array([t1[3],t2[3],t3[3],t4[3],t5[3],t6[3],t7[3]])
    F_xy = np.array([t1[4],t2[4],t3[4],t4[4],t5[4],t6[4],t7[4]])
    F_yy = np.array([t1[5],t2[5],t3[5],t4[5],t5[5],t6[5],t7[5]])
    F_yz = np.array([t1[6],t2[6],t3[6],t4[6],t5[6],t6[6],t7[6]])
    F_zz = np.array([t1[7],t2[7],t3[7],t4[7],t5[7],t6[7],t7[7]])
    F_xz = np.array([t1[8],t2[8],t3[8],t4[8],t5[8],t6[8],t7[8]])
    
    F_x_val = np.sum(F_x)
    F_y_val = np.sum(F_y)
    F_z_val = np.sum(F_z)
    F_xx_val = np.sum(F_xx)
    F_xy_val = np.sum(F_xy)
    F_yy_val = np.sum(F_yy)
    F_yz_val = np.sum(F_yz)
    F_zz_val = np.sum(F_zz)
    F_xz_val = np.sum(F_xz)
    
    return np.array([F_x_val, F_y_val, F_z_val, F_xx_val, F_xy_val, F_yy_val, F_yz_val, F_zz_val,F_xz_val])
    