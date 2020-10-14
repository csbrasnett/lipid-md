# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""
import numpy as np
import math
pi=math.pi

def term7(x,y,z,m):
    #t7_x
    cg = (-3 * pi * m * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
          5 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(5 * pi * m * x) - 
          math.cos(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * y) - 
          3 * pi * m * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) + 
          5 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(5 * pi * m * x) + 
          math.cos(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * y) - 
          3 * pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
          5 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.cos(5 * pi * m * x) + 
          math.sin(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * y) - 
          3 * pi * m * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) + 
          5 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi * m * math.sin(5 * pi * m * x) - 
          math.sin(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
          3 * pi * m * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) - 
          math.cos(3 * pi * m * y) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * z) - 
          5 * math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi * m * math.sin(5 * pi * m * x) - 
          3 * pi * m * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) + 
          math.cos(3 * pi * m * y) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * z) + 
          5 * math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi * m * math.cos(5 * pi * m * x) - 
          3 * pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) + 
          math.sin(3 * pi * m * y) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * z) - 
          5 * math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi * m * math.cos(5 * pi * m * x) - 
          3 * pi * m * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) - 
          math.sin(3 * pi * m * y) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * z) + 
          5 * math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi * m * math.sin(5 * pi * m * x))

    #t7_y
    cg0 = (-pi * m * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) - 
           3 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.cos(5 * pi * m * x) - 
           5 * math.cos(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * y) + 
           pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) -
           3 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi * m * math.sin(5 * pi * m * x) + 
           5 * math.cos(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * y) + 
           pi * m * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) - 
           3 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(5 * pi * m * x) - 
           5 * math.sin(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
           pi * m * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           3 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(5 * pi * m * x) + 
           5 * math.sin(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * y) - 
           5 * pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) - 
           3 * math.sin(3 * pi * m * y) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * z) - 
           math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi * m * math.cos(5 * pi * m * x) + 
           5 * pi * m * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) - 
           3 * math.sin(3 * pi * m * y) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * z) + 
           math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi * m * math.sin(5 * pi * m * x) - 
           5 * pi * m * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) - 
           3 * math.cos(3 * pi * m * y) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * z) + 
           math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi * m * math.sin(5 * pi * m * x) + 
           5 * pi * m * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) - 
           3 * math.cos(3 * pi * m * y) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * z) - 
           math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi * m * math.cos(5 * pi * m * x))

    #t7_z
    cg1 = (-5 * pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
           math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi * m * math.cos(5 * pi * m * x) - 
           3 * math.sin(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * y) + 
           5 * pi * m * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) +
           math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.sin(5 * pi * m * x) - 
           3 * math.sin(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * y) -
           5 * pi * m * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) +
           math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi * m * math.sin(5 * pi * m * x) - 
           3 * math.cos(3 * pi * m * z) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * y) + 
           5 * pi * m * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) - 
           math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi * m * math.cos(5 * pi * m * x) - 
           3 * math.cos(3 * pi * m * z) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * y) -
           pi * m * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) - 
           5 * math.cos(3 * pi * m * y) * pi * m * math.cos(pi * m * x) * math.sin(5 * pi * m * z) - 
           3 * math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi * m * math.cos(5 * pi * m * x) + 
           pi * m * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) + 
           5 * math.cos(3 * pi * m * y) * pi * m * math.sin(pi * m * x) * math.cos(5 * pi * m * z) - 
           3 * math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi * m * math.sin(5 * pi * m * x) + 
           pi * m * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) -
           5 * math.sin(3 * pi * m * y) * pi * m * math.cos(pi * m * x) * math.cos(5 * pi * m * z) -
           3 * math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi * m * math.sin(5 * pi * m * x) -
           pi * m * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) + 
           5 * math.sin(3 * pi * m * y) * pi * m * math.sin(pi * m * x) * math.sin(5 * pi * m * z) -
           3 * math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi * m * math.cos(5 * pi * m * x))

    #t7_xx
    cg2 = (-9 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           25 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) -
           math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
           9 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) -
           25 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) -
           math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) + 
           9 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) +
           25 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) + 
           9 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) + 
           25 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) +
           math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) - 
           9 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) -
           math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * z) - 
           25 * math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) -
           9 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) - 
           math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * z) - 
           25 * math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) +
           9 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) + 
           math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * z) + 
           25 * math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) +
           9 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) + 
           math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * z) +
           25 * math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x))

    #t7_xy
    cg3 = (3 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) +
           15 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) +
           5 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) -
           3 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) -
           15 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) +
           5 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) +
           3 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) -
           15 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) +
           5 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) -
           3 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) +
           15 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) +
           5 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) +
           15 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) +
           3 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * z) + 
           5 * math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           15 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) -
           3 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * z) + 
           5 * math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           15 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) +
           3 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * z) + 
           5 * math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           15 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) -
           3 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * z) + 
           5 * math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x))

    #t7_yy
    cg4 = (-pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) -
           9 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) -
           25 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) -
           pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) - 
           9 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) -
           25 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) +
           pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) + 
           9 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) +
           25 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) +
           pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) + 
           9 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) +
           25 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) -
           25 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) -
           9 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * z) - 
           math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           25 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) -
           9 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * z) - 
           math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           25 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) +
           9 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * z) +
           math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           25 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) +
           9 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * z) +
           math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x))

    #t7_yz
    cg5 = (5 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) + 
           3 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           15 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) + 
           5 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           3 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           15 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) +
           5 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) + 
           3 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           15 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) +
           5 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) - 
           3 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           15 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) + 
           5 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) + 
           15 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * z) +
           3 * math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           5 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) - 
           15 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * z) -
           3 * math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           5 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) - 
           15 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * z) +
           3 * math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           5 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) + 
           15 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * z) -
           3 * math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x))

    #t7_zz
    cg6 = (-25 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) - 
           math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           9 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) - 
           25 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) - 
           math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           9 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) + 
           25 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) +
           math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           9 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) +
           25 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) +
           math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           9 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) -
           pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) - 
           25 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * z) - 
           9 * math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) -
           pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) - 
           25 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * z) -
           9 * math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) + 
           25 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * z) + 
           9 * math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) + 
           25 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * z) + 
           9 * math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x))

    #t7_xz
    cg7 = (15 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * y) * math.sin(5 * pi * m * z) + 
           5 * math.cos(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) + 
           3 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * y) - 
           15 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * y) * math.cos(5 * pi * m * z) + 
           5 * math.cos(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           3 * math.sin(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * y) - 
           15 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * y) * math.cos(5 * pi * m * z) + 
           5 * math.sin(3 * pi * m * y) * math.sin(pi * m * z) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           3 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * y) + 
           15 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * y) * math.sin(5 * pi * m * z) + 
           5 * math.sin(3 * pi * m * y) * math.cos(pi * m * z) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           3 * math.cos(3 * pi * m * z) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * y) + 
           3 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.sin(pi * m * z) * math.cos(5 * pi * m * y) + 
           5 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.sin(5 * pi * m * z) + 
           15 * math.sin(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x) - 
           3 * pi ** 2 * m ** 2 * math.sin(3 * pi * m * x) * math.cos(pi * m * z) * math.sin(5 * pi * m * y) + 
           5 * math.cos(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.cos(5 * pi * m * z) - 
           15 * math.sin(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) + 
           3 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.sin(pi * m * z) * math.sin(5 * pi * m * y) + 
           5 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.sin(pi * m * x) * math.cos(5 * pi * m * z) - 
           15 * math.cos(3 * pi * m * z) * math.cos(pi * m * y) * pi ** 2 * m ** 2 * math.cos(5 * pi * m * x) - 
           3 * pi ** 2 * m ** 2 * math.cos(3 * pi * m * x) * math.cos(pi * m * z) * math.cos(5 * pi * m * y) + 
           5 * math.sin(3 * pi * m * y) * pi ** 2 * m ** 2 * math.cos(pi * m * x) * math.sin(5 * pi * m * z) + 
           15 * math.cos(3 * pi * m * z) * math.sin(pi * m * y) * pi ** 2 * m ** 2 * math.sin(5 * pi * m * x))

    return np.array([cg,cg0,cg1,cg2,cg3,cg4,cg5,cg6,cg7])