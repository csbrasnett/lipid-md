#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 15:46:55 2020

@author: Chris
"""

import ovito as ov
import numpy as np
import glob
from scipy.optimize import leastsq
import pickle

def get_coords(f):
    
    pipeline = ov.io.import_file(f)
    d = pipeline.compute()
    d.cell_.pbc = (False, False, False)
    
    p1 = ov.pipeline.Pipeline(source = ov.pipeline.StaticSource(data = d))
    
    p1.modifiers.append(ov.modifiers.SelectTypeModifier(types = {'W'}))
    p1.modifiers.append(ov.modifiers.ClusterAnalysisModifier(cutoff = 8.6, only_selected = True, sort_by_size = True))
    
    data = p1.compute()
    clusters = data.particle_properties.cluster[...]
    positions = data.particles.positions[...]

    cylinders = {}
    
    for i in range(1,8):
        cluster_positions = np.where(clusters == i)[0]
        d = np.array(positions[cluster_positions])

        # print(f[-10:], d,d.shape, type(d))

        cylinders[i] = d

    return cylinders

def get_solid_vol(file):
    pipeline = ov.io.import_file(file)

    #select W particles in the system.            
    pipeline.modifiers.append(ov.modifiers.SelectTypeModifier(property = 'Particle Type',
                                                              types = {'W'}))
    
    vol_frac = np.zeros(0)    

    mod1 = ov.modifiers.ConstructSurfaceModifier(only_selected = True, 
                                                  radius = 1, 
                                                  smoothing_level = 8, 
                                                  identify_regions = True)


    pipeline.modifiers.append(mod1)

    loops = np.linspace(9,12,5)

    for i in loops:
        mod1.radius = np.round(i, decimals = 1)
        
        data = pipeline.compute()
                
        solid_volume = data.attributes['ConstructSurfaceMesh.filled_volume']
        cell_volume = data.attributes['ConstructSurfaceMesh.cell_volume']
        
        fraction = solid_volume/cell_volume
        vol_frac = np.append(vol_frac, fraction)
    
    return vol_frac

    

def cylinderFitting(xyz):

    """
    This is a fitting for a vertical cylinder fitting
    Reference:
    http://www.int-arch-photogramm-remote-sens-spatial-inf-sci.net/XXXIX-B5/169/2012/isprsarchives-XXXIX-B5-169-2012.pdf

    xyz is a matrix contain at least 5 rows, and each row stores x y z of a cylindrical surface
    p is initial values of the parameter;
    p[0] = Xc, x coordinate of the cylinder centre
    P[1] = Yc, y coordinate of the cylinder centre
    P[2] = alpha, rotation angle (radian) about the x-axis
    P[3] = beta, rotation angle (radian) about the y-axis
    P[4] = r, radius of the cylinder

    th, threshold for the convergence of the least squares

    """   
    x = xyz[:,0]
    y = xyz[:,1]
    z = xyz[:,2]


    p = np.array([x.mean(), y.mean(), 0,0, 10])

    fitfunc = lambda p, x, y, z: (- np.cos(p[3])*(p[0] - x) - z*np.cos(p[2])*np.sin(p[3]) - np.sin(p[2])*np.sin(p[3])*(p[1] - y))**2 + (z*np.sin(p[2]) - np.cos(p[2])*(p[1] - y))**2 #fit function
    errfunc = lambda p, x, y, z: fitfunc(p, x, y, z) - p[4]**2 #error function 
    try:
        est_p , success = leastsq(errfunc, p, args=(x, y, z), maxfev=1000)
        
        return est_p, success
        
    except TypeError:
        return 0
        pass
        

if __name__=="__main__":

    
    
    '''
    NB: must be using water for both the cylinder fitting and the surface mesh
    calculation in order to calculate the lattice parameter of the hexagonal phase
    
    this script can be used to calculate the radius of the pivotal surface of the
    HII mesophase if only looking at the radius, and selecting the C1A/C1B (Martini CG) atoms
    '''
    
    
    files = glob.glob('*.pdb')
    
    all_radii = np.zeros(0)
    all_radii_std = np.zeros(0)
    vol_fracs = np.zeros(0)
    
    for i in files:
        print(i)
        xyzs = get_coords(i)
        
        vol_fracs = np.append(vol_fracs, get_solid_vol(i))
        
        radii = np.zeros(0)
        print(vol_fracs)
        for j in xyzs.keys():
            xyz = xyzs[j]
            
            est_p, success = cylinderFitting(xyz)
            
            if len(est_p) == 5:
                # print(est_p[4])
                radii = np.append(radii, est_p[4])
                
        print(radii)
        # print(radii.mean(), radii.std())
        all_radii = np.append(all_radii, radii.mean())
        all_radii_std = np.append(all_radii_std, radii.std())
    
    all_radii = all_radii[~np.isnan(all_radii)]
    all_radii_std = all_radii_std[~np.isnan(all_radii_std)]

    
    d_w = 2*all_radii.mean()
    d_w_err = all_radii_std.mean()/np.sqrt(len(all_radii_std))
        
    a = np.sqrt((1/vol_fracs.mean())*(np.pi/(2*np.sqrt(3)))*(d_w**2))

    a_err = np.sqrt((((a/(2*vol_fracs.mean()))**2)*((vol_fracs.std()/np.sqrt(len(vol_fracs)))**2) +
             (((a/d_w)**2)*(d_w_err**2))))
    
    d = {'radius': all_radii.mean(),
         'radius error':  all_radii_std.mean()/np.sqrt(len(all_radii_std)),
         'a': a,
         'a error': a_err             
         }
    
    pickle.dump(d, open('cylinder_results.p', 'wb'))
    
    
    
    
    
    