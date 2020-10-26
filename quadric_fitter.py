#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:30:15 2020

@author: Chris
"""

import pickle
import numpy as np

from scipy.spatial import KDTree
from scipy.optimize import shgo
from scipy import linalg

from mdma import atom
from skimage import measure
import glob
from itertools import product
from multiprocessing import get_context
import os
import pandas as pd
import argparse
import ovito as ov

def PCA(data):
    '''
    Perform Principal Component Analysis on a point cloud.
    Subsequently transform the point cloud to the origin and so that it lies 
    in the frame of principal components. 
    '''
    #centering the data
    data -= np.mean(data, axis = 0)  
    
    cov = np.cov(data, rowvar = False)
    try:
        evals , evecs = linalg.eigh(cov)
        
        idx = np.argsort(evals)[::-1]
        evecs = evecs[:,idx]
        evals = evals[idx]
        
        a = np.dot(data, evecs) 
    
        return a, evals, evecs
    
    except ValueError:
        return 0
    

def fun(paras, x, y, z):
    
    result = 0
    for i in range(len(x)):
        result += ((paras[0]*x[i]**2) + (paras[1]*y[i]**2) + (paras[2]*x[i]*y[i])
                   +(paras[3]*x[i]) + (paras[4]*y[i]) + paras[5] - (z[i]))**2

    v = result**0.5

    return v

def quadric(paras, x, y, z):

    t_1 = paras[0]*x**2
    t_2 = paras[1]*y**2 
    t_3 = paras[2]*x*y  
    t_4 = paras[3]*x 
    t_5 = paras[4]*y 
    t_6 = paras[5]
    t_7 = z
    
    t = t_1 + t_2 + t_3 + t_4 + t_5 + t_6 - t_7
    
    return t
    

def gaussian_curvature(paras):
    
    E = 1 + (paras[3]**2)
    F = paras[3]*paras[4]
    G = 1 + (paras[4]**2)
    
    L = 2*paras[0]
    M = paras[2]
    N = 2*paras[1]
    
    nom = (L*N) - (M**2)
    den = (E*G) - (F**2)
    
    K = nom/den
    
    return K

def mean_curvature(paras):
    
    E = 1 + (paras[3]**2)
    F = paras[3]*paras[4]
    G = 1 + (paras[4]**2)
    
    L = 2*paras[0]
    M = paras[2]
    N = 2*paras[1]

    nom = (E*N) - (2*F*M) + (G*L)
    den = (E*G) - (F**2)
    
    H = nom/den
    
    return H



def fit_writer(initial_points, fit_result, index, file, cut_off_radius, bead):
    
    #set up the mgrid in teh right spatial position
    d1 = initial_points[0:,0].mean()
    d2 = initial_points[0:,1].mean()
    d3 = initial_points[0:,2].mean()
    
    e1 = np.abs(initial_points[0:,0].max()-initial_points[0:,0].min()) /2
    e2 = np.abs(initial_points[0:,1].max()-initial_points[0:,1].min()) /2
    e3 = np.abs(initial_points[0:,2].max()-initial_points[0:,2].min()) /2
    
    
    X,Y,Z = np.mgrid[float(d1-e1):float(d1+e1):(2*float(e1)*1j), 
                     float(d2-e2):float(d2+e2):(2*float(e2)*1j), 
                     float(d3-e3):float(d3+e3):(2*float(e3)*1j)]
    
    #compute the value of the fitted result on the mgrid
    t = quadric(fit_result, X, Y, Z)
    
    try:
        #generate point solutions using marching cubes
        vertices, simplices,normals, values = measure.marching_cubes_lewiner(t)
        
        #sort out the vertices and append them to the list of coordinates to write to file
        Xp,Yp,Zp = zip(*vertices)
        
        surface_points = np.array([Xp,Yp,Zp]).T
        
        surface_points -= np.mean(surface_points, axis = 0)
        
        points_out = np.vstack((initial_points, surface_points))
        
        names_out = ['initial']*len(initial_points) + ['fitted']*len(surface_points)
        
        b = list(bead)
        c = ''.join(b)

        fname = os.path.abspath(file).split('.pdb')[0] + '_'+c+'_COradius_'+ str(cut_off_radius) +'_points_out_'+str(index)+'.atom'
        
        with open(fname, 'w') as f:
            atom.write(points_out, np.array([[points_out[0:,0].min(), points_out[0:,0].max()],
                                             [points_out[0:,1].min(), points_out[0:,1].max()],
                                             [points_out[0:,2].min(), points_out[0:,2].max()]]), f, names_out)

    except ValueError as e1:
        print(e1)
        pass
    except RuntimeError as e2:
        print(e2)
        pass
    
def fitting(a, index, file, cut_off_radius, bead):

        
    '''
    x0 is the initial guess array for the parameters of the quadric function
    
    x0[0] = P
    x0[1] = Q
    x0[2] = R
    x0[3] = S
    x0[4] = T
    x0[5] = C
    x0[6] = U
    
    '''
    b = [(-10, 10),(-10, 10),(-10, 10),(-10, 10),(-10, 10),(-10, 10)]

    x = a[0][0:,0]
    y = a[0][0:,1]
    z = a[0][0:,2]
    
    #perform a least squares fit of the quadric form to the point cloud
    res= shgo(fun, b, args = (x,y,z))
    
    # print(res_lsq)
    
    #calculate the gaussian curvature from the fit of the parameters
    valK = gaussian_curvature(res.x)
    valH = mean_curvature(res.x)
    success = res.success
    eval_val = res.fun
    
    r = np.random.random()
    if r>0.999:
        fit_writer(a[0], res.x, index, file, cut_off_radius, bead)
    
    return valK, valH, success, eval_val


def get_surrounding_coords(tree, coords, index, cut_off_radius):

    surrounding_indicies = tree.query_ball_point(coords[index], cut_off_radius)
    
    surrounding_coords = coords[surrounding_indicies]
    
    return surrounding_coords


def file_reader(file, bead, wrap = False):
    pipeline = ov.io.import_file(file)
    
    if wrap == True:
        pipeline.modifiers.append(ov.modifiers.WrapPeriodicImagesModifier())
    
    pipeline.modifiers.append(ov.modifiers.SelectTypeModifier(property = 'Particle Type', types = set(bead)))    
    
    data = pipeline.compute()
    
    a = np.where(data.particles.selection[:]==1)[0]

    pos = np.array(data.particles.positions[:][a])

    # b = list(bead)
    # c = ''.join(b)

    # fname = file.split('.pdb')[0]+'_'+c+'_coords.p'

    # pickle.dump(pos, open(fname, 'wb'))
    
    return pos

def coord_handling(file, cut_off_radius, bead):

    coords = file_reader(file, bead)
    
    tree = KDTree(coords)
    
    K_vals = []
    H_vals = []
    successes = []
    funs = []
    
    for index in range(coords.shape[0]):
        # print(file, index, coords.shape[0])
        
        #find the coordinates within a cutoff radius to form a point cloud.
        surrounding_coords = get_surrounding_coords(tree, coords, index, cut_off_radius)
        
        
        
        '''
        perform PCA on the patch in order to calculate the principal axes 
        of the point cloud. The points will then be transformed to lie in 
        the frame of the principal axes
        '''
        
        a = PCA(surrounding_coords)
        
        
        if type(a) == tuple:

            K_, H_, S_, F_ = fitting(a, index, file, cut_off_radius, bead)
            
            K_vals.append(K_)
            H_vals.append(H_)
            successes.append(S_)
            funs.append(F_)


    d = {'K': K_vals,
         'H': H_vals,
         'Status': successes,
         'Function Value': funs}
    
    
    '''
    for calculating bending modulus:
        the pivotal plane plays a role in determining the distance between lipid pairs.
        So: when using this data, find coords of terminal carbons, headgroup carbons, and C1* beads
        - lipid vector is then the HG-terminal coords,
        - distance between lipids is the distance at the pivotal plane.
        
        this is following the method of Johner et al.  J. Phys. Chem. Lett, (2014) (see the SI)
    
    '''
    
    df = pd.DataFrame(d)

    b = list(bead)
    c = ''.join(b)

    fname = os.path.abspath(file).split('.pdb')[0]+'_'+c+'_results_Rc_'+str(cut_off_radius)+'.p'
    pickle.dump(df, open(fname, 'wb'))
    



def argument_reader():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-b', '--bead', type = str, nargs = '+',  default = ['C5A'], help = 'The beads within the simulation frame to fit parabolas to')
    parser.add_argument('-r', '--radius', type = int, nargs = '+',  default = [30], help = 'Search radius cutoff length for creating point clouds')

    args = parser.parse_args()

    beads = np.array(args.bead)
    radius = np.array(args.radius)
    
    return beads, radius

def make_paramlist(files, ball_point_radii, bead):
    init_list = list(product(files, ball_point_radii))
    
    
    paramlist = []
    
    for i in init_list:
        paramlist.append(i+(bead,))
        
    return paramlist

if __name__ == '__main__':
    f = glob.glob('*.pdb')   
    
    bead, ball_point_radii = argument_reader()
    paramlist = make_paramlist(f, ball_point_radii, bead)
    
    k = len(paramlist)/14
    
    if k < 1:
        csize = 1
    else:
        csize = int(k)
    
    with get_context("spawn").Pool(processes = 14) as pool:
        pool.starmap(coord_handling, paramlist, chunksize = csize)
