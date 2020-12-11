#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 11:30:15 2020

@author: Chris
"""

import pickle
import numpy as np


from scipy.spatial import KDTree
from scipy.optimize import least_squares, differential_evolution, shgo, dual_annealing, minimize
from scipy import linalg

from mdma import atom
from skimage import measure
from itertools import product
from multiprocessing import get_context, current_process
import time
import pandas
import argparse
import glob
import ovito as ov
import os

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

def get_surrounding_coords(tree, C5Acoords, index, cut_off_radius):

    surrounding_indicies = np.array(tree.query_ball_point(C5Acoords[index], cut_off_radius))
    
    return C5Acoords[surrounding_indicies] 


def fit_writer(initial_points, fit_result, index, q,co,file):
    print(q, co, file, fit_result)
    
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
        
        fname = os.path.abspath(file).split('.pdb')[0] +'_cutoff_'+str(co)+'_points_out_'+str(index)+'_' + q +'.atom'
        
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

    
def fitting(a,index, cut_off,file):
        
    x = a[0:,0]
    y = a[0:,1]
    z = a[0:,2]
    
    '''
    x0 is the initial guess array for the parameters of the quadric function
    
    x0[0] = P
    x0[1] = Q
    x0[2] = R
    x0[3] = S
    x0[4] = T
    x0[5] = C
    x0[6] = U
    
    b are the bounds on the parameters of the fitting functions
    
    '''   
    # b = [(-1000, 1000),(-1000, 1000),(-1000, 1000),(-1000, 1000),(-1000, 1000),(-1000, 1000)]
    b = [(-1, 1),(-1, 1),(-1, 1),(-1, 1),(-1, 1),(-1, 1)]
    x0 = np.array([1,1,1,1,1,1,1])

    # print('starting fits')

    current = current_process()

    #perform a least squares fit of the quadric form to the point cloud
    start = time.perf_counter()
    
    res_div = differential_evolution(fun, b, args = (x,y,z), maxiter = 500)
    t1 = time.perf_counter()
    # print('t1', current.name, current._identity, t1-start, res_div.nit)
    
    res_shg = shgo(fun, b, args = (x,y,z))#, options = {'maxiter': 10})
    t2 = time.perf_counter()
    # print('t2', current.name, current._identity, t2-t1, res_shg.nit)
    
    res_dan = dual_annealing(fun, b, args = (x,y,z), maxiter = 500)
    t3 = time.perf_counter()
    # print('t3', current.name, current._identity, t3-t2, res_dan.nit)
    
    res_lsq = least_squares(fun, x0, args = (x,y,z), max_nfev = 500)
    t4 = time.perf_counter()
    # print('t4', current.name, current._identity, t4-t3, res_lsq.nfev)
    
    res_min = minimize(fun, x0, args = (x,y,z), options = {'maxiter': 500})
    t5 = time.perf_counter()
    # print('t5', current.name, current._identity, t5-t4, res_min.nit)
    
    # print('ending fits')
    
    #calculate the gaussian curvature from the fit of the parameters
    valKdiv = gaussian_curvature(res_div.x)
    valHdiv = mean_curvature(res_div.x)
    
    valKshg= gaussian_curvature(res_shg.x)
    valHshg = mean_curvature(res_shg.x)
    
    valKdan= gaussian_curvature(res_dan.x)
    valHdan = mean_curvature(res_dan.x)

    valKlsq = gaussian_curvature(res_lsq.x)
    valHlsq = mean_curvature(res_lsq.x)
    
    valKmin = gaussian_curvature(res_min.x)
    valHmin = mean_curvature(res_min.x)
    
    times = np.array([t1-start, t2-t1, t3-t2, t4-t3, t5-t4])
    func_vals = np.array([res_div.fun, res_shg.fun, res_dan.fun, res_lsq.fun[0], res_min.fun])
    n_its = np.array([res_div.nit, res_shg.nit, res_dan.nit, res_lsq.nfev, res_min.nit])
    n_fevs = np.array([res_div.nfev, res_shg.nfev, res_dan.nfev, res_lsq.nfev, res_min.nfev])
    successes = np.array([res_div.success, res_shg.success, res_dan.success, res_lsq.success, res_min.success])
    Ks = np.array([valKdiv, valKshg, valKdan, valKlsq, valKmin])
    Hs = np.array([valHdiv, valHshg, valHdan, valHlsq, valHmin])    
    
    d1 = {'Method': ['Differential Evolution', 'SHGO', 'Dual Annealing', 'Least Squares', 'Minimize'],
         'Time': times, 
         'Function Value': func_vals,
         'No. Iterations': n_its,
         'No. Function Evaluations': n_fevs, 
         'Status': successes,
         'K': Ks,
         'H': Hs}
    
    df = pandas.DataFrame(d1).set_index('Method')
    # print(df)
    # r = np.random.random()
    # print('r', r)
    # if r>0.5:
        # fit_writer(a, res_div.x, index, 'div',cut_off,file)
        # fit_writer(a, res_lsq.x, index, 'lsq',cut_off,file)
        # fit_writer(a, res_min.x, index, 'min',cut_off,file)
        # pickle.dump(df, open(os.path.abspath(file).split('.pdb')[0]+'_cutoff_'+str(cut_off)+'_'+str(index)+'.p', 'wb'))
        
    return df

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

def func(file, cut_off, bead):
        
    coords = file_reader(file, bead)
        
    tree = KDTree(coords)
    
    data_out = {}
    
    for index in range(len(coords)):
        # print(index, cut_off)
        surrounding_coords = get_surrounding_coords(tree, coords, index, cut_off)
        
        try:
            q = np.vstack(surrounding_coords)
            w = len(surrounding_coords[0])
            e = len(surrounding_coords[1])

            '''
            perform PCA on the patch in order to calculate the principal axes 
            of the point cloud. The points will then lie in the frame of the 
            '''
            
            pca_res = PCA(q)
            
            if type(pca_res) == tuple:
                
                a = pca_res[0]

                df = fitting(a, index, cut_off, file)
            
                data_out[index] = df
        
        except IndexError:
            pass
        
    fname = os.path.abspath(file).split('.pdb')[0]+'_comparison_results_cutoff_'+str(cut_off)+'.p'
        
    pickle.dump(data_out, open(fname, 'wb'))
            

def argument_reader():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-b', '--bead', type = str, nargs = '+',  default = ['C5A'], help = 'The beads within the simulation frame to fit parabolas to')
    parser.add_argument('-r', '--radius', type = int, nargs = '+',  default = [30], help = 'Search radius cutoff length for creating point clouds')

    args = parser.parse_args()

    beads = np.array(args.bead)
    radius = np.array(args.radius)
    
    return beads, radius

def make_paramlist(files, ball_point_radii, bead):
    init_list = list(product(f, ball_point_radii))
    
    paramlist = []
    
    for i in init_list:
        paramlist.append(i+(bead,))
        
    return paramlist

if __name__ == '__main__':
    __spec__ = None
    
    f = glob.glob('*.pdb')   

    bead, ball_point_radii = argument_reader()
    
    paramlist = make_paramlist(f, ball_point_radii, bead)
    
    k = len(paramlist)/14
    
    if k < 1:
        csize = 1
    else:
        csize = int(k)
    # print(csize)
    
    with get_context("spawn").Pool(processes = 2) as pool:
        pool.starmap(func, paramlist, chunksize = csize)

        
        
        
