#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 15:53:33 2020

@author: Chris
"""

import pickle
import pandas as pd
import numpy as np
from scipy import linalg
import ovito as ov
from scipy.spatial import KDTree
import glob
import os
import argparse
from itertools import product
from multiprocessing import get_context

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
    

def get_surrounding_coords(tree, coords, index, cut_off_radius):

    surrounding_indicies = tree.query_ball_point(coords[index], cut_off_radius)
    
    surrounding_coords = coords[surrounding_indicies]
    
    return surrounding_coords

def coord_handling(file, beads, cut_off_radius):
        
    
    head_bead = np.array(beads[0])
    neutral_bead = np.array(beads[1])
    tail_bead = np.array(beads[2])    
    
    # print('here', file, beads[0], type(beads[0]),len(beads[0]), head_bead, type(head_bead), head_bead.shape)
    
    #headgroup positions
    HG_pos = file_reader(file, head_bead)
    #neutral surface positions
    NS_pos_temp = file_reader(file, neutral_bead)
    #tail carbon positions
    TC_pos_temp = file_reader(file, tail_bead)
    
    
    #this section is to account for molecules with >1 tail so that a lipid vector can be constructed 
    #between the mean tail bead position and the headgroup.
    if tail_bead.shape[0]>1:
        TC_pos = np.zeros((0,3))
        NS_pos = np.zeros((0,3))
        for i in range(int(TC_pos_temp.shape[0]/2)):
            effective_tail_pos = np.mean(np.vstack((TC_pos_temp[(2*i)], TC_pos_temp[(2*i)+1])), axis = 0)
            effective_neutral_pos = np.mean(np.vstack((NS_pos_temp[(2*i)], NS_pos_temp[(2*i)+1])), axis = 0)
        
            TC_pos = np.vstack((TC_pos,effective_tail_pos))
            NS_pos = np.vstack((NS_pos,effective_neutral_pos))
    else:
        TC_pos = TC_pos_temp
        NS_pos = NS_pos_temp
    
    
    #create a KDTree from the terminal carbon pos 
    tree = KDTree(TC_pos)
    
    normals = []
    closest_indices = []
    ind = []
    vectors = []
    distances = []
    terminal_vectors = []
    
    #perform PCA on point clouds of terminal carbons and find the normal vector from the underlying surface
    for index in range(len(TC_pos)):
        # print(index)
        #find the coordinates within a cutoff radius to form a point cloud.
        surrounding_coords = get_surrounding_coords(tree, TC_pos, index, cut_off_radius)
        
        #in order to calculate the bending modulus, find the indices of the points within 10A
        #the query method will always return the point of interest at index 0
        # closest_index = tree.query(TC_pos[index],11)[1][1:]
        
        closest_index = np.array(tree.query_ball_point(TC_pos[index], 10))
        mask = np.ones(len(closest_index), dtype = bool)
        mask[np.where(closest_index==index)]=False
        closest_index = closest_index[mask]
        
        
        '''
        perform PCA on the patch in order to calculate the principal axes 
        of the point cloud. The points will then be transformed to lie in 
        the frame of the principal axes
        '''
        
        a = PCA(surrounding_coords)
        
        
        if type(a) == tuple:
            #get the normal vector of the point cloud as determined by the PCA
            normal = a[2][0:,2]
            
            #calculate the lipid vector
            lipid_vector = HG_pos[index] - TC_pos[index]
            
            #for every point within 10A, calculate 
            #the distance between the lipids at the neutral plane
            #and the vector between the terminal carbon beads
            cs = []
            vs = []
            for j in closest_index:
                c = np.linalg.norm(NS_pos[index]-NS_pos[j])
                cs.append(c)
                
                p = TC_pos[index] - TC_pos[j]
                vs.append(p)            
            
            ind.append(index)
            vectors.append(lipid_vector)
            closest_indices.append(closest_index)
            normals.append(normal)
            distances.append(cs)
            terminal_vectors.append(vs)
    
    d = {'Index': ind,
          'Lipid Vector': vectors,
          'Closest': closest_indices,
          'Normals': normals,
          'Distances': distances,
          'Terminal Vectors': terminal_vectors}
    
    df = pd.DataFrame(d)

    return df

def moduli_calculations(df):
    '''

    Parameters
    ----------
    df : Pandas DataFrame
        DESCRIPTION.

    Returns
    -------
    S_values : Numpy array
        An array of lipid splays

    '''
    S_values = np.zeros(0)
    angles = np.zeros(0)
    
    #for every point in the data...
    for i in range(df.shape[0]):
        # print(i, df.shape[0])
        #get the vector and normalise it.
        this_lipid_vector = df.iloc[i]['Lipid Vector']/np.linalg.norm(df.iloc[i]['Lipid Vector'])
        
        #get the surface vector (it is already normalised)
        this_surface_vector = df.iloc[i]['Normals']
        
        #calculate the tilt angle between the surface and the lipid for the tilt modulus.
        angles = np.append(angles, np.arccos(np.dot(this_lipid_vector, this_surface_vector)))
        
        #for each of the closest points found...
        for j in range(len(df.iloc[i]['Closest'])):
            
            jj = df.iloc[i]['Closest'][j]
            
            #get the lipid vector of the lipid
            #has to be done by this callable method because of error handling in the analysis: 
            #the DataFrame index and the Index column counting the numbering for the atoms may not match
            close_lipid_vector = df.loc[df['Index'] == jj]['Lipid Vector'].values[0]/np.linalg.norm(df.loc[df['Index'] == jj]['Lipid Vector'].values[0])
            
            #and get the closest lipid surface vector
            close_lipid_surface_vector = df.loc[df['Index'] == jj]['Normals'].values[0]
    
            #find the distance between them
            inter_lipid_distance = df.iloc[i]['Distances'][j]

            #and now calculate the S_i value for this lipid
            S_i = (close_lipid_vector + 
                    close_lipid_surface_vector - 
                    this_lipid_vector - 
                    this_surface_vector)/inter_lipid_distance
            
            #use the vector between them...
            inter_lipid_vector = df.iloc[i]['Terminal Vectors'][j]
            
            #to calculate the splay
            S = np.dot(S_i, inter_lipid_vector)
            S_values = np.append(S_values, S)
            
            
            
    return S_values, angles     

def calculation(file, radius, beads):
    # print(file, beads, radius)
    
    df = coord_handling(file, beads, radius)    
    S_values, angles = moduli_calculations(df)
    
    d = {'S_values': S_values,
          'angles': angles}
    
    # print(file, radius, type(file), type(radius))
    
    fname = os.path.abspath(file).split('.pdb')[0]+'_cutoff_'+str(radius)+'_moduli.p'
    # print(fname)
    pickle.dump(d, open(fname, 'wb'))
   


def argument_reader():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-e', '--head', type = str, nargs = '+',  default = ['ETH'], help = 'The headgroup atom of the lipid')
    parser.add_argument('-n', '--neutral', type = str, nargs = '+',  default = ['C1A'], help = 'The atom of the lipid at the neutral surface of the bilayer')
    parser.add_argument('-t', '--tail', type = str, nargs = '+',  default = ['C5A'], help = 'The terminal carbon atom of the lipid')
    parser.add_argument('-r', '--radius', type = int, nargs = '+',  default = [30], help = 'Search radius cutoff length for creating point clouds')

    args = parser.parse_args()
    
    
    beads = [list(args.head), list(args.neutral), list(args.tail)]
    
    radius = np.array(args.radius)
    
    return beads, radius 

if __name__=='__main__':

    files = glob.glob('*.pdb')
    
    beads, radius = argument_reader()
    # print(beads, type(beads))
    # for i in beads:
    #     print(i, type(i))
        
    prod = list(product(files, radius))
    
    paramlist = []
    for i in prod:
        paramlist.append(i+(beads,))

    k = len(paramlist)/14
        
    if k < 1:
        csize = 1
    else:
        csize = int(k)
    
    with get_context("spawn").Pool(processes = 14) as pool:
        pool.starmap(calculation, paramlist, chunksize = csize)


    
    
    
    
    
    

