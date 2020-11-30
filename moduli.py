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
import lmfit as lm
# import matplotlib.pyplot as plt

def file_reader(file, atoms, wrap = False):
    pipeline = ov.io.import_file(file)

    if wrap == True:
        pipeline.modifiers.append(ov.modifiers.WrapPeriodicImagesModifier())
    
    pipeline.modifiers.append(ov.modifiers.SelectTypeModifier(property = 'Particle Type', types = set(atoms)))    
    
    data = pipeline.compute()
    
    a = np.where(data.particles.selection[:]==1)[0]

    pos = np.array(data.particles.positions[:][a])
    
    mol_types = np.array(data.particles['Molecule Type'][:][a])
    
    return pos, mol_types



def PCA(data):
    '''
    Perform Principal Component Analysis on a point cloud.
    Subsequently transform the point cloud to the origin and so that it lies 
    in the frame of principal components. 
    '''
    #centering the data
    data -= np.mean(data, axis = 0)  
    
    #calculate the covariance
    cov = np.cov(data, rowvar = False)

    #find eigenvalues and eigenvectors of the covariance for the principal components
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

def coord_handling(file, atoms, cut_off_radius, multi = False):
    
    head_bead = np.array(atoms[0])
    neutral_bead = np.array(atoms[1])
    tail_bead = np.array(atoms[2])    
    
    #headgroup positions
    HG_pos, mol_types = file_reader(file, head_bead)
    #neutral surface positions
    NS_pos_temp, mol_types = file_reader(file, neutral_bead)
    #tail carbon positions
    TC_pos_temp, mol_types = file_reader(file, tail_bead)
    
    
    #this section is to account for molecules with >1 tail so that a lipid vector can be constructed 
    #between the mean tail bead position and the headgroup.
    if tail_bead.shape[0]>1 and multi[0] == False:
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
    mol_type = []
    vectors = []
    distances = []
    terminal_vectors = []
    closest_mol_types = []
    
    #perform PCA on point clouds of terminal carbons and find the normal vector from the underlying surface
    for index in range(len(TC_pos)):
        # print(index)
        #find the coordinates within a cutoff radius to form a point cloud.
        surrounding_coords = get_surrounding_coords(tree, TC_pos, index, cut_off_radius)
        
        #in order to calculate the bending modulus, find the indices of the points within 10A

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
            mts = []
            for j in closest_index:
                c = np.linalg.norm(NS_pos[index]-NS_pos[j])
                cs.append(c)
                
                p = TC_pos[index] - TC_pos[j]
                vs.append(p)            
                
                mt = mol_types[j]
                mts.append(mt)
            
            ind.append(index)
            mol_type.append(mol_types[index])
            vectors.append(lipid_vector)
            closest_indices.append(closest_index)
            normals.append(normal)
            distances.append(cs)
            terminal_vectors.append(vs)
            closest_mol_types.append(mts)
    
    d = {'Index': ind,
         'Molecule Type': mol_type,
         'Lipid Vector': vectors,
         'Closest': closest_indices,
         'Normals': normals,
         'Distances': distances,
         'Terminal Vectors': terminal_vectors,
         'Closest Molecule Types': closest_mol_types}
    
    df = pd.DataFrame(d)

    return df

def moduli_distributions(df):
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
    type_array = np.zeros(0)
    
    #for every point in the data...
    for i in range(df.shape[0]):
        # print(i, df.shape[0])
        #get the vector and normalise it.
        this_lipid_vector = df.iloc[i]['Lipid Vector']/np.linalg.norm(df.iloc[i]['Lipid Vector'])
        
        #get the surface vector (it is already normalised)
        this_surface_vector = df.iloc[i]['Normals']
        
        #calculate the tilt angle between the surface and the lipid for the tilt modulus.
        angles = np.append(angles, np.arccos(np.dot(this_lipid_vector, this_surface_vector)))
        
        #Find the type of this lipid
        this_lipid_type = df.iloc[i]['Molecule Type']
        
        #for each of the closest points found...
        for j in range(len(df.iloc[i]['Closest'])):
        
            jj = df.iloc[i]['Closest'][j]
            
            #get the lipid vector of the lipid
            #has to be done by this callable method because of error handling in the analysis: 
            #the DataFrame index and the Index column counting the numbering for the atoms may not match
            close_lipid_vector = df.loc[df['Index'] == jj]['Lipid Vector'].values[0]/np.linalg.norm(df.loc[df['Index'] == jj]['Lipid Vector'].values[0])
            
            #and get the closest lipid surface vector
            close_lipid_surface_vector = df.loc[df['Index'] == jj]['Normals'].values[0]
    
            #find the molecule type of this lipid closest to the central molecule
            close_lipid_type = df.iloc[i]['Closest Molecule Types'][j]
    
            #find the distance between them
            inter_lipid_distance = df.iloc[i]['Distances'][j]

            #and now calculate the S_i value for this lipid
            S_i = (close_lipid_vector + 
                    close_lipid_surface_vector - 
                    this_lipid_vector - 
                    this_surface_vector)/inter_lipid_distance
            
            #use the vector between them...
            inter_lipid_vector = df.iloc[i]['Terminal Vectors'][j]
            
            #a number to indicate the type of the interaction
            type_arr = np.sort(np.array([this_lipid_type, close_lipid_type]))
            type_val = int(str(type_arr[0])+str(type_arr[1]))
            
            #to calculate the splay
            S = np.dot(S_i, inter_lipid_vector)
            S_values = np.append(S_values, S)
            type_array = np.append(type_array, type_val)


        
    d = {'S_values': S_values,
         'angles': angles,
         'types' : type_array}
    
    return d


def func(x,a,b,c):
    return a*((x-b)**2) +c


def modulus_evaluation(data):
    
    try:
        S_values = data['S_values']
    except TypeError:
        return TypeError
    
    #find the number of types of pairs that need to be looped over
    types = data['types']
    unique_types = np.unique(types)
    
    
    angles = data['angles']


    chi = np.zeros(0)
    chi_error = np.zeros(0)
    
    n_pairs_total = np.zeros(0)
    
    '''
    the bending modulus
    '''
    try:
        for i in range(len(unique_types)):
            #find the indicies for this type of pair
            this_type = np.where(types==unique_types[i])[0]
            
            #find the number of pairs involved in this 
            n_pairs = len(this_type)
            n_pairs_total = np.append(n_pairs_total, n_pairs)
            
            #use the distribution of S to calculate chi for this pair type
            S_vals_init = S_values[this_type]
                    
            S_values_final = S_vals_init[np.where((S_vals_init>-1)&(S_vals_init<1))[0]]
            
            P = np.histogram(S_values_final, bins = 100)
            
            bin_mid_points = (P[1][:-1] + P[1][1:])/2
            bin_probs = P[0]/P[0].sum()
            
            gmodel = lm.models.VoigtModel()
            pars = gmodel.guess(bin_probs, x=bin_mid_points)
            out = gmodel.fit(bin_probs, pars, x=bin_mid_points)
            
    
            op_mu = out.best_values['center']
            op_sig = out.best_values['sigma']
            
            LHS = -2 * np.log(bin_probs)        
    
            
            n_sigma = 2
            a_ = bin_mid_points[np.where((bin_mid_points>(op_mu-(n_sigma*op_sig)))&(bin_mid_points<(op_mu+(n_sigma*op_sig))))[0]]
            b_ = LHS[np.where((bin_mid_points>(op_mu-(n_sigma*op_sig)))&(bin_mid_points<(op_mu+(n_sigma*op_sig))))[0]]
            
            a_ = a_[~np.isnan(a_)]
            b_ = b_[~np.isnan(b_)]
            
            hmodel = lm.Model(func, nan_policy='omit')
            
            
            result = hmodel.fit(b_,x=a_, a=10, b= -0.2, c=0)
            
            res_a = result.params['a'].value
            res_b = result.params['b'].value
            res_c = result.params['c'].value
            err_a = result.params['a'].stderr
            
            chi = np.append(chi, res_a)
            chi_error = np.append(chi_error, err_a)

        #calculate the total bending modulus for the system
        total_pairs = np.sum(n_pairs_total)
        X = (1/total_pairs)*np.sum(n_pairs_total/chi)
        Kc = 1/X
        
        #work out the error in X from chi and then finally error in Kc 
        X_err = np.sqrt(np.sum((np.abs((n_pairs_total/total_pairs)*(-1/(chi**2)))**2)*(chi_error**2)))
        Kc_err = X_err/(X**2)
    
    except TypeError as e1:
        Kc = 0
        Kc_err = 0
        pass
        
    '''
    the tilt modulus
    '''

    
    angles_conv = -np.abs(angles-(np.pi/2))+np.pi/2
            
    P1 = np.histogram(angles_conv, bins = 100)
    bin_mid_points1 = (P1[1][:-1] + P1[1][1:])/2
    bin_probs1 = P1[0]/P1[0].sum()
    
    bin_mid_points1_degree = bin_mid_points1*(180/np.pi)
    
    
    gmodel1 = lm.models.SkewedGaussianModel()
    pars1 = gmodel1.guess(bin_probs1, x=bin_mid_points1)
    out1 = gmodel1.fit(bin_probs1, pars1, x=bin_mid_points1)
    
    
    op_sig1 = out1.best_values['sigma']
    op_alp1 = out1.best_values['gamma']
    
    delta = op_alp1/np.sqrt(1+(op_alp1**2))
    mu_z = np.sqrt(2/np.pi)*delta
    sigma_z = np.sqrt(1-(mu_z**2))
    
    gamma1_factor = ((4-np.pi)/2)
    gamma1_nom = (delta*np.sqrt(2/np.pi))**3
    gamma1_den = (1- 2*delta**2 /np.pi)**(3/2)
    gamma1 = gamma1_factor*(gamma1_nom/gamma1_den)
    
    Mo = mu_z - (gamma1*sigma_z)/2 - ((np.sign(op_alp1)/2)*np.exp((-2*np.pi)/np.abs(op_alp1)))
        
    sig = np.sqrt((op_sig1**2)*(1-(2*(delta**2))/np.pi))
        
    LHS_KbT = np.log(bin_probs1/np.sin(bin_mid_points1)) * -2
    
    Mo_deg = Mo*(180/np.pi)
    sig_deg = sig*(180/np.pi)
    
    a1_ = bin_mid_points1[np.where((bin_mid_points1>(Mo-(sig)))&(bin_mid_points1<(Mo+(sig))))[0]]
    b1_ = LHS_KbT[np.where((bin_mid_points1>(Mo-(sig)))&(bin_mid_points1<(Mo+(sig))))[0]]
    
    
    a1_ = a1_[~np.isnan(a1_)]
    b1_ = b1_[~np.isnan(b1_)]
    
    hmodel1 = lm.Model(func)
    result1 = hmodel1.fit(b1_,x=a1_, a=10, b= -0.2, c=0)
    
    res_a1 = result1.params['a'].value
    res_b1 = result1.params['b'].value
    res_c1 = result1.params['c'].value
    err_a1 = result1.params['a'].stderr
    
    Kt = res_a1
    Kt_err = err_a1
        
        
        
        # '''
        # make the plot
        # '''
        # fig, axarr = plt.subplots(2,2, figsize = (10,10),sharex = 'col')
        
        # axarr[0,0].scatter(bin_mid_points, bin_probs)
        # axarr[0,0].set_ylabel('P(S)')
        # axarr[0,0].set_ylim(0,bin_probs.max()+0.005)
        # axarr[0,0].axvline(op_mu-(n_sigma*op_sig), ls = '--', c = '#e79d24')
        # axarr[0,0].axvline(op_mu+(n_sigma*op_sig), ls = '--', c = '#e79d24')
        
        # axarr[1,0].scatter(bin_mid_points, LHS)
        # xp = np.linspace(op_mu-(4*op_sig), op_mu+(4*op_sig), 100)
        
        # axarr[1,0].plot(xp, func(xp, res_a, res_b, res_c), c= '#b6253a')
        # axarr[1,0].axvline(op_mu-(n_sigma*op_sig), ls = '--', c = '#e79d24')
        # axarr[1,0].axvline(op_mu+(n_sigma*op_sig), ls = '--', c = '#e79d24')
        # axarr[1,0].text(0.05,0.05, 'K$_{c}$ = %.3f [k$_{B}$T]' %res_a, transform = axarr[1,0].transAxes)
        # axarr[1,0].set_ylabel('PMF(S) [k$_{B}$T]')
        # axarr[1,0].set_xlabel('S')
        
        
        # axarr[0,1].scatter(bin_mid_points1_degree, bin_probs1)
        # axarr[0,1].set_ylabel(r'P($\theta$)')
        # axarr[0,1].set_ylim(0,bin_probs1.max()+0.005)
        # axarr[0,1].axvline(Mo_deg-(sig_deg), ls = '--', c = '#e79d24')
        # axarr[0,1].axvline(Mo_deg+(sig_deg), ls = '--', c = '#e79d24')
        # axarr[0,1].yaxis.set_label_position("right")
        # axarr[0,1].yaxis.tick_right()
        
        # axarr[1,1].scatter(bin_mid_points1_degree, LHS_KbT)
        # xp1 = np.linspace(0,np.pi/2, 100)
        # axarr[1,1].plot(bin_mid_points1_degree, func(xp1, res_a1, res_b1, res_c1), c= '#b6253a')
        # axarr[1,1].axvline(Mo_deg-(sig_deg), ls = '--', c = '#e79d24')
        # axarr[1,1].axvline(Mo_deg+(sig_deg), ls = '--', c = '#e79d24')
        # axarr[1,1].text(0.7,0.05, r'$\chi_{c}$ = %.3f [k$_{B}$T]' %res_a1, transform = axarr[1,1].transAxes)
        # axarr[1,1].set_ylabel(r'PMF($\theta$) [k$_{B}$T]')
        # axarr[1,1].set_xlabel(r'$\theta$')
        # axarr[1,1].yaxis.set_label_position("right")
        # axarr[1,1].yaxis.tick_right()
        
        
        # fig.savefig(file.split('.p')[0]+'.png', dpi =200)
        
        # plt.close(fig)
        
    return np.array([Kc, Kc_err, Kt, Kt_err])







def calculation(file, radius, atoms, multi):
    
    df = coord_handling(file, atoms, radius, multi)    

    distribution_results = moduli_distributions(df)

    results = modulus_evaluation(distribution_results)
    
    fname = os.path.abspath(file).split('.pdb')[0]+'_cutoff_'+str(radius)+'_moduli.p'
 
    pickle.dump(results, open(fname, 'wb'))
   


def argument_reader():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-e', '--head', type = str, nargs = '+',  default = ['ETH'], help = 'The headgroup atom of the lipid')
    parser.add_argument('-n', '--neutral', type = str, nargs = '+',  default = ['C1A'], help = 'The atom of the lipid at the neutral surface of the bilayer')
    parser.add_argument('-t', '--tail', type = str, nargs = '+',  default = ['C5A'], help = 'The terminal carbon atom of the lipid')
    parser.add_argument('-r', '--radius', type = int, nargs = '+',  default = [30], help = 'Search radius cutoff length for creating point clouds')
    parser.add_argument('-m', '--multi', type = bool, nargs = '+',  default = False, help = 'Boolean flag for indicating whether multiple lipid species are being investigated')

    args = parser.parse_args()
    
    
    atoms = [list(args.head), list(args.neutral), list(args.tail)]
    
    radius = np.array(args.radius)
    
    multi = args.multi
    
    return atoms, radius, multi

if __name__=='__main__':

    files = glob.glob('*.pdb')
    
    atoms, radius, multi = argument_reader()
        
    prod = list(product(files, radius))
    
    paramlist = []
    for i in prod:
        paramlist.append(i+(atoms,)+ (multi,))

    k = len(paramlist)/14
        
    if k < 1:
        csize = 1
    else:
        csize = int(k)
        
    with get_context("spawn").Pool(processes = 14) as pool:
        pool.starmap(calculation, paramlist, chunksize = csize)


    
    
    
    
    
    

