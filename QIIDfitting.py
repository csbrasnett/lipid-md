#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from mdma import atom
from skimage import measure
from QIIDcurvature import main as curvature
import ovito as ov
import pickle
import glob
from itertools import product
from multiprocessing import get_context

'''these initial functions are all associated with fitting minimal surface nodal approximations'''
def trig_choose(a,var):
    if a == 's':
        return np.sin(var)
    if a == 'c':
        return np.cos(var)
        
def O_pqr(p,q,r,h,k,l,x,y,z,lamb):
    X=lamb*np.pi*x
    Y=lamb*np.pi*y
    Z=lamb*np.pi*z
    
    term_1=trig_choose(p,(h*X))*trig_choose(q,(k*Z))*trig_choose(r,(l*Y))
    term_2=trig_choose(p,(h*Y))*trig_choose(q,(k*X))*trig_choose(r,(l*Z))
    term_3=trig_choose(p,(h*Z))*trig_choose(q,(k*Y))*trig_choose(r,(l*X))    
    
    return term_1+term_2+term_3
    
def E_pqr(p,q,r,h,k,l,x,y,z,lamb):
    X=lamb*np.pi*x
    Y=lamb*np.pi*y
    Z=lamb*np.pi*z
    
    term_1=trig_choose(p,(h*X))*trig_choose(q,(k*Y))*trig_choose(r,(l*Z))
    term_2=trig_choose(p,(h*Y))*trig_choose(q,(k*Z))*trig_choose(r,(l*X))
    term_3=trig_choose(p,(h*Z))*trig_choose(q,(k*X))*trig_choose(r,(l*Y))    
    
    return term_1+term_2+term_3

def permE1(h,k,l,x,y,z,lamb):
    return E_pqr('c','c','c',h,k,l,x,y,z,lamb)

def permE2(h,k,l,x,y,z,lamb):
    return E_pqr('c','s','s',h,k,l,x,y,z,lamb)

def permE3(h,k,l,x,y,z,lamb):
    return E_pqr('s','c','s',h,k,l,x,y,z,lamb)

def permE4(h,k,l,x,y,z,lamb):
    return E_pqr('s','s','c',h,k,l,x,y,z,lamb)

def permE5(h,k,l,x,y,z,lamb):
    return E_pqr('s','c','c',h,k,l,x,y,z,lamb)

def permE6(h,k,l,x,y,z,lamb):
    return E_pqr('c','s','c',h,k,l,x,y,z,lamb)

def permE7(h,k,l,x,y,z,lamb):
    return E_pqr('s','s','s',h,k,l,x,y,z,lamb)

def permE8(h,k,l,x,y,z,lamb):
    return E_pqr('c','c','s',h,k,l,x,y,z,lamb)

def permO1(h,k,l,x,y,z,lamb):
    return O_pqr('c','c','c',h,k,l,x,y,z,lamb)
    
def permO2(h,k,l,x,y,z,lamb):
    return O_pqr('c','s','s',h,k,l,x,y,z,lamb)
    
def permO3(h,k,l,x,y,z,lamb):
    return O_pqr('s','c','s',h,k,l,x,y,z,lamb)
    
def permO4(h,k,l,x,y,z,lamb):
    return O_pqr('s','s','c',h,k,l,x,y,z,lamb)

def permO5(h,k,l,x,y,z,lamb):
    return O_pqr('c','s','c',h,k,l,x,y,z,lamb)

def permO6(h,k,l,x,y,z,lamb):
    return O_pqr('c','c','s',h,k,l,x,y,z,lamb)

def permO7(h,k,l,x,y,z,lamb):
    return O_pqr('s','s','s',h,k,l,x,y,z,lamb)

def permO8(h,k,l,x,y,z,lamb):
    return O_pqr('s','c','c',h,k,l,x,y,z,lamb)

def D_nodal_approx(x,y,z,lamb):
    
    term1=permE1(1,1,1,x,y,z,lamb)+permE2(1,1,1,x,y,z,lamb)+permE3(1,1,1,x,y,z,lamb)+permE4(1,1,1,x,y,z,lamb)
   
    term2=0.1407*(permE1(1,1,1,x,y,z,lamb)+permE2(1,1,1,x,y,z,lamb)+
                  permE3(1,1,1,x,y,z,lamb)+permE4(1,1,1,x,y,z,lamb)+
                  permO1(1,1,1,x,y,z,lamb)+permO2(1,1,1,x,y,z,lamb)+
                  permO3(1,1,1,x,y,z,lamb)+permO4(1,1,1,x,y,z,lamb))
    
    term3=0.0200*(permE1(3,1,1,x,y,z,lamb)+permE2(3,1,1,x,y,z,lamb)-
                  permE3(3,1,1,x,y,z,lamb)+permE4(3,1,1,x,y,z,lamb)+
                  permO1(3,1,1,x,y,z,lamb)+permO2(3,1,1,x,y,z,lamb)-
                  permO3(3,1,1,x,y,z,lamb)-permO4(3,1,1,x,y,z,lamb))
    
    term4=0.0138*(permE1(3,1,3,x,y,z,lamb)-permE2(3,1,3,x,y,z,lamb)+
                  permE3(3,1,3,x,y,z,lamb)-permE4(3,1,3,x,y,z,lamb)+
                  permO1(3,1,3,x,y,z,lamb)-permO2(3,1,3,x,y,z,lamb)+
                  permO3(3,1,3,x,y,z,lamb)+permO4(3,1,3,x,y,z,lamb))

    term5=0.0028*(permE1(1,1,5,x,y,z,lamb)+permE2(1,1,5,x,y,z,lamb)+
                  permE3(1,1,5,x,y,z,lamb)+permE4(1,1,5,x,y,z,lamb)+
                  permO1(1,1,5,x,y,z,lamb)+permO2(1,1,5,x,y,z,lamb)+
                  permO3(1,1,5,x,y,z,lamb)+permO4(1,1,5,x,y,z,lamb))
    
    term6=0.0021*(permE1(3,3,3,x,y,z,lamb)+permE2(3,3,3,x,y,z,lamb)+
                  permE3(3,3,3,x,y,z,lamb)+permE4(3,3,3,x,y,z,lamb)+
                  permO1(3,3,3,x,y,z,lamb)+permO2(3,3,3,x,y,z,lamb)+
                  permO3(3,3,3,x,y,z,lamb)+permO4(3,3,3,x,y,z,lamb))
    
    term7=0.0001*(permE1(3,1,5,x,y,z,lamb)+permE2(3,1,5,x,y,z,lamb)-
                  permE3(3,1,5,x,y,z,lamb)-permE4(3,1,5,x,y,z,lamb)+
                  permO1(3,1,5,x,y,z,lamb)+permO2(3,1,5,x,y,z,lamb)-
                  permO3(3,1,5,x,y,z,lamb)-permO4(3,1,5,x,y,z,lamb))
    
    return term1+term2+term3-term4-term5-term6+term7
    

def P_nodal_approx(x,y,z,lamb):
    
    term1=permE1(1,0,0,x,y,z,lamb)
    
    term2=0.2260*(permE1(1,0,0,x,y,z,lamb)+permO1(1,0,0,x,y,z,lamb))
    
    term3=0.0516*(permE1(1,1,1,x,y,z,lamb)+permO1(1,1,1,x,y,z,lamb))
    
    term4=0.0196*(permE1(2,1,0,x,y,z,lamb)+permO1(2,1,0,x,y,z,lamb))
    
    term5=0.0027*(permE1(3,0,0,x,y,z,lamb)+permO1(3,0,0,x,y,z,lamb))
    
    return term1+term2-term3-term4-term5


def G_nodal_approx(x,y,z,lamb):
    
    term1=permE5(1,1,0,x,y,z,lamb)+permO5(1,1,0,x,y,z,lamb)
    
    term2=0.3435*(permE4(1,1,0,x,y,z,lamb)+permO5(1,1,0,x,y,z,lamb))
    
    term3=0.3435*(permE6(0,3,1,x,y,z,lamb)+permO6(0,3,1,x,y,z,lamb))
    
    term4=0.3435*(permE7(2,2,2,x,y,z,lamb)+permO7(2,2,2,x,y,z,lamb))
    
    term5=0.3435*(permE6(2,1,3,x,y,z,lamb)+permO6(2,1,3,x,y,z,lamb))
    
    term6=0.3435*(permE8(3,0,3,x,y,z,lamb)+permO8(3,0,3,x,y,z,lamb))
    
    term7=0.3435*(permE4(1,1,4,x,y,z,lamb)+permO5(1,1,4,x,y,z,lamb))
    
    return term1+term2-term3-term4+term5+term6-term7
    

'''performing a rotational and transformational coordinate transform on points x,y,z by 6 values om the vals array'''
def translations(x,y,z,vals):
    #vals[0:3] are x,y,z translational translations
    translational_matrix=np.matrix(np.vstack((vals[0]*np.ones_like(x),vals[1]*np.ones_like(x),vals[2]*np.ones_like(x))))

    #vals[3:] are angles to rotate about x y z axes.
    rot_x=np.matrix([[1,0,0],[0,np.cos(vals[3]),-np.sin(vals[3])],[0,np.sin(vals[3]), np.cos(vals[3])]])
    rot_y=np.matrix([[np.cos(vals[4]),0,np.sin(vals[4])],[0,1,0],[-np.sin(vals[4]),0,np.cos(vals[4])]])
    rot_z=np.matrix([[np.cos(vals[5]),-np.sin(vals[5]),0],[np.sin(vals[5]),np.cos(vals[5]),0],[0,0,1]])

    rot=rot_z*rot_y*rot_x
    
    #perform the translation.
    pos_transform=rot*(np.vstack((x,y,z))+translational_matrix)
    
    #nb this is returning a numpy matrix object.
    return pos_transform
    
'''function to minimize in nodal approximation'''
def sigmasq(vals,surf,x,y,z,final=False):

    #do a translation
    transformation=translations(x,y,z,vals)

    #try to fit a surface depending on an initial guess
    if surf == 'D':
        C_i=D_nodal_approx(transformation[0].A[0], transformation[1].A[0],  transformation[2].A[0], 1/vals[6])
    
    elif surf == 'P':
        C_i=P_nodal_approx(transformation[0].A[0], transformation[1].A[0],  transformation[2].A[0], 1/vals[6])
    
    elif surf == 'G':
        C_i=G_nodal_approx(transformation[0].A[0], transformation[1].A[0],  transformation[2].A[0], 1/vals[6])
    
    #evaluate how well the function has been fitted
    C_i_squared=np.square(C_i)
    
    sigma_squared=np.sum(C_i_squared)/len(C_i)
        
    if final==True:
        return C_i, transformation, vals
    
    return sigma_squared

def point_generator(x_mean,y_mean,z_mean,l,lamb):
    
    #define the mgrid from which to generate the points
    X,Y,Z = np.mgrid[x_mean-l:x_mean+l:(2*l * 1j), 
                     y_mean-l:y_mean+l:(2*l * 1j), 
                     z_mean-l:z_mean+l:(2*l * 1j)]

    term1=permE1(1,1,1,X,Y,Z,lamb)+permE2(1,1,1,X,Y,Z,lamb)+permE3(1,1,1,X,Y,Z,lamb)+permE4(1,1,1,X,Y,Z,lamb)
   
    term2=0.1407*(permE1(1,1,1,X,Y,Z,lamb)+permE2(1,1,1,X,Y,Z,lamb)+
                  permE3(1,1,1,X,Y,Z,lamb)+permE4(1,1,1,X,Y,Z,lamb)+
                  permO1(1,1,1,X,Y,Z,lamb)+permO2(1,1,1,X,Y,Z,lamb)+
                  permO3(1,1,1,X,Y,Z,lamb)+permO4(1,1,1,X,Y,Z,lamb))
    
    term3=0.0200*(permE1(3,1,1,X,Y,Z,lamb)+permE2(3,1,1,X,Y,Z,lamb)-
                  permE3(3,1,1,X,Y,Z,lamb)+permE4(3,1,1,X,Y,Z,lamb)+
                  permO1(3,1,1,X,Y,Z,lamb)+permO2(3,1,1,X,Y,Z,lamb)-
                  permO3(3,1,1,X,Y,Z,lamb)-permO4(3,1,1,X,Y,Z,lamb))
    
    term4=0.0138*(permE1(3,1,3,X,Y,Z,lamb)-permE2(3,1,3,X,Y,Z,lamb)+
                  permE3(3,1,3,X,Y,Z,lamb)-permE4(3,1,3,X,Y,Z,lamb)+
                  permO1(3,1,3,X,Y,Z,lamb)-permO2(3,1,3,X,Y,Z,lamb)+
                  permO3(3,1,3,X,Y,Z,lamb)+permO4(3,1,3,X,Y,Z,lamb))
    
    term5=0.0028*(permE1(1,1,5,X,Y,Z,lamb)+permE2(1,1,5,X,Y,Z,lamb)+
                  permE3(1,1,5,X,Y,Z,lamb)+permE4(1,1,5,X,Y,Z,lamb)+
                  permO1(1,1,5,X,Y,Z,lamb)+permO2(1,1,5,X,Y,Z,lamb)+
                  permO3(1,1,5,X,Y,Z,lamb)+permO4(1,1,5,X,Y,Z,lamb))
    
    term6=0.0021*(permE1(3,3,3,X,Y,Z,lamb)+permE2(3,3,3,X,Y,Z,lamb)+
                  permE3(3,3,3,X,Y,Z,lamb)+permE4(3,3,3,X,Y,Z,lamb)+
                  permO1(3,3,3,X,Y,Z,lamb)+permO2(3,3,3,X,Y,Z,lamb)+
                  permO3(3,3,3,X,Y,Z,lamb)+permO4(3,3,3,X,Y,Z,lamb))
    
    term7=0.0001*(permE1(3,1,5,X,Y,Z,lamb)+permE2(3,1,5,X,Y,Z,lamb)-
                  permE3(3,1,5,X,Y,Z,lamb)-permE4(3,1,5,X,Y,Z,lamb)+
                  permO1(3,1,5,X,Y,Z,lamb)+permO2(3,1,5,X,Y,Z,lamb)-
                  permO3(3,1,5,X,Y,Z,lamb)-permO4(3,1,5,X,Y,Z,lamb))
    
    surf_eq = term1+term2+term3-term4-term5-term6+term7
    
    #find vertices on the surface
    # print('generating surface vertices')
    vertices, simplices,normals, values = measure.marching_cubes_lewiner(surf_eq, level=0)

    #sort out the vertices and append them to the list of coordinates to write to file
    Xp,Yp,Zp = zip(*vertices)
    
    surface_points = np.array([Xp,Yp,Zp]).T
    
    return surface_points, normals


def file_reader(file, fit_atoms):

    pipeline = ov.io.import_file(file)
    pipeline.modifiers.append(ov.modifiers.WrapPeriodicImagesModifier())
    pipeline.modifiers.append(ov.modifiers.SelectTypeModifier(property = 'Particle Type', types = {'C5A'}))    
    data = pipeline.compute()
    b = np.where(data.particles.selection[:]>0)[0]    
    c = data.particles.positions[:][b]    
    pos = np.array(c)
    
    lp = data.cell[:][0][0]
    
    all_pos = data.particles.positions[:]
    
    names = []
    
    tprop = data.particles['Particle Type']
    for i, t in enumerate(tprop):
        names.append(tprop.type_by_id(t).name)    
    
    return lp, pos, all_pos, names

def main(file, ntrials=100, fit_atom = 'C5A'):
    
    print('starting', file)
    
    #read the file using the pdb_parsing file reader
    pdb_data = file_reader(file, fit_atom)
    
    lp = pdb_data[0]
    x = pdb_data[1][0:, 0]
    y = pdb_data[1][0:, 1]
    z = pdb_data[1][0:, 2]
    all_coords_init = pdb_data[2]

    #perform the fit by randomly guessing the mesophase and an initial translation of the terminal carbon coordinates
    result_dict = {}
    
    # print('starting randomised fit testing')
    
    for i in range (0,ntrials):
        # if i%(int(ntrials/10)) == 0:
            # print("\tfile: %s trial: %d/%d" %(file.split('/')[-1], i+1,ntrials))
        #print("trial: %d/%d" %(i+1,ntrials))
        #generate some random numbers between -pi and pi as initial guesses for the angle transformations
        a = (np.pi*np.random.random()) - (np.pi/2)
        b = (np.pi*np.random.random()) - (np.pi/2)
        c = (np.pi*np.random.random()) - (np.pi/2)

        #only bother looking for the fit of the QIID nodal surface        
        mesophase_guess = 'D'
        
        #give an initial guess of the translations and minimize sigma squared from there.
        x0=np.array([-x.mean(),-y.mean(),-z.mean(),a,b,c,lp])
        res=minimize(sigmasq,x0,args=(mesophase_guess,x,y,z))

        result_dict[i] = res.fun, a, b, c, mesophase_guess, res.x[6]
    
    #gather the results of the fit
    results = np.zeros(0)
    keys = np.zeros(0)
    for key in result_dict.keys():
        #so that badly fitted results won't appear
        if result_dict[key][0]>0.1:
            results = np.append(results, result_dict[key][0])
            keys = np.append(keys, key)
    
    #find the best result
    optimised_result = result_dict[keys[results.argmin()]]
    
    #best guess array for finding final coefficients.    
    xf = np.array([-x.mean(), -y.mean(), -z.mean(), optimised_result[1], optimised_result[2],optimised_result[3],lp])
    final_result= minimize(sigmasq, xf, args=(optimised_result[4],x,y,z))
    
    #find the final transformed coordinates and values to return       
    final_C_i, final_transformed_coords, final_matrix_coeffs = sigmasq(final_result.x, optimised_result[4],x,y,z, final =True)
    
    #need this condition otherwise something has gone wrong in the fit
    if final_result.fun<2 and np.abs(final_matrix_coeffs[6])<200:
        try:
            # print('lattice parameter was %f, fit result was %.3f, mesophase was %s' %(final_matrix_coeffs[6], final_result.fun, optimised_result[4]))

                        
            #find the transformation of all the coordinates in the input file and save them if required.

            # print('correcting data')
    
            #transform all the coordinates from the initial positions into the final ones.
            initial_transformed=np.zeros(0)
            for i in all_coords_init:
                initial_transformed = np.append(initial_transformed,translations(i[0],i[1],i[2],final_matrix_coeffs))
            initial_transformed = initial_transformed.reshape(all_coords_init.shape)
                    
            '''
            use marching cubes to generate points on the surface from the lattice parameter
            '''
            
            #get the size of the unit cell
            lamb = 1/final_matrix_coeffs[6]
            l = final_matrix_coeffs[6]/2
    
            #translate the data points to the origin.
            x_mean = np.mean(initial_transformed[0:,0])
            y_mean = np.mean(initial_transformed[0:,1])
            z_mean = np.mean(initial_transformed[0:,2])
            
            #find the points on the surface
            surface_points, normals = point_generator(x_mean,y_mean,z_mean,l,lamb)
            
            #translating points preparation
            a = np.ones_like(initial_transformed)
            corrected_data_points = np.ones_like(initial_transformed)
            
            #translate the data points by the difference between half the lattice parameter and the mean positions so that the origin of the box is at 0.
            corrected_data_points[0:,0] = initial_transformed[0:,0] - a[0:,0]*(x_mean-l)
            corrected_data_points[0:,1] = initial_transformed[0:,1] - a[0:,1]*(y_mean-l)
            corrected_data_points[0:,2] = initial_transformed[0:,2] - a[0:,2]*(z_mean-l)

            
            #translate the points on the surface to their correct spatial positions for curvature to be calculated.
            curvature_positions = np.zeros_like(surface_points)
            
            curvature_positions[0:,0] = surface_points[0:,0]+(x_mean-l)
            curvature_positions[0:,1] = surface_points[0:,1]+(y_mean-l)
            curvature_positions[0:,2] = surface_points[0:,2]+(z_mean-l)
        
            curvatures = np.zeros(0)
            

            for i in range(len(curvature_positions)):
                #for an update on progress
                #don't know how to actually make this regular, but it doesn't really matter
                if i/len(curvature_positions)*100 %10 <0.001:
                    prog = int(i/len(curvature_positions)*100)
                    # print('curvature calculation progress: %d%%' %prog)
        
                K = curvature(curvature_positions[i][0],curvature_positions[i][1],curvature_positions[i][2],lamb)
                curvatures = np.append(curvatures, K)
    
            '''
            '''
    
            curvature_bins = np.linspace(curvatures.min(), 0, 20)
            mids_curvature_bins = (curvature_bins[:-1]+curvature_bins[1:])/2
            
            inds = np.digitize(curvatures, curvature_bins)
            
            alpha = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k',
                     'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
                     'w', 'x', 'y', 'z']
            
            
            point_inds = np.zeros(0, dtype = int)
            opstr = []
            for i in range(1, inds.max()+1):
                find = np.where(inds == i)[0].astype(int)
                point_inds = np.append(point_inds, find)
                opstr+= [alpha[i]]*len(find)
                
            
            #by default the surface points are now distributed from 0 to the size of the box. 
            #Need to translate points by the difference in the means for correct output.
            x_mean_diff = np.mean(corrected_data_points[0:, 0]) - surface_points[0:,0].mean()
            y_mean_diff = np.mean(corrected_data_points[0:, 1]) - surface_points[0:,1].mean()
            z_mean_diff = np.mean(corrected_data_points[0:, 2]) - surface_points[0:,2].mean()
            
            t = np.ones_like(surface_points)
            corrected_surface_points = np.ones_like(surface_points)
            
            corrected_surface_points[0:,0] =  surface_points[0:,0] + t[0:,0]*x_mean_diff
            corrected_surface_points[0:,1] =  surface_points[0:,1] + t[0:,1]*y_mean_diff
            corrected_surface_points[0:,2] =  surface_points[0:,2] + t[0:,2]*z_mean_diff
            
            sp_final = corrected_surface_points[point_inds]
            
            all_points = np.concatenate((corrected_data_points,sp_final))
            
            #return the points to 0
            to_zero = np.zeros_like(all_points)
            
            to_zero[0:,0] = all_points[0:,0]-sp_final[0:,0].min()
            to_zero[0:,1] = all_points[0:,1]-sp_final[0:,1].min()
            to_zero[0:,2] = all_points[0:,2]-sp_final[0:,2].min()
            
            
            # return_names = pdb_data[3] + names
            return_names = pdb_data[3] + opstr
            
            name = file.split('/')[-1].split('.pdb')[0] + '_atoms_surface.atom'
            
            #write files out into a lammps format
            with open(name,'w') as f:
                atom.write(to_zero,np.array([[0,surface_points[0:,0].max()],
                                             [0,surface_points[0:,1].max()],
                                             [0,surface_points[0:,2].max()]]),f,return_names)
                    
                
            return_dict = {'final coefficients': final_C_i,
                           'objective function value': final_result.fun,
                           'lattice parameter': final_matrix_coeffs[6],
                           'mesophase': optimised_result[4],
                           'all points': to_zero,
                           'data points': corrected_data_points,
                           'surface points': sp_final,
                           'curvatures': curvatures[point_inds],
                           'point names': return_names,
                           'surface normals': normals[point_inds]}
            
            
            pickle.dump(return_dict, open(file.split('.pdb')[0]+'.p', 'wb'))
            
            print('finished', file)
            
            #return the fitting coefficients, the fit value, the final matrix coefficients, and the mesophase
            #the matrix coefficients are ordered as x,y,z,3 rotations, lattice parameter.
            return return_dict
        
        except UnboundLocalError:
            return 0    
    else:
        print('lattice parameter was %f, fit result was %.3f, mesophase was %s' %(final_matrix_coeffs[6], final_result.fun, optimised_result[4]))
        print('fit unsuccessful')
        return 0
    
    
    
if __name__ == '__main__':
    
    files = glob.glob('*.pdb')
    
    k = len(files)/14
    
    if k < 1:
        csize = 1
    else:
        csize = int(k)
    
    with get_context("spawn").Pool(processes = 14) as pool:
        pool.map(main, files, chunksize = csize)
    
    
    
    
