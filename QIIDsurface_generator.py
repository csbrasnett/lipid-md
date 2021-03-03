#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

author: Chris Brasnett, University of Bristol, christopher.brasnett@bristol.ac.uk

"""

import numpy as np
from skimage import measure
from mdma import atom
from QIIDcurvature import main as curvature

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
    
    #term2=0.2260*(permE1(1,0,0,x,y,z,lamb)+permO1(1,0,0,x,y,z,lamb))
    
    #term3=0.0516*(permE1(1,1,1,x,y,z,lamb)+permO1(1,1,1,x,y,z,lamb))
    
    #term4=0.0196*(permE1(2,1,0,x,y,z,lamb)+permO1(2,1,0,x,y,z,lamb))
    
    #term5=0.0027*(permE1(3,0,0,x,y,z,lamb)+permO1(3,0,0,x,y,z,lamb))
    
    return term1#+term2-term3-term4-term5


def G_nodal_approx(x,y,z,lamb):
    
    term1=permE5(1,1,0,x,y,z,lamb)+permO5(1,1,0,x,y,z,lamb)
    
    # term2=0.3435*(permE4(1,1,0,x,y,z,lamb)+permO5(1,1,0,x,y,z,lamb))
    
    # term3=0.0202*(permE6(0,3,1,x,y,z,lamb)+permO6(0,3,1,x,y,z,lamb))
    
    # term4=0.0106*(permE7(2,2,2,x,y,z,lamb)+permO7(2,2,2,x,y,z,lamb))
    
    # term5=0.0298*(permE6(2,1,3,x,y,z,lamb)+permO6(2,1,3,x,y,z,lamb))
    
    # term6=0.0016*(permE8(3,0,3,x,y,z,lamb)+permO8(3,0,3,x,y,z,lamb))
    
    # term7=0.0011*(permE4(1,1,4,x,y,z,lamb)+permO5(1,1,4,x,y,z,lamb))
    
    return term1#+term2-term3-term4+term5+term6-term7
    

def generator(a_1,a_2,a_3,l,D = False, P = False, G = False):

    lamb = (1/l)/2

    #define the mgrid from which to generate the points
    X,Y,Z = np.mgrid[a_1-l:a_1+l:(2*l * 1j), 
                     a_2-l:a_2+l:(2*l * 1j), 
                     a_3-l:a_3+l:(2*l * 1j)]
    
    if D:
        surf_eq = D_nodal_approx(X,Y,Z, lamb)
    elif P:
        surf_eq = P_nodal_approx(X,Y,Z, lamb)
    elif G:
        surf_eq = G_nodal_approx(X,Y,Z, lamb)
    
    #find vertices on the surface
    vertices, simplices,normals, values = measure.marching_cubes_lewiner(surf_eq, level=0)
      
    #sort out the vertices and append them to the list of coordinates to write to file
    Xp,Yp,Zp = zip(*vertices)
    
    table = np.array([Xp,Yp,Zp]).T
    
    correct_positions = np.zeros_like(table)
    
    correct_positions[0:,0] = table[0:,0]+(a_1-l)
    correct_positions[0:,1] = table[0:,1]+(a_2-l)
    correct_positions[0:,2] = table[0:,2]+(a_3-l)
    
    if D==True:
        curvatures = np.zeros(0)
        for i in range(len(correct_positions)):
            #for an update on progress
            #don't know how to actually make this regular, but it doesn't really matter
            if i/len(correct_positions)*100 %10 <0.001:
                prog = int(i/len(correct_positions)*100)
                print('curvature calculation progress: %d%%' %prog)
            K = curvature(correct_positions[i][0],correct_positions[i][1],correct_positions[i][2],lamb)
            curvatures = np.append(curvatures, K)
        
        return correct_positions, curvatures
        
    else:
        return correct_positions



if __name__=='__main__':

    
    
    points, curvatures = generator(40,40,40,40, D = True)    
        
    curvature_bins = np.linspace(curvatures.min(), 0, 30)
    mids_curvature_bins = (curvature_bins[:-1]+curvature_bins[1:])/2
    
    inds = np.digitize(curvatures, curvature_bins)
    
    point_inds = np.zeros(0, dtype = int)
    opstr = []
    for i in range(1, inds.max()+1):
        find = np.where(inds == i)[0].astype(int)
        point_inds = np.append(point_inds, find)
        opstr+= [str(i)]*len(find)
    
    points_out = points[point_inds]
                
    with open('../D.atom', 'w') as f:
        atom.write(points_out, np.array([[points_out[0:,0].min(), points_out[0:,0].max()],
                                          [points_out[0:,1].min(), points_out[0:,1].max()],
                                          [points_out[0:,2].min(), points_out[0:,2].max()]]), f, opstr)
     


    # points = generator(25,25,25,50, P = True)    
    # points_out = points
        
                
    # with open('../temp1.atom', 'w') as f:
    #     atom.write(points_out, np.array([[points_out[0:,0].min(), points_out[0:,0].max()],
    #                                      [points_out[0:,1].min(), points_out[0:,1].max()],
    #                                      [points_out[0:,2].min(), points_out[0:,2].max()]]), f)
     


    
    
    







