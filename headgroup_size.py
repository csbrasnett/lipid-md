#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import ovito as ov
import glob
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy import optimize
import pickle
from itertools import product
from multiprocessing import get_context

def func(x, phi_l):
    a = (4/3)*np.pi*-2
    b = 2*1.919
    q = (a*(x**3)) + (b*x) - phi_l
    return q
    
    
    
def main(file):
    
    # hgs = np.zeros(0)
    # time = np.zeros(0)
    # vols = np.zeros(0)

    try:
        
        pipeline = ov.io.import_file(file)

        #select W particles in the system.            
        pipeline.modifiers.append(ov.modifiers.SelectTypeModifier(property = 'Particle Type',
                                                                  types = {'W'}))
        
        
        solid_vols = np.zeros(0)    
        cell_vols = np.zeros(0)    
        sa = np.zeros(0)
        calc_area = np.zeros(0)
        
        mod1 = ov.modifiers.ConstructSurfaceModifier(only_selected = True, 
                                                      radius = 1, 
                                                      smoothing_level = 8, 
                                                      identify_regions = True)


        pipeline.modifiers.append(mod1)

        loops = np.linspace(10,16,15)

        for i in loops:
            mod1.radius = np.round(i, decimals = 1)
            
            data = pipeline.compute()
            
            area = data.attributes['ConstructSurfaceMesh.surface_area']
            
            solid_volume = data.attributes['ConstructSurfaceMesh.filled_volume']
            cell_volume = data.attributes['ConstructSurfaceMesh.cell_volume']
            
            fraction = solid_volume/cell_volume
            
            tprop = data.particles['Particle Type']

            #get the c5a id
            c5a_id = tprop.type_by_name('C5A').id
            c2_id = tprop.type_by_name('C2').id
            
            
            #count the number of terminal carbons
            n_lip = np.count_nonzero(tprop == c5a_id) + np.count_nonzero(tprop == c2_id) 
            
            phi_l = 1-fraction
            
            #need a in nm not Angstroms
            a = data.cell.matrix[0,0]
            
            sigma = 1.919
            chi = -2
            
            root = optimize.fsolve(func, x0 = [0], args = phi_l)

            l = root[0]*a
            A_L = (sigma*(a**2))+((2*np.pi*chi)*(l**2))
            
            a_0 = ((2*A_L)/(n_lip))
                    
            calc_area = np.append(calc_area, a_0)
            solid_vols = np.append(solid_vols, solid_volume)
            cell_vols = np.append(cell_vols, cell_volume)
            sa = np.append(sa, area/n_lip)


        # print('MAKING PLOT NOW')
        # fig, (ax0,ax1) = plt.subplots(2,1,sharex = True)

        # ax0.scatter(loops, v)
        # ax1.scatter(loops, sa, label = 'measured')
        # ax1.scatter(loops, calc_area, label = 'calculated')
        # ax1.legend()

        # ax0.set_ylabel('Surface Volume\nFraction in Unit Cell')
        # ax1.set_ylabel('Surface Area\nper molecule (Å$^2$)')

        # ax0.axhline(v.mean())
        # ax1.axhline(sa.mean())
        # ax1.axhline(calc_area.mean())
    
        # ax0.text(loops[1],v.mean(),'Mean = '+str(v.mean())[:4])
        # ax1.text(loops[1],sa.mean(),'Mean = '+str(sa.mean())[:4] +' Å')
        # ax1.text(loops[1],calc_area.mean(),'Mean = '+str(calc_area.mean())[:4] +' Å')
    
        # ax1.set_xlabel('Probe Sphere Radius')    

        # fig.subplots_adjust(hspace=0.1)
    
        # name = files[f].split('.pdb')[0] + ' headgroup analysis.png'
        
        # fig.savefig(name, dpi =200)
        
        
        
        d = {'Solid Volume': solid_vols,
             'Cell Volume': cell_vols,
             'Surface Area per Lipid': sa,
             'Calculated Area per Lipid': calc_area}
        
        
        dname = file.split('.pdb')[0] + '_headgroup_analysis.p'
        pickle.dump(d, open(dname, 'wb'))
        
        # t = file.split('md')[1].split('-')[0]

        # hgs = np.append(hgs, sa.mean())
        # vols = np.append(vols, v.mean())
        # time = np.append(time, int(t))
        
    except RuntimeError:
        print('error!', file)
        pass
    
    
    # fig1, ax2 = plt.subplots(1,1)
    # ax2.scatter(time/100, hgs)
    # ax2.set_xlabel('Simulation Time ($\mu$s)')
    # ax2.set_ylabel('Mean Head Group Area (Å$^{2}$)')
    # ax2.axhline(hgs.mean())
    # ax2.set_xlim(0,time.max()/100+0.1)
    # ax2.text(0,hgs.mean(), 'mean = %.3f, std = %.3f' %(hgs.mean(), hgs.std()))
    # fig1.savefig(folder+'/head group areas.png', dpi =200)
    
    # fig2, ax3 = plt.subplots(1,1)
    # ax3.scatter(time/100, vols)
    # ax3.set_xlabel('Simulation Time ($\mu$s)')
    # ax3.set_ylabel('Fractional Volume of Surface')
    # ax3.axhline(vols.mean())
    # ax3.set_xlim(0,time.max()/100+0.1)
    # ax3.text(0,vols.mean(), 'mean = %.3f, std = %.3f' %(vols.mean(), vols.std()))
    # fig2.savefig(folder+'/volumes.png', dpi =200)
    
    


if __name__ == '__main__':
    
    folder = os.getcwd()
    
    files = glob.glob(folder+'/*.pdb')

    paramlist = list(product(files))
    


    k = len(paramlist)/14
        
    if k < 1:
        csize = 1
    else:
        csize = int(k)
    
    print(paramlist, csize)
    
    with get_context("spawn").Pool(processes = 14) as pool:
        pool.starmap(main, paramlist, chunksize = csize)


    