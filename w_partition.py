#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

# from matplotlib import pyplot as plt

import glob
import sys
import pickle
from natsort import natsorted

from ovito.io import import_file
from ovito.modifiers import ConstructSurfaceModifier, SelectTypeModifier, ClearSelectionModifier, WrapPeriodicImagesModifier

from itertools import product
from multiprocessing import get_context

def bd_position_testing(f, bd_types_init):

    bd_types = set(bd_types_init)
    
    loops = np.linspace(6, 15, 15)

    #open a file
    pipeline = import_file(f)

    pipeline.modifiers.append(WrapPeriodicImagesModifier())

    #select all lipid molecules
    pipeline.modifiers.append(SelectTypeModifier(property = 'Particle Type', types = {'WN'}))

    #construct a surface from the water molecules
    mod1 = ConstructSurfaceModifier(radius = 10,only_selected = True,smoothing_level =8)

    #and append it to the pipeline
    pipeline.modifiers.append(mod1)

    props = np.zeros(0)


    for i in range(len(loops)):
        try:
            # print('loops: %d/%d' %(i+1, len(loops)))
            mod1.radius = loops[i]
            
            #compute the pipeline
            data = pipeline.compute()
            print(f, data)

            #get the surface mesh
            mesh = data.surfaces['surface']
            #clear the pipeline of the selected water molecules
            pipeline.modifiers.append(ClearSelectionModifier())
            #now select the central butanediol molecule
            pipeline.modifiers.append(SelectTypeModifier(property = 'Particle Type', types = bd_types))
            #recompute the pipeline
            data1 = pipeline.compute()
            #find the particle indices of the selected butanediol particles
            sel = np.where(data1.particles.selection[:]>0)[0]
            #now find the positions of the butanediol particles
            BD_position = data1.particles.positions[:][sel]
            
            inside_surface = np.zeros(0)
            outside_surface = np.zeros(0)
            
            #for every position of a butanediol...
            for j in BD_position:
                #test whether the position is inside the water surface or not.
                pos = mesh.locate_point(j)
                print(pos)
                #keep count of whether it is or not.
                try:
                    if pos == -1:
                        outside_surface = np.append(outside_surface, pos)              
                    elif pos >= 0:
                        inside_surface = np.append(inside_surface, pos)
                except TypeError:
                    pass
            
            tot = BD_position.shape[0]
            
            #find the proportion of BD molecules that are inside the W surface
            p = len(inside_surface)/tot

            props = np.append(props, p)

        except RuntimeError:
            pass
    

    fname = f.split('.pdb')[0]+'_bd_pos_analysis.p'
    pickle.dump(props, open(fname, 'wb'))
    

    # try:
    #     figname = f.split('.pdb')[0] + '_w bd proportions.png'
    #     fig, ax = plt.subplots(1,1)
    #     ax.scatter(loops, props, label = 'inside')
    #     ax.set_xlabel('Probe Sphere radius')
    #     ax.set_ylabel('Proportion of butanediol inside water surface')

    #     ax.axhline(np.mean(props),c='r')
    #     ax.text(ax.get_xlim()[0]+ 0.5, ax.get_ylim()[1]*0.75,'mean = %.3f' %np.mean(props), c = 'r')

    #     fig.savefig(figname, bbox_inches = 'tight')

    #     #return the mean and the std of the proportion inside the surface.
    #     return props
    
    # except ValueError as e1:
    #     print(e1)
    #     pass
    
    return 0

if __name__ == '__main__':

    files = natsorted(glob.glob('*.pdb'))
    
    bd_types_init = sys.argv[1:]
    
    if len(bd_types_init) == 0:
        print('you must indicate which bd beads to find the positions of!')

    else:    
        
        print('butanediol beads are: ', bd_types_init)          
        
        
        props = np.zeros(0)
        errs = np.zeros(0)
        times = np.zeros(0)
        
        paramlist = []
        
        for i in files:
            paramlist.append([i,bd_types_init])
        
        
        
        
            
        k = len(paramlist)/14
            
        if k < 1:
            csize = 1
        else:
            csize = int(k)
        
        print(paramlist, csize)
        
        with get_context("spawn").Pool(processes = 14) as pool:
            pool.starmap(bd_position_testing, paramlist, chunksize = csize)
    
            
            
        
        
        
        # data_dict = {}
        
        # for k in range(len(files)):

        #     print('progress: %d/%d' %(k+1, len(files)))
            
        #     # plt.clf()
            
        #     inside_proportions = bd_position_testing(files[k], bd_types, loops)

        #     if type(inside_proportions)!= int:
        #         t = files[k].split('md')[1].split('-ns.pdb')[0]
                
        #         data_dict[int(t)] = inside_proportions
                
        #         times = np.append(times, int(t))
        #         props = np.append(props, np.mean(inside_proportions))
        #         errs = np.append(errs, np.std(inside_proportions))

        # pickle.dump((bd_types, loops, data_dict), open('W_BDTypes_times_proportions_error.p', 'wb'))







        # fig, ax = plt.subplots(1,1)

        # ax.errorbar(times, props, yerr = errs, linestyle = '', marker = '.', markersize =12)
        # ax.set_xlabel('Simulation Time (ns)')
        # ax.set_ylabel('Proportion of butanediol beads in water region')
        # ax.axhline(np.mean(props),c='r', linestyle='--')

        # #NB: this std is of the values only - should be weighted by the std of each data point somehow.
        # ax.set_title('mean proportion = %.3f$\pm$%.3f' %(np.average(props, weights = errs), np.std(props)))

        # fig.savefig('W_BD_occupancy.png', bbox_inches = 'tight')






















