#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 26 16:56:26 2019

@author: Chris
"""

import pandas as pd
import numpy as np
import re
import argparse


def reader(file):
    MOx = []
    MOy = []
    MOz = []
    MO_atomname = []
    MO_atomnumber = []
    MO_residuenumber = []
    
    
    Wx = []
    Wy = []
    Wz = []
    W_atomname = []
    W_atomnumber = []
    W_residuenumber = []
    
    
    with open(file, 'r') as f:
        for line in f:   
            if line[:7] == 'Lattice':
                # this regex finds the part of the line between quotes
                quotes = re.findall(r'"[^"^\u201d]*["\u201d]',line)
                
                vals_list = list(quotes[0][1:-1].split(' '))
                lp_arr = np.array(vals_list, dtype = float)
                
                # print(lp_arr, lp_arr.shape, type(lp_arr))
                # for i in lp_arr:
                #     print(i, type(i))
                
                lp_arr = lp_arr[lp_arr != 0]
                
            if len(line.split())==7:
                
                #these have to be converted back to nm for gro files
                x = float(line.split()[0])/10
                y = float(line.split()[1])/10
                z = float(line.split()[2])/10
                
                atom_name = line.split()[3]
                #atom_number = int(line.split()[4])
                residue_number = int(line.split()[5])
                
                if line.split()[6] == '1':
                    #residue_name = 'MAG99'
                    
                    MOx.append(x)
                    MOy.append(y)
                    MOz.append(z)
                    MO_atomname.append(atom_name)
    
                    MO_residuenumber.append(residue_number)
                    
                    
                elif line.split()[6] == '2':
                    #residue_name = 'W'
                    
                    Wx.append(x)
                    Wy.append(y)
                    Wz.append(z)
                    W_atomname.append(atom_name)
                    
                    W_residuenumber.append(residue_number)
                    
                    
    MAG99_list = ['MO']* len(MOx)
    W_list = ['CHOL']* len(Wx)
    
    MO_atomnumber = list(range(1, len(MOx)+1))
    W_atomnumber = list(range(len(MOx)+1, len(MOx)+1+len(Wx)))
    
    n_atoms = int(len(MO_atomnumber) + len(W_atomnumber))
    
    
    
    a = np.zeros(0)
    
    for i in range(1,int(len(MO_atomnumber)/7)+1):
        for j in range(7):
            a = np.append(a, i)

    
    b = np.zeros(0)
    
    for i in range(1,int(len(W_atomnumber)/9)+1):
        for j in range(9):
            b = np.append(b, i+a[-1])

    


    # b = list(range(int(a[-1]+1), int(a[-1]+1+len(W_atomnumber))))
    
    # print(len(a), len(MAG99_list), len(MO_atomname), len(MO_atomnumber), len(MOx), len(MOy), len(MOz))
    
    MO_df = pd.DataFrame({'Residue Number': a, 'MAG99': MAG99_list, 'Atom Name': MO_atomname,
                          'Atom Number': MO_atomnumber, 'x': MOx, 'y': MOy, 'z': MOz})
    

    # print(W_residuenumber,b)
    
    W_df = pd.DataFrame({'Residue Number': b, 'CHOL': W_list, 'Atom Name': W_atomname,
                          'Atom Number': W_atomnumber, 'x': Wx, 'y': Wy, 'z': Wz})
    
    data_dict = {'MO_df': MO_df,
                'W_df': W_df,
                'n_atoms': n_atoms,
                'lp': lp_arr}
    
    return data_dict
    

def writer(data_dict,file):
    
    MO_df = data_dict['MO_df']
    W_df = data_dict['W_df']
    n_atoms = data_dict['n_atoms']
    lp = data_dict['lp']/10
    
    lp_str = '%.3f %.3f %.3f' %(lp[0], lp[1], lp[2])
    
    fname = file.split('.xyz')[0]
    
    with open(fname+'.gro', 'w') as f1:
        
        f1.write('this is the title string\n')
    
        f1.write('%d\n' %n_atoms)
        
        for i in range(MO_df.shape[0]):
            
            a = MO_df.iloc[i]['Residue Number']
            b = MO_df.iloc[i]['MAG99']
            c = MO_df.iloc[i]['Atom Name']
            d = MO_df.iloc[i]['Atom Number']
            e = MO_df.iloc[i]['x']
            f = MO_df.iloc[i]['y']
            g = MO_df.iloc[i]['z']
            
            l = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(a,b,c,d,e,f,g)
                    
            f1.write(l)
            
        for j in range(W_df.shape[0]):
            
            a1 = W_df.iloc[j]['Residue Number']
            b1 = W_df.iloc[j]['CHOL']
            c1 = W_df.iloc[j]['Atom Name']
            d1 = W_df.iloc[j]['Atom Number']
            e1 = W_df.iloc[j]['x']
            f11 = W_df.iloc[j]['y']
            g1 = W_df.iloc[j]['z']
            
            l1 = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" %(a1,b1,c1,d1,e1,f11,g1)
                    
            f1.write(l1)
        
        f1.write(lp_str)


def file_reader():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-f', '--file', type = str, nargs = '+',  help = 'The file to convert from .xyz to .gro format')

    args = parser.parse_args()
    
    
    files = args.file
        
    return files 


if __name__ == '__main__':
    a = file_reader()
    
    for i in a:
        b = reader(i)
        writer(b,i)



    