#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.ndb to *.cndb format
"""

__author__ = "Vinícius Contessoto / Matheus Mello / Antonio B. Oliveira Junior"
__date__   = "Sep/2020"

################################################################
# 
# Trajectories file *.ndb to Compacted Nucleome Data Bank format .cndb
#
# usage:
#  ./ndb2cndb.py -f file.ndb -n name_CNDB_file
#
################################################################

import time
import argparse
import numpy as np
import h5py

parser = argparse.ArgumentParser(description='Converting *.ndb to *.cndb format')
parser.add_argument('-f', metavar='input-file-ndb-frames',help='ndb file with frames',type=argparse.FileType('rt'))
parser.add_argument('-n', action='store', default='chromatin', dest='arg_name',  help='Name of output file')

try:
    arguments = parser.parse_args()
    print('################################################')
    print('Chosen file: {:}'.format(arguments.f.name))

except IOError as msg:
    parser.error(str(msg))                    

Main_chrom      = ['ChrA','ChrB','ChrU'] # Type A B and Unknow
Chrom_types     = ['ZA','OA','FB','SB','TB','LB','UN']
Chrom_types_NDB = ['A1','A2','B1','B2','B3','B4','UN']
Res_types_PDB   = ['ASP', 'GLU', 'ARG', 'LYS', 'HIS', 'HIS', 'GLY']
Type_conversion = {'A1': 0,'A2' : 1,'B1' : 2,'B2' : 3,'B3' : 4,'B4' : 5,'UN' : 6}

##################################################################################################

b_time = time.time()

# pdb formats              1         2         3
#                 123456789012345678901234567890
model          = "MODEL     {0:4d}"
atom           = "ATOM  {0:5d} {1:^4s}{2:1s}{3:3s} {4:1s}{5:4d}{6:1s}   {7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}{11:6.2f}          {12:>2s}{13:2s}"
ter            = "TER   {0:5d}      {1:3s} {2:1s}{3:4d}{4:1s}"

file_ndb = arguments.f
name     = arguments.arg_name

cndbf = h5py.File(name + '.cndb', 'w')

print('Converting file...')

title_options = ['HEADER','OBSLTE','TITLE ','SPLT  ','CAVEAT','COMPND','SOURCE','KEYWDS','EXPDTA','NUMMDL','MDLTYP','AUTHOR','REVDAT','SPRSDE','JRNL  ','REMARK']

loop = 0
types = []
types_bool = True
loop_list = []
x = []
y = [] 
z = []

frame = 0

for line in file_ndb:
    
    entry = line[0:6]

    info = line.split()

    # [ Primary Structure Section ]

    # [ Coodinate Section ]

    if 'MODEL' in entry:
        frame += 1

        inModel = True
    
    elif 'CHROM' in entry:

        subtype = line[16:18]
        #print(subtype, line)

        types.append(subtype)
        x.append(float(line[40:48]))
        y.append(float(line[49:57]))
        z.append(float(line[58:66]))

    elif 'ENDMDL' in entry:
        if types_bool:
            typelist = [Type_conversion[x] for x in types]
            cndbf['types'] = typelist
            types_bool = False
        
        positions = np.vstack([x,y,z]).T
        cndbf[str(frame)] = positions
        x = []
        y = []
        z = []

    # [ Loops file ]

    elif 'LOOPS' in entry:
        loop_list.append([int(info[1]), int(info[2])])
        loop += 1
        
if loop > 0:
    cndbf['loops'] = loop_list

cndbf.close()
print('Finished!')


e_time = time.time()
elapsed = e_time - b_time
print('Ran in %.3f sec' % elapsed)
