#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.ndb to *.pdb format
"""

__author__ = "Vinícius Contessoto / Matheus Mello"
__date__   = "Nov/2019"

################################################################
# 
# Trajectories file *.ndb to Nucleome Data Bank format .pdb
#
# usage:
#  ./ndb2pdb.py -f file.ndb -n name_PDB_file
#
################################################################

import time
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Converting *.ndb to *.pdb format')
parser.add_argument('-f', metavar='input-file-Grofile-frames',help='ndb file with frames',type=argparse.FileType('rt'))
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


##################################################################################################

b_time = time.time()

# pdb formats              1         2         3
#                 123456789012345678901234567890
model          = "MODEL     {0:4d}"
atom           = "ATOM  {0:5d} {1:^4s}{2:1s}{3:3s} {4:1s}{5:4d}{6:1s}   {7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}{11:6.2f}          {12:>2s}{13:2s}"
ter            = "TER   {0:5d}      {1:3s} {2:1s}{3:4d}{4:1s}"

file_ndb = arguments.f
name     = arguments.arg_name
pdbf = open(name + '.pdb', "w+")

print('Converting file...')

title_options = ['HEADER','OBSLTE','TITLE ','SPLT  ','CAVEAT','COMPND','SOURCE','KEYWDS','EXPDTA','NUMMDL','MDLTYP','AUTHOR','REVDAT','SPRSDE','JRNL  ','REMARK']

loop = 0

for line in file_ndb:
    
    entry = line[0:6]

    info = line.split()

    # [ Title Section ]

    if entry in title_options:
        pdbf.write(line)
    elif 'ASMBLY' in entry:
        pdbf.write(' '.join(['REMARK   ', 'Assembly: ' + line[6:]]))
        pdbf.write('\n')
    elif 'LENGTH' in entry:
        pdbf.write(' '.join(['REMARK   ', 'Length scale of ' + line[6:]]))
        pdbf.write('\n')
    elif 'CYCLE' in entry:
        pdbf.write(' '.join(['REMARK   ', 'Cell cycle: ' + line[6:]]))
        pdbf.write('\n')

    # [ Primary Structure Section ]

    # [ Coodinate Section ]

    elif 'MODEL' in entry:
        pdbf.write(model.format(int(info[1])))
        pdbf.write('\n')

        inModel = True
    
    elif 'CHROM' in entry:

        subtype = line[16:18]

        residue = Res_types_PDB[Chrom_types_NDB.index(subtype)]

        pdbf.write(atom.format(int(line[7:15]), 'CA', '', residue, '', int(line[31:39]), '', float(line[40:48]), float(line[49:57]), float(line[58:66]), 1.00, 0.00, 'C', ''))
        pdbf.write('\n')

        atom_serial  = int(line[7:15])
        res_sequence = int(line[31:39])

    elif 'TER   ' == entry:

        pdbf.write(ter.format(int(line[7:15]), residue, '', res_sequence, ''))
        pdbf.write('\n')

    elif 'ENDMDL' in entry:
        pdbf.write('ENDMDL\n')

    # [ Loops file ]

    elif 'LOOPS' in entry:
        if loop == 0:
            loops = open(name + '.loops', 'w+')
        
        loops.write('{0:d} {1:d}\n'.format(int(info[1]), int(info[2])))
        loop += 1
        
pdbf.write('END\n')
pdbf.close()

print('Finished!')

if not loop == 0:
    loops.close
    print('Generated files: {:} and {:}'.format(name + '.pdb', name + '.loops'))
else:
    print('Generated files: {:}'.format(name + '.pdb'))

e_time = time.time()
elapsed = e_time - b_time
print('Ran in %.3f sec' % elapsed)
