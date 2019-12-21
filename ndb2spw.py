#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.ndb to *.spw format
"""

__author__ = "Vinícius Contessoto / Matheus Mello"
__date__   = "Nov/2019"

################################################################
# 
# Trajectories file *.ndb to Nucleome Data Bank format .spw
#
# usage:
#  ./ndb2spw.py -f file.ndb -n name_SPW_file
#
################################################################

import re
import time
import argparse
import numpy as np

parser = argparse.ArgumentParser(description='Converting *.ndb to *.spw format')
parser.add_argument('-f', metavar='input-file-Grofile-frames',help='ndb file with frames',type=argparse.FileType('rt'))
parser.add_argument('-n', action='store', default='chromatin', dest='arg_name',  help='Name of output file')
parser.add_argument('-t', action='store', dest='header_name',  help='Value of name field at spw header')

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
model          = "trace {0:4d}"
atom           = "ATOM  {0:5d} {1:^4s}{2:1s}{3:3s} {4:1s}{5:4d}{6:1s}   {7:8.3f}{8:8.3f}{9:8.3f}{10:6.2f}{11:6.2f}          {12:>2s}{13:2s}"
ter            = "TER   {0:5d}      {1:3s} {2:1s}{3:4d}{4:1s}"

file_ndb = arguments.f
name     = arguments.arg_name

hname    = arguments.header_name

spwf     = open(name + '.spw', "w+")

print('Converting file...')

loop = 0

spwf.write('##format=sw1 ')

try:
    spwf.write('name={:}'.format(hname))
except:
    spwf.write('name={:} '.format(name))

first = True
for line in file_ndb:
    
    entry = line[0:6]

    info = line.split()

    if 'ASMBLY' in entry:
        spwf.write('genome=' + line[6:])

    if 'MODEL' in entry:
        if first:
            spwf.write('\nchromosome	start	end	x	y	z\n')
            first = False
        spwf.write('trace {:}\n'.format(int(line[6:].replace(' ',''))-1))
    
    elif 'CHROM' in entry:

        temp = re.findall(r'\d+', info[3])
        chro = [ int(x) for x in temp ][0]  # Getting chro number from 1st spw field

        spwf.write('chr{0:} {1:} {2:} {3:} {4:} {5:}\n'.format(chro, info[8], info[9], info[5], info[6], info[7]))

    # [ Loops file ]

    elif 'LOOPS' in entry:
        if loop == 0:
            loops = open(name + '.loops', 'w+')
        
        loops.write('{0:d} {1:d}\n'.format(int(info[1]), int(info[2])))
        loop += 1
        
spwf.write('END\n')
spwf.close()

print('Finished!')

if not loop == 0:
    loops.close
    print('Generated files: {:} and {:}'.format(name + '.spw', name + '.loops'))
else:
    print('Generated files: {:}'.format(name + '.spw'))

e_time = time.time()
elapsed = e_time - b_time
print('Ran in %.3f sec' % elapsed)
