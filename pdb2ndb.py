#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.pdb to Nucleome Data Bank format .ndb
"""

__author__ = "Vinícius Contessoto / Matheus Mello"
__date__   = "Nov/2019"

################################################################
# 
# Trajectories file *.pdb to Nucleome Data Bank format .ndb
#
# usage:
#  ./pdb2ndb.py -f file.pdb -n name_NDB_file
#
################################################################


import time
import argparse
import numpy as np
from   itertools     import islice

parser = argparse.ArgumentParser(description='Chromatin Gromacs simulation result plot')
parser.add_argument('-f', metavar='input-file-Grofile-frames',help='grofile with frames',type=argparse.FileType('rt'))
parser.add_argument('-res', action='store', default=50000, dest='arg_res', help='Resolution for each simulation bead')
parser.add_argument('-chroID', action='store', default='C', dest='arg_chroID', help='Chain ID')
parser.add_argument('-n', action='store', default='chromatin', dest='arg_name',  help='Name of output file')
parser.add_argument('-loops', metavar='input-file-Loops',help='Loops contact pair between i and j',type=argparse.FileType('rt'), required=False)
parser.add_argument('-sigma', action='store', default=0.000, dest='arg_sigma', help='Distance fluctuation')
parser.add_argument('-scale', action='store', default=1.000, dest='arg_scale', help='Distance scale')

try:
    arguments = parser.parse_args()
    print('################################################')
    print('Chosen file: {:}'.format(arguments.f.name))
    print('Resolution: {:}'.format(arguments.arg_res))
    print('Chrom ID: {:}'.format(arguments.arg_chroID))
    print('Distance scale: {:}'.format(arguments.arg_scale))

except IOError as msg:
    parser.error(str(msg))                    

Main_chrom      = ['ChrA','ChrB','ChrU'] # Type A B and Unknow
Chrom_types     = ['ZA','OA','FB','SB','TB','LB','UN']
Chrom_types_NDB = ['A1','A2','B1','B2','B3','B4','UN']

##################################################################################################

b_time = time.time()

gro_string     = "{0:5d}{1:5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}"
gro_box_string = "{0:10.5f}{1:10.5f}{2:10.5f}"
ndb_string     = "{0:6s} {1:8d} {2:2s} {3:6s} {4:4s} {5:8d} {6:8.3f} {7:8.3f} {8:8.3f} {9:10d} {10:10d} {11:8.3f}"
pdb_string     = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
header_string  = "{0:6s}    {1:40s}{2:9s}   {3:4s}"
title_string   = "{0:6s}  {1:2s}{2:80s}"
author_string  = "{0:6s}  {1:2s}{2:79s}"
expdata_string = "{0:6s}  {1:2s}{2:79s}"
model_string   = "{0:6s}     {1:4d}"
seqchr_string  = "{0:6s} {1:3d} {2:2s} {3:5d}  {4:69s}" 
ter_string     = "{0:6s} {1:8d} {2:2s}        {3:2s}" 
loops_string   = "{0:6s}{1:6d} {2:6d}"
master_string  = "{0:6s} {1:8d} {2:6d} {3:10d} {4:>8s}"

file_pdb   = arguments.f
file_loops = arguments.loops
res        = np.int(arguments.arg_res)
chroID     = np.str(arguments.arg_chroID)
sigma      = np.float(arguments.arg_sigma)
scale      = np.float(arguments.arg_scale)

ndbf       = open(arguments.arg_name+'.ndb', "w+")

for line in file_pdb:
    info = line.split()

    if 'MODEL' in info[0] or 'ATOM' in info[0]:
        break

    if 'CRYST1' in info[0]:
        ndbf.write(' '.join(['REMARK   ', ' '.join(info[1:])]))
        ndbf.write('\n')
    else:
        ndbf.write(line)

# [ Creating SEQCHR ]

file_pdb.seek(0)

def chunk_list(lista, size):
    lista = iter(lista)
    return iter(lambda: tuple(islice(lista, size)), ())

SEQCHROM = []
bead     = []

chain = 0

for line in file_pdb:
  
    info = line[0:6]

    if 'ATOM' in info:
        if not 'PL' in line[12:16]:
            
            index = int(line[22:26].replace(" ", ""))  # Getting index number for 1st pdb field

            if index == 1:
                chain += 1
                SEQCHROM.append([])
                bead.append(0)
            
            try:
                SEQCHROM[chain-1].append(str(Chrom_types_NDB[np.int(np.where(np.char.find(Chrom_types, str(line[12:16].replace(" ", ""))) == 0)[0][0])]))
            except:
                SEQCHROM[chain-1].append('UN')
            
            bead[chain-1] += 1
    
    if 'TER' in info:
        break

for i in range(chain):

  seq_chunk_23 = list(chunk_list(list(SEQCHROM[i]), 23))
  
  for j in range(len(seq_chunk_23)):
      
    seq_str = " ".join(seq_chunk_23[j])
    ndbf.write(seqchr_string.format('SEQCHR', j+1, chroID + str(i + 1), bead[i], seq_str))
    ndbf.write("\n")

# [ Writing .ndb body ]

file_pdb.seek(0)

model = 1
chain = 0

isAtom = False

for line in file_pdb:
    
    info = line[0:6]

    if 'ATOM' in info:
        if not 'PL' in line[12:16]:
            
            index_c     = int(line[6:11].replace(" ", ""))

            subtype_pdb = line[12:16].replace(" ", "")
            
            try:
                subtype_ndb  = str(Chrom_types_NDB[np.int(np.where(np.char.find(Chrom_types, subtype_pdb) == 0)[0][0])])
            except:
                subtype_ndb  = 'UN'

            index = int(line[22:26].replace(" ", ""))  # Getting index number for 1st pdb field

            if index == 1:
                if chain == 0:
                    ndbf.write(model_string.format('MODEL ', int(model)))
                    ndbf.write("\n")
                chain += 1

            X = np.float(line[30:38].replace(" ", ""))
            Y = np.float(line[38:46].replace(" ", ""))
            Z = np.float(line[46:54].replace(" ", ""))

            start = np.int((index-1) * res)+1
            end   = np.int( index * res )
            
            ndbf.write(ndb_string.format('CHROM ', index_c, subtype_ndb, " ", chroID + str(chain), index, X, Y, Z, start, end, sigma)) # Aqui a gente escreve as coordenadas e os campos coloridos
            ndbf.write("\n")
        
        isAtom = True
    
    else:
        if isAtom:
            ndbf.write(ter_string.format('TER   ', index_c + 1, subtype_ndb, chroID + str(chain)))
            ndbf.write("\n")
            model += 1
            chain = 0
            
            ndbf.write('ENDMDL\n')
    
        isAtom = False


print("Last traj frame:", model-1)

mloops=0
if file_loops is not None:
  for line in file_loops:

    loop = line.split()

    i = np.int(loop[0])
    j = np.int(loop[1])
    ndbf.write(loops_string.format('LOOPS ', i, j))
    ndbf.write("\n")

    mloops += 1

ndbf.write(master_string.format('MASTER', sum(bead), mloops, res, str(scale)))
ndbf.write("\n")
ndbf.write('END')
ndbf.write("\n")

ndbf.close()

e_time = time.time()
elapsed = e_time - b_time
print("Ran in %.3f sec" % elapsed)
