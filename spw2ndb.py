#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.spw to Nucleome Data Bank format .ndb
"""

__author__ = "Vinícius Contessoto / Matheus Mello"
__date__   = "Nov/2019"

################################################################
# 
# SpaceWalk file *.spw to Nucleome Data Bank format .ndb
#
# usage:
#  ./spw2ndb.py -f file.spw -n name_NDB_file
#
################################################################

import re
import time
import argparse 
import numpy as np
from   datetime  import date
from   itertools import islice

parser = argparse.ArgumentParser(description='Converting from *.spw to *.ndb file')
parser.add_argument('-f', metavar='input-file-spwFile-frames',help='spw file with frames',type=argparse.FileType('rt'))
parser.add_argument('-res', action='store', default=50000, dest='arg_res', help='Resolution for each simulation bead')
parser.add_argument('-chroID', action='store', default='C', dest='arg_chroID', help='Chain ID')
parser.add_argument('-n', action='store', default='traj', dest='arg_name',  help='Name of output file')
parser.add_argument('-loops', metavar='input-file-Loops',help='Loops contact pair between i and j',type=argparse.FileType('rt'), required=False)
parser.add_argument('-sigma', action='store', default=0.000, dest='arg_sigma', help='Distance fluctuation')

try:
    arguments = parser.parse_args()
    print('################################################')
    print('Chosen file: {:}'.format(arguments.f.name))
    print('Resolution: {:}'.format(arguments.arg_res))
    print('Chrom ID: {:}'.format(arguments.arg_chroID))

except IOError as msg:
    parser.error(str(msg))                    

Main_chrom      = ['ChrA','ChrB','ChrU'] # Type A B and Unknow
Chrom_types     = ['ZA','OA','FB','SB','TB','LB','UN']
Chrom_types_NDB = ['A1','A2','B1','B2','B3','B4','UN']

##################################################################################################

b_time = time.time()

gro_string     = "{0:5d}{1:5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}"
gro_box_string = "{0:10.5f}{1:10.5f}{2:10.5f}"
ndb_string     = "{0:6s} {1:8d} {2:2s} {3:6s} {4:4s} {5:8d} {6:8.1f} {7:8.1f} {8:8.1f} {9:10d} {10:10d} {11:8.3f}"
pdb_string     = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
header_string  = "{0:6s}    {1:40s}{2:9s}   {3:4s}"
title_string   = "{0:6s}  {1:2s}{2:80s}"
author_string  = "{0:6s}  {1:2s}{2:79s}"
expdata_string = "{0:6s}  {1:2s}{2:79s}"
model_string   = "{0:6s}     {1:4d}"
seqchr_string  = "{0:6s} {1:3d} {2:2s} {3:5d}  {4:69s}" 
ter_string     = "{0:6s} {1:8d} {2:2s}        {3:2s}" 
loops_string   = "{0:6s} {1:8d} {2:8d}"
master_string  = "{0:6s} {1:8d} {2:6d} {3:6d} {4:10d}"   # MASTER BEADS TER LOOPS RES

file_spw   = arguments.f
file_loops = arguments.loops
res        = np.int(arguments.arg_res)
chroID     = np.str(arguments.arg_chroID)
sigma      = np.float(arguments.arg_sigma)

ndbf       = open(arguments.arg_name+'.ndb', "w+")

# [ Writing header ]

today = date.today()

ndbf.write(header_string.format('HEADER','File converted from .spw file',today.strftime("%b %d %y"),''))
ndbf.write("\n")

# [ Creating SEQCHR ]

def chunk_list(lista, size):
    lista = iter(lista)
    return iter(lambda: tuple(islice(lista, size)), ())

SEQCHROM = []

bead  = []
chain = 0

chro_nr = []

for line in file_spw:
  
    info = line.split()
    
    if 'chr' in info[0] and info[0] != 'chromosome':
        
        temp = re.findall(r'\d+', info[0])
        chro = [ int(x) for x in temp ][0]  # Getting chro number from 1st spw field
        
        if not chro in chro_nr:
            chro_nr.append(chro)
            chain += 1
            SEQCHROM.append([])
            bead.append(0)

        bead[chain-1] += 1
        SEQCHROM[chain-1].append('UN')

    if info[0] == 'trace' and info[1] == '1':
        break

for i in range(chain):

  seq_chunk_23 = list(chunk_list(list(SEQCHROM[i]), 23))
  
  for j in range(len(seq_chunk_23)):
      
    seq_str = " ".join(seq_chunk_23[j])
    ndbf.write(seqchr_string.format('SEQCHR', j+1, chroID + str(chro_nr[i]), bead[i], seq_str))
    ndbf.write("\n")

# [ Writing ndb file body ]

file_spw.seek(0)

chain   = 0
index_c = 0
model   = 0

chro_nr = []

for line in file_spw:
    
    info = line.split()

    if 'trace' in info[0]:
        model += 1 

        ndbf.write(model_string.format('MODEL ', model))
        ndbf.write("\n")

        chain = 0

    elif 'chr' in info[0] and info[0] != 'chromosome':

        temp = re.findall(r'\d+', info[0])
        chro = [ int(x) for x in temp ][0]  # Getting chro number from 1st spw field

        if not chro in chro_nr:
            if chain != 0:
                index_c += 1
                ndbf.write(ter_string.format('TER   ', index_c, subtype_ndb, chroID + str(chro_nr[len(chro_nr)-1])))
                ndbf.write("\n")

            chro_nr.append(chro)
            chain += 1
            index  = 0

        index_c += 1

        subtype_ndb = 'UN'

        index += 1

        X = np.float(info[3])
        Y = np.float(info[4])
        Z = np.float(info[5])

        start = int(info[1])
        end   = int(info[2])

        ndbf.write(ndb_string.format('CHROM ', index_c, subtype_ndb, " ", chroID + str(chro), index, X, Y, Z, start, end, sigma)) # Aqui a gente escreve as coordenadas e os campos coloridos
        ndbf.write("\n")

        if index_c == sum(bead):
            index_c += 1

            ndbf.write(ter_string.format('TER   ', index_c, subtype_ndb, chroID + str(chro)))
            ndbf.write("\n")

            index_c = 0
            chro_nr = []

            ndbf.write('ENDMDL\n')

print("Last traj frame:", model)

# [ Writing loops ]

mloops = 0
if file_loops is not None:
  for line in file_loops:

    loop = line.split()

    i = np.int(loop[0])
    j = np.int(loop[1])
    ndbf.write(loops_string.format('LOOPS ', i, j))
    ndbf.write("\n")

    mloops += 1

ndbf.write(master_string.format('MASTER', sum(bead), chain, mloops, res))
ndbf.write("\n")
ndbf.write('END')
ndbf.write("\n")

ndbf.close()

e_time = time.time()
elapsed = e_time - b_time
print("Ran in %.3f sec" % elapsed)
