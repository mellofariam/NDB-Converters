#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.gro to Nucleome Data Bank format .ndb
"""

__author__ = "Vinícius Contessoto / Matheus Mello"
__date__   = "Nov/2019"

################################################################
# 
# Gromacs file *.gro to Nucleome Data Bank format .ndb
#
# usage:
#  ./gro2ndb.py -f file.gro -n name_NDB_file
#
################################################################

import re
import time
import argparse
import linecache
import numpy as np
from   datetime  import date
from   itertools import islice

parser = argparse.ArgumentParser(description='Converting from *.gro to *.ndb file')
parser.add_argument('-f', metavar='input-file-Grofile-frames',help='grofile with frames',type=argparse.FileType('rt'))
parser.add_argument('-res', action='store', default=50000, dest='arg_res', help='Resolution for each simulation bead')
parser.add_argument('-chro', nargs='+', action='store', dest='arg_chro', help='Chromosome number', required=False, type=int)
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
master_string  = "{0:6s} {1:8d} {2:6d} {3:6d} {4:10d}"   # MASTER BEADS TER LOOPS RES

file_gro   = arguments.f
file_loops = arguments.loops
res        = np.int(arguments.arg_res)
chroID     = np.str(arguments.arg_chroID)
sigma      = np.float(arguments.arg_sigma)
scale      = np.float(arguments.arg_scale)

try:
    chro   = arguments.arg_chro
except:
    pass
    
b_time     = time.time()

# [ Getting beads number ]

all_beads  = int(linecache.getline(file_gro.name, 2).strip('\n'))
testPL     = linecache.getline(file_gro.name, all_beads + 2).strip('\n')

if 'PL' in testPL.split()[1]:
    beads = all_beads - 1
else:
    beads = all_beads

print('Number of beads:', beads)

# [ Writing header ]

ndbf = open(arguments.arg_name+'.ndb', "w+")

today = date.today()

ndbf.write(header_string.format('HEADER','NDB File for Chromosome {:}'.format(chro),today.strftime("%b %y"),''))
ndbf.write("\n")
ndbf.write(title_string.format('TITLE ','  ','The Nucleome Data Bank: Web-based Resources to'))
ndbf.write("\n")
ndbf.write(title_string.format('TITLE ',' 2',' Simulate and Analyze the Three-Dimensional Genome'))
ndbf.write("\n")
ndbf.write(expdata_string.format('EXPDTA','  ','Cell Line A549 @ 50k bp resolution'))
ndbf.write("\n")
ndbf.write(expdata_string.format('EXPDTA','  ','Simulation - MEGABASE - MiChroM'))
ndbf.write("\n")
ndbf.write(author_string.format('AUTHOR','  ','Contessoto et. al. - 2019'))
ndbf.write("\n")

# [ Creating SEQCHR ]

def chunk_list(lista, size):
    lista = iter(lista)
    return iter(lambda: tuple(islice(lista, size)), ())

SEQCHROM = []
bead     = []

chain    = 0

for num, line in enumerate(file_gro):
    if num > 1 and num < beads + 2:
        info = line.split()

        if "ChrA" in info[0] or "ChrB" in info[0] or "ChrU" in info[0]:
        
            temp  = re.findall(r'\d+', info[0])
            index = [ int(x) for x in temp ][0]  # Getting index number for 1st gro field

            if index == 1:
                chain += 1
                SEQCHROM.append([])
                bead.append(0)

            SEQCHROM[chain-1].append(Chrom_types_NDB[Chrom_types.index(str(info[1]))]) 
            bead[chain-1] += 1

for i in range(chain):

    seq_chunk_23 = list(chunk_list(list(SEQCHROM[i]), 23))
    
    for j in range(len(seq_chunk_23)):
        
        seq_str = " ".join(seq_chunk_23[j])
        try:
            ndbf.write(seqchr_string.format('SEQCHR', j+1, chroID + str(chro[i]), bead[i], seq_str))
        except:        
            ndbf.write(seqchr_string.format('SEQCHR', j+1, chroID + str(i + 1), bead[i], seq_str))
        ndbf.write("\n")

# [ Writing .ndb body ]

file_gro.seek(0)

frame_size = all_beads + 3

model = 1
chain = 0

index_c = 0

for num, line in enumerate(file_gro):

    info = line.split()

    if num % frame_size == 0:
        ndbf.write(model_string.format('MODEL ', model))
        ndbf.write("\n")

        index_c = 0
        numTer  = 0

    if "ChrA" in info[0] or "ChrB" in info[0] or "ChrU" in info[0]:
        
        temp  = re.findall(r'\d+', info[0])
        index = [ int(x) for x in temp ][0]  # Getting index number for 1st gro field

        if index == 1:

            if chain != 0:
                index_c    += 1
                try:
                    ndbf.write(ter_string.format('TER   ', index_c, subtype_ndb, chroID + str(chro[chain-1])))
                except:
                    ndbf.write(ter_string.format('TER   ', index_c, subtype_ndb, chroID + str(chain)))
                ndbf.write("\n")

                numTer += 1

            chain += 1

        X = np.float(info[3])
        Y = np.float(info[4])
        Z = np.float(info[5])

        index_c    += 1
        subtype_gro = str(info[1])

        subtype_ndb = Chrom_types_NDB[Chrom_types.index(str(info[1]))]
        
        start = np.int((index-1) * res)+1
        end   = np.int( index * res )
        
        try:
            ndbf.write(ndb_string.format('CHROM ', index_c, subtype_ndb, " ", chroID + str(chro[chain-1]), index, X, Y, Z, start, end, sigma)) 
        except:
            ndbf.write(ndb_string.format('CHROM ', index_c, subtype_ndb, " ", chroID + str(chain), index, X, Y, Z, start, end, sigma)) 
        ndbf.write("\n")

        if np.int(info[2]) == beads:
            index_c    += 1

            try:
                ndbf.write(ter_string.format('TER   ', index_c, subtype_ndb, chroID + str(chro[chain-1])))
            except:
                ndbf.write(ter_string.format('TER   ', index_c, subtype_ndb, chroID + str(chain)))
            ndbf.write("\n")

            numTer += 1

            model += 1
            chain = 0
            
            if model % 100 == 0:
                print("Traj frame:", model)
                e_time = time.time()
                elapsed = e_time - b_time
                print("Ran in %.3f sec" % elapsed)

    if num % frame_size == frame_size - 1:
        ndbf.write('ENDMDL')
        ndbf.write("\n")

print("Last traj frame:", model-1)

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

ndbf.write(master_string.format('MASTER', beads, numTer, mloops, res))
ndbf.write("\n")
ndbf.write('END')
ndbf.write("\n")

ndbf.close()

e_time  = time.time()
elapsed = e_time - b_time

print("Ran in %.3f sec" % elapsed)
