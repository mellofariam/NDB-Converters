#!/usr/bin/env python
#coding: utf8 

__description__ = \
"""
Converting *.csv to Nucleome Data Bank format .ndb
"""

__author__ = "Vinícius Contessoto / Matheus Mello"
__date__   = "Nov/2019"

################################################################
# 
# Trajectories file *.csv to Nucleome Data Bank format .ndb
#
# usage:
#  ./csv2ndb.py -f file.csv -n name_NDB_file
#
################################################################

# Ps.: this script was written for Bintu et al. 2018 csv files

import time
import argparse
import numpy as np
from   datetime  import date
from   itertools import islice

parser = argparse.ArgumentParser(description='Chromatin Gromacs simulation result plot')
parser.add_argument('-f', metavar='input-file-csvFile-frames',help='csv file with frames',type=argparse.FileType('rt'))
parser.add_argument('-res', action='store', default=30000, dest='arg_res', help='Resolution for each simulation bead')
parser.add_argument('-chro', action='store', dest='arg_chro', help='Chro Number')
parser.add_argument('-chroID', action='store', default='C', dest='arg_chroID', help='Chain ID')
parser.add_argument('-n', action='store', default='chromatin', dest='arg_name',  help='Name of output file')
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
# title_string   = "{0:6s}  {1:>2s}{2:80s}"
title_string   = "{0:6s}  {1:>2s}"
author_string  = "{0:6s}  {1:2s}{2:79s}"
expdata_string = "{0:6s}  {1:2s}{2:79s}"
cycle_string   = "{0:6s}  {1:14s} {2:57s}"
model_string   = "{0:6s}     {1:4d}"
seqchr_string  = "{0:6s} {1:3d} {2:2s} {3:5d}  {4:69s}" 
ter_string     = "{0:6s} {1:8d} {2:2s}        {3:2s}" 
loops_string   = "{0:6s} {1:8d} {2:8d}"
master_string  = "{0:6s} {1:8d} {2:6d} {3:6d} {4:10d}"   # MASTER BEADS TER LOOPS RES

file_csv   = arguments.f
file_loops = arguments.loops
res        = np.int(arguments.arg_res)
chro       = np.int(arguments.arg_chro)
chroID     = np.str(arguments.arg_chroID)
sigma      = np.float(arguments.arg_sigma)

ndbf       = open(arguments.arg_name+'.ndb', "w+")

# [ Writing header ]

today = date.today()

ndbf.write(header_string.format('HEADER','NDB File for {:}'.format(file_csv.name.split('.')[0]),'',''))
ndbf.write("\n")

def chunkstring(string, length):
    return (string[0+i:length+i] for i in range(0, len(string), length))    

len_title = 0

ndbf.write(title_string.format('TITLE ','  '))
t = 1

for line in file_csv:
    if 'Chromosome' in line.split()[0]:
        break

    title = line.split()

    for i in range(len(title)):

        if len_title + len(title[i]) > 70:
            t += 1
            ndbf.write('\n')
            ndbf.write(title_string.format('TITLE ',str(t)))
            len_title = 1
        
        if t == 1 and i == 0:
            len_title += len(title[i])
            ndbf.write(title[i])
        else:
            len_title += len(title[i]) + 1
            ndbf.write(' ' + title[i])

ndbf.write('\n')

ndbf.write(author_string.format('AUTHOR','  ','Bintu et al., Science 362, 419 (2018)'))
ndbf.write("\n")

ndbf.write(author_string.format('ASMBLY','  ','hg38'))
ndbf.write("\n")

# [ Creating SEQCHR ]

file_csv.seek(0)

def chunk_list(lista, size):
    lista = iter(lista)
    return iter(lambda: tuple(islice(lista, size)), ())

SEQCHROM = []

beads = 0
chain = 1

for line in file_csv:
  
    if line[0].isdigit():

        info = line.split(',')
        
        if int(info[0]) ==  1:

            SEQCHROM.append('UN')
            beads += 1
        
        elif int(info[0]) != 1:
            break

seq_chunk_23 = list(chunk_list(list(SEQCHROM), 23))

for j in range(len(seq_chunk_23)):
    
    seq_str = " ".join(seq_chunk_23[j])
    ndbf.write(seqchr_string.format('SEQCHR', j+1, chroID + str(chro), beads, seq_str))
    ndbf.write("\n")


# [ Writing .ndb body ]

file_csv.seek(0)

numTer = 0

for line in file_csv:
    if line[0].isdigit():

        info = line.split(',')
        
        model = info[0]

        if int(info[1]) == 1:
            ndbf.write(model_string.format('MODEL ', int(model)))
            ndbf.write("\n")
            numTer = 0
        
        index_c = index = int(info[1])

        Z = np.float(info[2])
        X = np.float(info[3])
        Y = np.float(info[4])

        start = 18000000 + (index - 1) * res
        end   = 18000000 + (  index  ) * res

        ndbf.write(ndb_string.format('CHROM ', index_c, 'UN', " ", chroID + str(chro), index, X, Y, Z, start, end, sigma)) # Aqui a gente escreve as coordenadas e os campos coloridos
        ndbf.write("\n")

        if int(info[1]) == beads:
            ndbf.write(ter_string.format('TER   ', index_c + 1, 'UN', chroID + str(chro)))
            ndbf.write("\n")
            numTer += 1
            
            ndbf.write('ENDMDL\n')

# [ Writing loops]

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
