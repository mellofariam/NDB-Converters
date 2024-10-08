#!/usr/bin/env python
# coding: utf8

__description__ = """
Converting *.cndb to Nucleome Data Bank format .ndb
"""

__author__ = "Vinícius Contessoto / Matheus Mello / Antonio B. Oliveira Junior"
__date__ = "Set/2020"

################################################################
#
# Trajectories file *.cndb to Nucleome Data Bank format .ndb
#
# usage:
#  ./cndb2ndb.py -f file.cndb -n name_NDB_file
#
################################################################


import time
import argparse
import numpy as np
import h5py
from itertools import islice

parser = argparse.ArgumentParser(description="Converting from *.cndb to *.ndb file")
parser.add_argument(
    "-f",
    metavar="input-file-CNDB-frames",
    help="cndb file",
    action="store",
    dest="arg_fname",
    required=True,
)
parser.add_argument(
    "-res",
    action="store",
    default=50000,
    dest="arg_res",
    help="Resolution for each simulation bead",
)
parser.add_argument(
    "-chroID", action="store", default="C", dest="arg_chroID", help="Chain ID"
)
parser.add_argument(
    "-n",
    action="store",
    default="chromatin",
    dest="arg_name",
    help="Name of output file",
)
parser.add_argument(
    "-loops",
    metavar="input-file-Loops",
    help="Loops contact pair between i and j",
    type=argparse.FileType("rt"),
    required=False,
)
parser.add_argument(
    "-sigma",
    action="store",
    default=0.000,
    dest="arg_sigma",
    help="Distance fluctuation",
)
parser.add_argument(
    "-scale", action="store", default=1.000, dest="arg_scale", help="Distance scale"
)

try:
    arguments = parser.parse_args()
    print("################################################")
    print("Chosen file: {:}".format(arguments.arg_fname))
    print("Resolution: {:}".format(arguments.arg_res))
    print("Chrom ID: {:}".format(arguments.arg_chroID))
    print("Distance scale: {:}".format(arguments.arg_scale))

except IOError as msg:
    parser.error(str(msg))

Main_chrom = ["ChrA", "ChrB", "ChrU"]  # Type A B and Unknow
Chrom_types = ["ZA", "OA", "FB", "SB", "TB", "LB", "UN"]
Chrom_types_NDB = ["A1", "A2", "B1", "B2", "B3", "B4", "UN"]

##################################################################################################

b_time = time.time()

gro_string = "{0:5d}{1:5s}{2:5s}{3:5d}{4:8.3f}{5:8.3f}{6:8.3f}"
gro_box_string = "{0:10.5f}{1:10.5f}{2:10.5f}"
ndb_string = "{0:6s} {1:8d} {2:2s} {3:6s} {4:4s} {5:8d} {6:8.3f} {7:8.3f} {8:8.3f} {9:10d} {10:10d} {11:8.3f}"
pdb_string = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}"
header_string = "{0:6s}    {1:40s}{2:9s}   {3:4s}\n"
title_string = "{0:6s}  {1:2s}{2:80s}\n"
author_string = "{0:6s}  {1:2s}{2:79s}\n"
expdata_string = "{0:6s}  {1:2s}{2:79s}\n"
model_string = "{0:6s}     {1:4d}"
seqchr_string = "{0:6s} {1:3d} {2:2s} {3:5d}  {4:69s}\n"
ter_string = "{0:6s} {1:8d} {2:2s}        {3:2s}"
loops_string = "{0:6s} {1:8d} {2:8d}"
master_string = "{0:6s} {1:8d} {2:6d} {3:6d} {4:10d}"  # MASTER BEADS TER LOOPS RES

filename = arguments.arg_fname
file_loops = arguments.loops
res = int(arguments.arg_res)
chroID = str(arguments.arg_chroID)
sigma = float(arguments.arg_sigma)
scale = float(arguments.arg_scale)

ndbf = open(arguments.arg_name + ".ndb", "w+")

mode = "r"
cndbf = h5py.File(filename, mode)

type_list = list(cndbf["types"])

if isinstance(type_list[0], (int, np.int32, np.int64)):
    type_list = [Chrom_types_NDB[i] for i in type_list]


### [ HEADER ]
ndbf.write(
    header_string.format("HEADER", "NDB File generated by Open-MiChroM", " ", " ")
)
ndbf.write(
    title_string.format("TITLE ", "  ", "A Scalable Computational Approach for ")
)
ndbf.write(
    title_string.format("TITLE ", "2 ", "Simulating Complexes of Multiple Chromosomes")
)
ndbf.write(expdata_string.format("EXPDTA", "  ", "Cell Line  @50k bp resolution"))
ndbf.write(expdata_string.format("EXPDTA", "  ", "Simulation - Open-MiChroM"))
ndbf.write(author_string.format("AUTHOR", "  ", "Antonio B. Oliveira Junior - 2020"))


# [ Creating SEQCHR ]
def chunks(l, n):
    n = max(1, n)
    return [l[i : i + n] for i in range(0, len(l), n)]


seq_chunks= chunks(type_list, 23)

for num, line in enumerate(seq_chunks):
    ndbf.write(
        seqchr_string.format("SEQCHR", num + 1, "C1", len(type_list), " ".join(line))
    )

# [ Writing .ndb body ]
model = 1
chain = 1
index_c = 0
try:
    loop = np.array(cndbf["loops"])
    totalModel = len(cndbf) - 1
    write_loops = True
except KeyError:
    totalModel = len(cndbf)
    write_loops = False

for mdl in range(1, totalModel, 1):
    ndbf.write(model_string.format("MODEL ", int(mdl)))
    ndbf.write("\n")

    data = np.array(cndbf[str(mdl)])

    for i, line in zip(list(range(len(data))), data):
        ndbf.write(
            ndb_string.format(
                "CHROM ",
                i + 1,
                type_list[i],
                " ",
                "C1",
                i + 1,
                line[0],
                line[1],
                line[2],
                int((i) * 50000) + 1,
                int(i * 50000 + 50000),
                0,
            )
        )  
        ndbf.write("\n")

    ndbf.write(ter_string.format("TER", len(data) + 1, type_list[-1], "C1"))
    ndbf.write("\nENDMDL\n")
    model += 1
print("Last traj frame:", model)

if write_loops:
    loop = np.sort(loop, axis=0)
    for p in loop:
        ndbf.write(loops_string.format("LOOPS", p[0], p[1]))
        ndbf.write("\n")
ndbf.write("END")
ndbf.write("\n")

ndbf.close()

e_time = time.time()
elapsed = e_time - b_time
print("Ran in %.3f sec" % elapsed)
