# NDB-Converters
 Python scripts to convert from and to .ndb file format

## To *.ndb

### gro2ndb.py

Usage: gro2ndb.py -f file.gro -n name_NDB_file

Other options:

-res {int}          : Resolution for each bead
-chroID {string}    : String up to 2 characters with the chain identifier 
                    : Default = C
-chro {int}         : Chromosome number in the *.gro file to be converted
-loops {file}       : Text file that contains loop annotations
-sigma {float}      : Value for the Variance of DNA Distribution

### pdb2ndb.py

Usage: pdb2ndb.py -f file.pdb -n name_NDB_file

Other options:

-res {int}          : Resolution for each bead
-chroID {string}    : String up to 2 characters with the chain identifier 
                    : Default = C
-loops {file}       : Text file that contains loop annotations
-sigma {float}      : Value for the Variance of DNA Distribution


### csv2ndb.py

Usage: csv2ndb.py -f file.csv -n name_NDB_file

This corverter was originally wrote to convert data from Bintu et. al. (2018) Science.

Other options:

-res {int}          : Resolution for each bead
-chroID {string}    : String up to 2 characters with the chain identifier 
                    : Default = C
-loops {file}       : Text file that contains loop annotations
-sigma {float}      : Value for the Variance of DNA Distribution


### spw2ndb.py

Usage: spw2ndb.py -f file.spw -n name_NDB_file

Other options:

-res {int}          : Resolution for each bead
-chroID {string}    : String up to 2 characters with the chain identifier 
                    : Default = C
-loops {file}       : Text file that contains loop annotations
-sigma {float}      : Value for the Variance of DNA Distribution

## From *.ndb:

### ndb2pdb.py

Usage: ndb2pdb.py -f file.ndb -n name_PDB_file

### ndb2spw.py

Usage: ndb2spw.py -f file.ndb -n name_SPW_file

Other options:

-t {string}         : Value of field "name" at spw header
-g {string}         : Assembly of genome
-c {string}         : Number of chromosome
