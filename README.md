## NDB-Converters
 Python scripts to convert from and to .ndb file format

# gro2ndb.py

Usage: gro2ndb.py -f file.gro -n name_NDB_file

Other options:

-res {int}          : Resolution for each bead
-chroID {string}    : String up to 2 characters with the chain identifier 
                    : Default = C
-chro {int}         : Chromosome number in the *.gro file to be converted
-loops {file}       : Text file that contains loop annotations
-sigma {float}      : Value for the Variance of DNA Distribution
-scale {string}     : Scale (units) of distances in the file
                    : For simulations, this value would be 1.0 sigma

# pdb2ndb.py

Usage: pdb2ndb.py -f file.pdb -n name_NDB_file

Other options:

-res {int}          : Resolution for each bead
-chroID {string}    : String up to 2 characters with the chain identifier 
                    : Default = C
-loops {file}       : Text file that contains loop annotations
-sigma {float}      : Value for the Variance of DNA Distribution
-scale {string}     : Scale (units) of distances in the file
                    : For simulations, this value would be 1.0 sigma

# csv2ndb.py

Usage: csv2ndb.py -f file.csv -n name_NDB_file

This corverter was originally wrote to convert data from Bintu et. al. (2018) Science.

Other options:

-res {int}          : Resolution for each bead
-chroID {string}    : String up to 2 characters with the chain identifier 
                    : Default = C
-loops {file}       : Text file that contains loop annotations
-sigma {float}      : Value for the Variance of DNA Distribution
-scale {string}     : Scale (units) of distances in the file
                    : For simulations, this value would be 1.0 sigma

# sw2ndb.py

Usage: sw2ndb.py -f file.sw -n name_NDB_file

Other options:

-res {int}          : Resolution for each bead
-chroID {string}    : String up to 2 characters with the chain identifier 
                    : Default = C
-loops {file}       : Text file that contains loop annotations
-sigma {float}      : Value for the Variance of DNA Distribution
-scale {string}     : Scale (units) of distances in the file
                    : For simulations, this value would be 1.0 sigma

