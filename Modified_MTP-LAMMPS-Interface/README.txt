To compile the modified MTP/LAMMPS interface, the source code file within this directory
are used as a direct replacement for source code file from the original package. 

The file "pair_MLIP.cpp" and "pair_MLIP.h" should be replacing the file within the "LAMMPS-MLIP interface-v2" package at https://gitlab.com/ashapeev/interface-lammps-mlip-2.
Overwrite and replace the files from the package at the follow location:
replace ./pair_MLIP.cpp -> interface-lammps-mlip-2/LAMMPS/USER-MLIP/pair_MLIP.cpp
replace ./pair_MLIP.h -> interface-lammps-mlip-2/LAMMPS/USER-MLIP/pair_MLIP.h

The file "interface.cpp" should be replacing the file within the "MLIP-v2" package at https://gitlab.com/ashapeev/mlip-2.
Overwrite and replace the file from the package at the follow location:
replace ./interface.cpp -> mlip-2/src/external/interface.cpp

After the replacement, re-compile the mlip interface libirary with make libinterface.
Then, the MLIP compatible LAMMPS execution can be re-compile with the new lib_mlip_interface.a following the installation procdure of the "LAMMPS-MLIP interface-v2" package.