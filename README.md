# Radiation
A CUDA/C++ program to calculate the radiation emission from a high energy electron in a laser field using the Semi-classical approach of Baier et al. 
This code was used to produce the results shown in Phys. Rev. D 100, 116001 and any use of this code must be accompanied with a reference to this paper. 
The output is the file spectrum.txt where the first column is the energy of the emitted photon and the proceeding 8 columns 
the emission rates of the different combinations of initial and final spin and polarization between laser and emitted photon. 
Parameters which you may wish to change are in constants.h
