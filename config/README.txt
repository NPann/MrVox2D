This package model the MR signal of an MRI voxel taking into account the diffusion
of the protons within a local pertubated magnetic field produced by magnetic inclusions.
The voxel is considered as 2D plane and the magnetic inclusions are disks randomly spread
in the 2D plane.
Different MRI sequence can be acquired and 2 examples is provided 
(n RF excitations sequence and diffusion single refocused narrow pulse sequence)

All the input parameters must be define in a .txt file.

## SINGLE USE
The main function is:

[Sa, Sphi] = VoxelSim2D_do_one(FileName,LineNum)
Input:
    - FileName: Path to the file containing the structure of the Model and
              the Seququence
    - LineNum:  Line Index of the set of input parameters that has to be simulated (optional, default = 1)
Output:
    - Sa:   Magnitude of the MR signal
    - Sphi: Phase of the MR signal

Ex: [Sa, Sphi] = VoxelSim2D_do_one('Param1.txt',1)

2 files example are provided (Param1.txt and Param2.txt)
Input parameters are described within these files.


## DICTIONARY
This package can also be used to build up a dictionary of MR signals
Any parameters X in the structure Model.phy Model.vox or Model.geo with size(X,2) > 1 will
be considered for dictionary building.

The Main function is:
Dico = GenLookUp2D(FileIn,PathOut,ServerInfo)
Input: - FileIn: path to the file where the Model and the MR Sequence are defined
               (Must be a a cell of strings for multiple files)
       - PathOut: path to the repertory where the dictionaries will be saved

 Optional: - ServerInfo: Path to the file defining the cluster info for
 parallel jobs (required MCDS + Matlab Runtime Library installed)

Ex: GenLookUp2D('LT.txt','~','ClusterInfo.txt')


## DEPLOYEMENT

The Code can be deployed using Matlab parallel toolbox but can also be compiled using
the Matlab compiler and deployed with another scheduler (or the Matlab Scheduler).
When deployed, the clusterinfo.txt file must defined accordingly.
The compiled code provide has been compiled with MatlabR2011b under Linux64b.


## Log
v1.0 2013, Feb 13rd   
    - Release
v1.1 2013, April 4th  
    - Bug fixed
    - Adds the B0 orientation as an option in the parameter file
    - Changes management of input dictionary
    - Allows to generate the dicitonary with what ever input X with size(X,2) > 1
    - Adds the option to cast central part of the voxel only (avoid side effect on applying gradient lower than 2pi over the voxel)
    - Support basic "imaging" kind of sequence
    - Support hindered diffusion in the cell compartments
v1.2 2013, June 6th
    - Sequence and Voxel Model are now separated into files
    - RF pulse structure fixed and updated