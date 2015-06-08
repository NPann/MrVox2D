# README #

MrVox is a toolkit writen in Matlab to simulate the Magnetic Resonance 
signal


### What is MrVox for? ###

* MrVox is for simulating the MR complex signal given:
    1. a voxel with a user-defined microstructure
    2. a MR pulse sequence  
* MrVox is designed with versatility in mind and user can control:
    1. the geometry of the voxel (size, number of blood vessels and cells, their size, their spacing, etc.)
    2. the MR-related properties of the different compartments (T1, T2, M0, susceptibility)
    3. the magnetic field and its variation over the voxel dimension
    4. the water diffusion  
* MrVox is particularly designed to generate dictionaries of MR signals with varying input properties (blood volume, vessel size). 
For example, it has been used succesfully in the vascluar fingerprinting framework 
[(Christen et al. NeuroImage 20014)](http://www.sciencedirect.com/science/article/pii/S1053811913012019)
* **Specifically**: The voxel is considered as a 2D plane and the magnetic inclusions (e.g. vessels) 
are disks randomly spread in the 2D plane. Voxel and sequence parameters are defined in a text file.
MR sequence is passed as a function of the simulator.

### How do I get set up? ###

* Example of configurations files for the voxel and the MR sequences are found in config/
* Single/simple use:
```
#!matlab
[Sa, Sphi] = VoxelSim2D_do_one(FileName, LineNum)
```

Input:
    * FileName: Path to the file containing the structure of the Model and the Seququence
    * LineNum:  Line Index of the set of input parameters that has to be simulated (optional, default = 1)
Output:
    * Sa:   Magnitude of the MR signal
    * Sphi: Phase of the MR signal

Ex: 
```
#!matlab
[Sa, Sphi] = VoxelSim2D_do_one('config/Param1.txt',1)
```

## Dictionary use ##
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



### Contribution guidelines ###

* 

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact