# README #

MrVox2D is a toolkit written in Matlab to simulate the magnetic resonance 
signal considering the effect of susceptibility inclusions and wather diffusion.

### What is MrVox2D for? ###

* MrVox2D is for simulating the MR complex signal given:
    1. a voxel with a user-defined microstructure
    2. a MR pulse sequence  
* MrVox is designed with versatility in mind and user can control:
    1. the geometry of the voxel (number of blood vessels and cells, size, spacing, etc.)
    2. the MR-related properties of the different compartments (T1, T2, M0, susceptibility)
    3. the magnetic field and its variation over the voxel dimension (B0 orientation, linear gradient)
    4. the water diffusion (diffusivity, hindered or not)
* MrVox2D is particularly designed to generate dictionaries of MR signals with varying input properties (e.g. blood volume, vessel size, oxygenation). 
For example, it has been used for the vascular fingerprinting framework 
[(Christen et al. NeuroImage 20014)](http://www.sciencedirect.com/science/article/pii/S1053811913012019)
* **Specifically**: The voxel is considered as a 2D plane and the magnetic inclusions (e.g. vessels) are disks randomly spread in this 2D plane. Voxel and sequence parameters are defined in a text file. The MR sequence is passed as a function handle of the simulator. Perturbations of the magnetic field by the susceptibility inclusions are considered. Diffusion is modeled by convolution of a gaussian kernel and can be hindered by cell and vessel compartment walls. Even though the simulation is 2D, the MR signal is similar to the one obtained in a 3D voxel with isotropic vessel orientation when the number of vessels is "high enough" [(Pannetier et al. Plos. 2014)](http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0057636).
* **Limitations**: 
    * The simulation of the diffusion is modeled by a convolution kernel and the computation is performed in the Fourier domain. The related aliasing effect have some consequences: 
        1. The geometry lattice must be periodic
        2. When applying gradient, their intensity must be such that the dephasing over the voxel extent during the simulation step time dt is modulo 2 x pi.
    * IMPORTANT: It has been reported that the approach for generating the field offset map (average of field maps generated with B0 in 3 orthogonal directions) was not the most accurate when comparing to a standard 3D approach. However, excellent agreement with the 3D was reported when using a random B0 direction for each cylinder (as [here](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.21690.)), which is not implemented in this repo.
 
### How do I get set up? ###

* Add mrvox/ to your matlab path
* Check out the example in example.m
* Start tinkering with your own configuration files
* Examples of configuration files for the voxel and MR sequences can be found in config/.
* The two main use cases are:

    * **Single/simple use**: to run the simulator on a voxel geometry and MR sequence:
```
[Sa, Sphi] = VoxelSim2D_do_one('config/voxpar_single.txt','config/seqpar_GESFIDE.txt')
```

* or:
    * **Dictionary use**: to generate a dictionary of MR signals by varying *any* parameters defined as an array in the Model structure of the configuration file (Model.phy, Model.vox or Model.geo):
```
Dico = GenLookUp2D('config/voxpar_dico.txt','config/seqpar_GESFIDE.txt')
```

* **Deployment**: Dictionary generation is highly parrallelizable and the code
provides a support for using Matlab Distributed Computing Server (MDCS). Server
configuration must be defined in the configuration file (see config/cluster_info.txt for an example).
The example file is compatible with Matlab R2011b and the code will likely need updates
to be used with newer version. The code can also be compiled using the Matlab compiler and deployed with another scheduler


### Contribution guidelines ###

Through the issues tracker 

### Who do I talk to? ###

Nicolas Pannetier  
Thomas Christen  
Clement Debacker  

### License ###

MrVox2D is licensed under the terms of the BSD license. Please see the License file in the repository

### Citation ###

If you use this simulator and need to reference it, please cite:  

* MR Vascular Fingerprinting: A New Approach to Compute Cerebral Blood Volume, Mean Vessel Radius, and Oxygenation Maps in the Human Brain T. Christen, NA. Pannetier, W. Ni, D. Qiu, M. Moseley, N. Schuff, and G. Zaharchuk. Neuroimage. 2014 Apr 1; 89: 262–270.

* or: A Simulation Tool for Dynamic Contrast Enhanced MRI. NA. Pannetier, CS. Debacker, F. Mauconduit, T. Christen, EL. Barbier. Plos One 2014. 10.1371/journal.pone.0057636

![mrvox2d_illustration.png](https://bitbucket.org/repo/G9qXRy/images/877007194-mrvox2d_illustration.png)
