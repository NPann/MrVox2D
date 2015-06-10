% The following peace of code simulate the MR signal from a voxel as defined
% in 'config/voxpar_single.txt' when running a GESFIDE MR pulse sequence as
% defined in 'config/seqpar_GESFIDE.txt'. A very simplistic display will pop
% up that shows the MR signal and the some voxel configurations (and SLOW
% DOWN the simulation)

addpath(genpath(fullfile(pwd, 'mrvox')));

%% FID example
fprintf('Simulating a FID on simple voxel...\n');
[sa, sphi] = VoxelSim2D_do_one('config/voxpar_simple.txt', 'config/seqpar_fid.txt');

%% Spin echo example
fprintf('Simulating a spin echo on simple voxel...\n');
[sa, sphi] = VoxelSim2D_do_one('config/voxpar_simple.txt', 'config/seqpar_se.txt');

%% Spin echo example with strong susceptibility and reduced diffusion
fprintf('Simulating a spin echo on simple voxel...\n');
[sa, sphi] = VoxelSim2D_do_one('config/voxpar_simple_highsusc_lowdiff.txt', 'config/seqpar_se.txt');

%% Stimulated echo example
fprintf('Simulating a spin echo on simple voxel...\n');
[sa, sphi] = VoxelSim2D_do_one('config/voxpar_simple_highsusc_lowdiff.txt', 'config/seqpar_stim.txt');

%% Narrow diffusion gradient example with hindered diffusion and cells
fprintf('Simulating a narrow diffusion gradient on simple voxel...\n');
[sa, sphi] = VoxelSim2D_do_one('config/voxpar_simple_hindered_cell.txt', 'config/seqpar_narrowpulse_singlerefoc.txt');

%% Dictionary generation example
fprintf('Simulating a dictionary on a simplistic voxel with varying vessel radius...\n');
dico = GenLookUp2D('config/voxpar_dico.txt', 'config/seqpar_GESFIDE.txt');

% Just to get input parameters for plotting legend
Model   = ReadModel('config/voxpar_dico.txt');
figure,plot(abs(dico')), xlabel('Time (ms)'), ylabel('MR signal magnitude')
title('Dictionary')
legend(cellstr(num2str(Model.geo.vasc.R'*1e6, 'R=%-dum')))
