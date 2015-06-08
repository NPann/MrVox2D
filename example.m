% The following peace of code simulate the MR signal from a voxel as defined
% in 'config/voxpar_single.txt' when running a GESFIDE MR pulse sequence as
% defined in 'config/seqpar_GESFIDE.txt'. A very simplistic display will pop
% up that shows the MR signal and the some voxel configurations (and SLOW
% DOWN the simulation)

addpath(genpath(fullfile(pwd, 'mrvox')));

%% Test a simple example
fprintf('Simulating a GESFIDE sequence on a simplistic voxel...\n');
[sa, sphi] = VoxelSim2D_do_one('config/voxpar_single.txt', 'config/seqpar_GESFIDE.txt');

%% Test a dictionary generation 
fprintf('Simulating a dictionary on a simplistic voxel with varying vessel radius...\n');
dico = GenLookUp2D('config/voxpar_dico.txt', 'config/seqpar_GESFIDE.txt');

% Just to get input parameters for plotting legend
Model   = ReadModel('config/voxpar_dico.txt');
figure,plot(abs(dico')), xlabel('Time (ms)'), ylabel('MR signal magnitude')
title('Dictionary')
legend(cellstr(num2str(Model.geo.vasc.R'*1e6, 'R=%-dum')))
