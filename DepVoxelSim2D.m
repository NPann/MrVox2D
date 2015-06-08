function [Sa, Sphi] =  DepVoxelSim2D(FileIn,Cluster,ModelLine)

% Execute the compiled version of VoxelSim2D

Path2Script = which('run_VoxelSim2D.sh');
command = sprintf('system(''%s %s %s %d '')',Path2Script,Cluster.MRL.path,FileIn,ModelLine);
s = evalc(command);
expressn1 = 'Sphi =(.*?)END';
expressn2 = 'Sa =(.*?)END';
[matchstr ~] = regexp(s,expressn1,'tokens','match');
Sphi = str2num(matchstr{1}{1});
[matchstr ~] = regexp(s,expressn2,'tokens','match');
Sa = str2num(matchstr{1}{1});