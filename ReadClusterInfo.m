function [Cluster] = ReadClusterInfo(File)

% Read and parse the parameter file
% Input: 
%   - File  : Path to the parameter file
% Output:

%
% Code by N.Pannetier, Feb 14th 2013

Cluster = struct;

%% read lines
fid = fopen(File,'rt');
C = textscan(fid, '%s', 'Delimiter',''); C = C{1};
fclose(fid);

%% Search for Model and Seq structure in the file
for a=1:numel(C)
    ClusterLine(a) = strncmp(C{a},'Cluster.',8);
end
CLine = find(ClusterLine);

%% Read Cluster strucutre
for a=1:numel(CLine)
   eval(C{CLine(a)}); 
end