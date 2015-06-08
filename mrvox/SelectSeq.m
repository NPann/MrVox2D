% Copyright 2012, 2013 Nicolas Pannetier, Cl�ment Debacker, Franck Mauconduit, Thomas Christen, Emmanuel Barbier
% This file is part of DCESIM.
%
%   DCESIM is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%    DCESIM is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with DCESIM.  If not, see <http://www.gnu.org/licenses/>.

function [Seq] = SelectSeq(File,Line)

% Read and parse the parameter file
% Input:
%   - File  : Path to the parameter file
% Output:
%   - Model     : Structure containing the voxel related parameters
%   - Seq       : Structure containing the sequence related parameters
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

Seq = struct;

%% read lines
fid = fopen(File,'rt');
if fid > 0
    C = textscan(fid, '%s', 'Delimiter',''); C = C{1};
else
    fprintf('Cannot open Parameter file\n');
end
fclose(fid);

%% Search for Seq structure in the file or variable tmp
for a=1:numel(C)
    SeqLine(a) = (strncmp(C{a},'Seq.',4) || strncmp(C{a},'tmp',3));
end
SLine = find(SeqLine);

%% Read Seq structure
for a=1:numel(SLine)
    eval(C{SLine(a)});
end

%% Select Sequence
cname = 'Seq';
DicInp = [];
if isstruct(eval(cname))
    DicInp = GoDownStruc(cname,DicInp,Seq);
else
    fprintf('Input is not a Seq structure\n');
end
if ~isempty(DicInp)
    eval(sprintf('%s = %s(%d,:);',DicInp{1}.name,DicInp{1}.name,Line));
else
    Seq =Seq;
end

end


% Recursive function
function DicInp = GoDownStruc(cname_o,DicInp,Seq)
a=1;
if isstruct(eval(cname_o))
    Fields = fieldnames(eval(cname_o));
    while a <= numel(Fields)
        cname = sprintf('%s.%s',cname_o,Fields{a});
        DicInp = GoDownStruc(cname,DicInp,Seq);
        a = a+1;
    end
else
    if size(eval(cname_o),1) > 1 && ~ischar(eval(cname_o))
        DicInp{numel(DicInp)+1}.name = cname_o;
    end
    return
end
end