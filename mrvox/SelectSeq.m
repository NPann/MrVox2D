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