function [Label, par] = GenerateSeqTable(Seq)

% Explore the structure to find arrays in it

cname = 'Seq';
DicInp = [];
if isstruct(eval(cname))
    DicInp = GoDownStruc(cname,DicInp,Seq);
else
    fprintf('Input is not a Seq structure\n');
end

if ~isempty(DicInp)
    Label{1} = DicInp{1}.name;
    strtmp = DicInp{1}.name;
    strtmp2 = 'A';
    strtmp3 = 'A(:)';
    abc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    for a=2:numel(DicInp)
        strtmp = sprintf('%s, %s',strtmp,DicInp{a}.name);
        strtmp2 = sprintf('%s, %s',strtmp2,abc(a));
        strtmp3 = sprintf('%s, %s(:)',strtmp3,abc(a));
        Label{a} = DicInp{a}.name;
    end
    eval(sprintf('[ %s ] = ndgrid(%s);',strtmp2, strtmp));
    eval(sprintf('par = [%s];',strtmp3));
else
    Label = [];
    par = [];
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
    if size(eval(cname_o),2) > 1 && ~ischar(eval(cname_o))
        DicInp{numel(DicInp)+1}.name = cname_o;
    end
    return
end
end
