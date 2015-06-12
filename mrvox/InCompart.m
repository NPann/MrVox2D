function compart=InCompart(M,q,s)
% function compart=InCompart(M,q,s)
%  
% Renvoie le numéro de compartiment d'un point ou d'un ensemble de points M
% dans une géométrie de matière blanche définie par un ensemble de cylindres.
%   Input parameters:
%       M    : Coordinate vector of points (Npts x 2)
%       q    : Coordinate vector of the center of the N cylinders generated (Nx2)
%       s    : Radius vector of the N spheres (1xN)
%
%   Output parameters:
%    compart : Vector indicating the compartiment the Npts points
%
% Code adapted F. Mauconduit, GIN, Grenoble, 2010

r1=[];r2=[];r3=[];
npart=size(M,1);

dist = ...
    (( M(:,1)*ones(1,length(s)) - ones(npart,1)*q(:,1)' ).^2 +...
     ( M(:,2)*ones(1,length(s)) - ones(npart,1)*q(:,2)' ).^2 ).^(1/2);

pos = dist<ones(npart,1)*s;

[r1 c1] = find(pos(sum(pos,2)==1,:)'); % c1 est necessaire pour la sortie sur r1 (voir "help find")

if any( find(pos(sum(pos,2)==2,:)) )
    [row col] = find(pos(sum(pos,2)==2,:)');
    signe = VectorialCalc(M, sum(pos,2)==2, q(row,:), s(row)');
    r2 = row((1:2:length(row))+signe');
end

if any( find(pos(sum(pos,2)==3,:)) )
    [row col] = find(pos(sum(pos,2)==3,:)');
    row2 = row; row2(3:3:end)=[];
    signe1 = VectorialCalc(M, sum(pos,2)==3, q(row2,:), s(row2)');
    row([1:3:length(row)]+~signe1')=[];
    signe2 = VectorialCalc(M, sum(pos,2)==3, q(row,:), s(row)');
    r3 = row((1:2:length(row))+signe2');
end

if any( find(pos(sum(pos,2)>3,:)) )
    warning('A particle has found to be in more than 3 compartments')
    compart=NaN;
    return
end

compart=zeros(npart,1);
compart(sum(pos,2)==1) = r1;
compart(sum(pos,2)==2) = r2;
compart(sum(pos,2)==3) = r3;



function signe = VectorialCalc(M, sum_pos, O, Rayons)
% calcul vectoriel pour trouver le cote où se trouve le point M lorsqu'il
% est à l'interieur de deux sphères.

u12 = O(2:2:end,:) - O(1:2:end,:);
O1_O2 = sqrt(sum(u12.^2,2));
I_O1 = ( O1_O2.^2 + Rayons(1:2:end).^2 - Rayons(2:2:end).^2 ) ./ (2*O1_O2);
I = O(1:2:end,:) + I_O1*ones(1,size(u12,2)) .* ( u12 ./ (O1_O2*ones(1,size(u12,2))) );

signe = sign( sum( ( M(sum_pos,:) - I ) .* u12 , 2) )/2 + 1/2;
