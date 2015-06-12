function model = VesselGenerator(BV,R,n,Inter,verb)

%  
% Generates a vessel geometry. The geometry is composed of
% cylinders randomly positioned and duplicated when they overlap the square
% sides. The geometry is longitudinally invariant. A exclusive surface
% arround each vessel might be specified where none other vessel can be
% poistionned.
%   Input parameters:
%       BV          : Blood Volume fraction
%       R           : Radius of the vessels
%       n           : Number of vessels
%       Inter       : Radius arround the vessels for exclusion
%
%   Output parameters:
%       model       : structure containing the generated geometry
%
% Adapted by N. Pannetier, 2011
% from F. Mauconduit, GIN, Grenoble, 2009 (whose code do many other things)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define Parameters
progress_str = '';
ax = sqrt(n*pi*R^2/BV);
ay = ax;
if numel(R) == 1
    rexc = ones([1 n])*(R + Inter);
    Rvec = ones([1 n])*R;
elseif numel(R) == n
    rexc = (R + Inter);
    Rvec = R;
else
    disp('Problem in the defined vessel geometry')
    return
end

p = ax * rand(n,2);
a = [ax ay];

k=1;                %% Iteration sur les n cylindre a� positionner
num=0;              %% Indice de la sphere placee
count=0;            %% Nombre total d'iterations
lim=0;              %% Nombre d'iteration limite avant de passer
pass=0;             %% Nombre de sphères passées

rep_marge = 0;      %% Marge de replication des spheres

%%% Output parameters
q=[];               %% Coordonnees des centres des n cylindre
s=[];               %% Rayons des n cylindre
ind=[];             %% Numero d'indice des cylindres (les retitions correspondent
%%% a� des spheres identiques repliquees sur les bords)

%% Loop on the n vessels
while (k <= n && num < n)
    
    count=count+1;
    lim=lim+1;
    
    if mod(count,100000)==0
        
        if lim>=5e4
            %% Try to much time without success
            warndlg('Cannot position any more vessel, watch out your input parameters, you are maybe to greedy')
            return
        end
    end
    
    %% If two vessels overlap
    if collapsed(q,s,p(k,:),rexc(k)),
        p(k,:)=[ax ay].*rand(1,2);
        continue;
    end
    
    new_p=[];
    flag=[];
    
    %% Is the vessel inside the square? %%%
    if p(k,1)+Rvec(k)>ax-rep_marge
        flag=[1 flag];
        p1=[p(k,1)-ax p(k,2)];
        if collapsed(q,s,p1,rexc(k)),
            p(k,1)=Inf; continue; end
    end
    if p(k,2)+Rvec(k)>ay-rep_marge
        flag=[2 flag];
        p2=[p(k,1) p(k,2)-ay];
        if collapsed(q,s,p2,rexc(k)),
            p(k,1)=Inf; continue; end
    end
    
    if p(k,1)-Rvec(k)<0+rep_marge
        flag=[4 flag];
        p4=[p(k,1)+ax p(k,2)];
        if collapsed(q,s,p4,rexc(k)),
            p(k,1)=Inf; continue; end
    end
    if p(k,2)-Rvec(k)<0+rep_marge
        flag=[5 flag];
        p5=[p(k,1) p(k,2)+ay];
        if collapsed(q,s,p5,rexc(k)),
            p(k,1)=Inf; continue; end
    end
    
    %% Copy vessel on the opposite side if it intersects with the border
    if ~isempty(flag)
        
        if numel(flag)==1
            eval(sprintf('new_p = p%d;',flag))
            q=[q; new_p];
            s=[s Rvec(k)];
            num=num+1;
            ind=[ind num];
            
        elseif numel(flag)==2
            p7=p(k,:);
            if sum(flag==1)
                p7(1)=p(k,1)-ax;
                if collapsed(q,s,p7,rexc(k)), p(k,1)=Inf; continue; end
            end
            if sum(flag==2)
                p7(2)=p(k,2)-ay;
                if collapsed(q,s,p7,rexc(k)), p(k,1)=Inf; continue; end
            end
            if sum(flag==3)
                p7(3)=p(k,3)-az;
                if collapsed(q,s,p7,rexc(k)), p(k,1)=Inf; continue; end
            end
            if sum(flag==4)
                p7(1)=p(k,1)+ax;
                if collapsed(q,s,p7,rexc(k)), p(k,1)=Inf; continue; end
            end
            if sum(flag==5)
                p7(2)=p(k,2)+ay;
                if collapsed(q,s,p7,rexc(k)), p(k,1)=Inf; continue; end
            end
            if sum(flag==6)
                p7(3)=p(k,3)+az;
                if collapsed(q,s,p7,rexc(k)), p(k,1)=Inf; continue; end
            end
            q1=eval(sprintf('p%d',flag(1)));
            q2=eval(sprintf('p%d',flag(2)));
            
            if collapsed(q,s,q1,rexc(k)), p(k,1)=Inf; continue; end
            if collapsed(q,s,q2,rexc(k)), p(k,1)=Inf; continue; end
            
            q=[q; p7; q1; q2];
            s=[s Rvec(k) Rvec(k) Rvec(k)];
            num=num+1;
            ind=[ind num num num];
            
        elseif numel(flag)==3
            q1=eval(sprintf('p%d',flag(1)));q12=p(k,:);
            q2=eval(sprintf('p%d',flag(2)));q13=p(k,:);
            q3=eval(sprintf('p%d',flag(3)));q23=p(k,:);q123=p(k,:);
            
            q12(p(k,:)~=q1)=q1(p(k,:)~=q1);
            q12(p(k,:)~=q2)=q2(p(k,:)~=q2);
            
            q13(p(k,:)~=q1)=q1(p(k,:)~=q1);
            q13(p(k,:)~=q3)=q3(p(k,:)~=q3);
            
            q23(p(k,:)~=q2)=q2(p(k,:)~=q2);
            q23(p(k,:)~=q3)=q3(p(k,:)~=q3);
            
            q123(p(k,:)~=q1)=q1(p(k,:)~=q1);
            q123(p(k,:)~=q2)=q2(p(k,:)~=q2);
            q123(p(k,:)~=q3)=q3(p(k,:)~=q3);
            
            if collapsed(q,s,q12,rexc(k)), p(k,1)=Inf; continue; end
            if collapsed(q,s,q13,rexc(k)), p(k,1)=Inf; continue; end
            if collapsed(q,s,q23,rexc(k)), p(k,1)=Inf; continue; end
            if collapsed(q,s,q123,rexc(k)), p(k,1)=Inf; continue; end
            
            q=[q; q1; q2; q3; q12; q13; q23; q123];
            s=[s Rvec(k) Rvec(k) Rvec(k) Rvec(k) Rvec(k) Rvec(k) Rvec(k)];
            num=num+1;
            ind=[ind num num num num num num num];
        end
        
    else
        num=num+1;
    end
    
    %% Define the center of the new placed vessel
    q=[q; p(k,:)];
    s=[s Rvec(k)];
    ind=[ind num];
    k=k+1;
    lim=0;
    
    if isempty(progress_str)
        prevLength = 0;
    else
        prevLength = numel(progress_str);
    end
    progress_str = sprintf('Number of vessels positioned: %d/%d \n',k-1,n);
    if verb, fprintf([repmat('\b',1,prevLength) '%s'],progress_str);end
    drawnow;
end

model.Name = 'Vessel_Generator';
model.Centers = q;
model.Radius = s;
model.Indices = ind;
model.CubeSize = a;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% SUBFUNCTIONS %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = collapsed(Cn,rn,C1,r1)
%%  Verifie si le cylindre C1 se superpose sur l'un des cylindres Cn

if isinf(C1(1)), res = 1; return; end
if isempty(Cn),  res = 0; return; end

X = C1(:,1) - Cn(:,1);
Y = C1(:,2) - Cn(:,2);
%Z = C1(:,3) - Cn(:,3);

dist = sqrt( X.^2 + Y.^2 );

sumR = r1 + rn';

if min(dist - sumR)<0
    res=1;
elseif min(dist - sumR)>=0
    res=0;
else
    res=[];
end

function out = FracVol(a,q,s,ind)
%% Calcul la fraction volumique intra-cellulaire de la distribution

[dum indice]=unique(ind);
r = s(indice);
out = sum(pi*r.^2)/prod(a);

function Fractions(a,q,s)
% Calcul la fraction volumique intra-cellulaire de la distribution
% Méthode MonteCarlo (calcul long)
inCount=0;
N=0;
fp=Inf;
fn=Inf;
precision=0;

while N<5 || precision<=30
    
if abs(1-fp/fn)<1e-5
    precision=precision+1;
else
    precision=0;
end

P = a*rand(1,2);
[in box]=IsInCylinder(P,q,s);
if in==1, inCount=inCount+1; end
N=N+1;
fp=fn;
fn=inCount/N;
%disp(sprintf('%.4g',abs(1-fp/fn)))
end
% disp(sprintf('Volume intra cellaulaire: %.3g\nNombre d''iterations: %d',fn,N))

function [res box]=IsInCylinder(P,q,s)

X = P(:,1) - q(:,1);
Y = P(:,2) - q(:,2);
%Z = P(:,3) - q(:,3);

dist = sqrt( X.^2 + Y.^2 ) -s';
[value box]=min(dist);

if value<0
    res=1;
elseif value>=0
    res=0;
else
    res=[];
end
