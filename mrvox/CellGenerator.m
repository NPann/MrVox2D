function model = CellGenerator(R0,Inter,fextra,radii,center,ax,vasc_ind,verb)

% Generates a White Matter model geometry. The geometry is composed of
% cylinders randomly positioned and duplicated when they overlap the square
% sides. The geometry is longitudinally invariant.
%   Input parameters:
%       a    : Square size
%       n    : Number of sphere to be generated
%
%   Output parameters:
%       q    : Coordinate vector of the center of the N spheres generated (Nx3)
%       s    : Radius vector of the N spheres (1xN)
%       ind  : Indice vector of the spheres. Same indices mean a same sphere
%              duplicated on the other side of the cibe.
%
% Code adapted from F. Mauconduit, GIN, Grenoble, 2009

% Distribution des cylindres
progress_str = '';
NbTir = round((1 - fextra) * ax^2/(pi * R0^2) * 4);
Tirage = 2*rand([1 NbTir]);

% Distribution en considerant tirage sur theta uniforme
DistriRayCell = 2*R0 * exp(-Tirage);

DiametreAxonal = [2*radii sort(DistriRayCell,'descend')];
DiametreAxonal(DiametreAxonal < 1)=[];
n = numel(DiametreAxonal);          %% Nombre de sphères
rmyelin = DiametreAxonal/2;         %% Rayon gaine myélinique

ay = ax;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calcul du reste des parametres %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = repmat(rmyelin,[1 n]);
Radius = 'constant';
Myelin = 'No';
r_Intra = r*0.6;
r_Myelin = r;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Definitions de parametres locaux %%%
%%% p: Definitions initale des centres des spheres
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = [ center ;ax * rand(n,2)];

% p = ax * rand(n,2);
Centers = 'defined';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% FIN %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


a = [ax ay];

k=1;                %% Itération sur les n sphères à positionner
num=0;              %% Indice de la sphère placée
count=0;            %% Nombre total d'itérations
lim=0;              %% Nombre d'itération limite avant de passer
pass=0;             %% Nombre de sphères passées

rep_marge = 0;   %% Marge de réplication des sphères

%%% Paramètres en sortie %%%
q=center;               %% Coordonnées des centres des n sphères
s=radii;               %% Rayons des n sphères
ind=vasc_ind;             %% Numero d'indice des sphères (les répétitions correspondent 
                    %%% à des sphères identiques répliquées sur les bords)

if fextra <= (1 - FracVol(a,q,s,ind))
    while (k <= n && fextra <= (1-FracVol(a,q,s,ind)))
        %% Loop to positione cells
        count=count+1;
        lim=lim+1;


        if mod(count,10000)==0
            if lim>=5e4
                %% Iteration limite depassee: abandon de la sphere k

                %splot(q,s)
                %warndlg('Impossible de reach the number of cells in the geometry')
                k=k+1;
                pass=pass+1;
                lim=0;
                continue
            end
        end

        % If the cell collapsed with cells or vessels
        if collapsed(q,s,p(k,:),r(k) + Inter),
            p(k,:)=[ax ay].*rand(1,2);
            continue;
        end

        new_p=[];
        flag=[];

        %% La sphère k est elle à l'intérieur du cube ? %%%
        if p(k,1)+r(k)>ax-rep_marge
            flag=[1 flag];
            p1=[p(k,1)-ax p(k,2)];
            if collapsed(q,s,p1,r(k)+Inter),
                p(k,1)=Inf; continue; end
        end
        if p(k,2)+r(k)>ay-rep_marge
            flag=[2 flag];
            p2=[p(k,1) p(k,2)-ay];
            if collapsed(q,s,p2,r(k)+Inter),
                p(k,1)=Inf; continue; end
        end
        %     if p(k,3)+r(k)>az-rep_marge
        %         flag=[3 flag];
        %         p3=[p(k,1) p(k,2) p(k,3)-az];
        %         if collapsed(q,s,p3,r(k)),
        %             p(k,1)=Inf; continue; end
        %     end
        if p(k,1)-r(k)<0+rep_marge
            flag=[4 flag];
            p4=[p(k,1)+ax p(k,2)];
            if collapsed(q,s,p4,r(k)+Inter),
                p(k,1)=Inf; continue; end
        end
        if p(k,2)-r(k)<0+rep_marge
            flag=[5 flag];
            p5=[p(k,1) p(k,2)+ay];
            if collapsed(q,s,p5,r(k)+Inter),
                p(k,1)=Inf; continue; end
        end
        %     if p(k,3)-r(k)<0+rep_marge
        %         flag=[6 flag];
        %         p6=[p(k,1) p(k,2) p(k,3)+az];
        %         if collapsed(q,s,p6,r(k)),
        %             p(k,1)=Inf; continue; end
        %     end

        %% Duplication de la sphère k si elle dépasse les cotés du cube
        % Différent cas selon qu'elle dépasse selon une seule direction (x, y ou z)
        % ou selon plusieurs directions.

        if ~isempty(flag)

            if numel(flag)==1
                eval(sprintf('new_p = p%d;',flag))
                q=[q; new_p];
                s=[s r(k)];
                num=num+1;
                ind=[ind num];

            elseif numel(flag)==2
                p7=p(k,:);
                if sum(flag==1)
                    p7(1)=p(k,1)-ax;
                    if collapsed(q,s,p7,r(k)+Inter), p(k,1)=Inf; continue; end
                end
                if sum(flag==2)
                    p7(2)=p(k,2)-ay;
                    if collapsed(q,s,p7,r(k)+Inter), p(k,1)=Inf; continue; end
                end
                if sum(flag==3)
                    p7(3)=p(k,3)-az;
                    if collapsed(q,s,p7,r(k)+Inter), p(k,1)=Inf; continue; end
                end
                if sum(flag==4)
                    p7(1)=p(k,1)+ax;
                    if collapsed(q,s,p7,r(k)+Inter), p(k,1)=Inf; continue; end
                end
                if sum(flag==5)
                    p7(2)=p(k,2)+ay;
                    if collapsed(q,s,p7,r(k)+Inter), p(k,1)=Inf; continue; end
                end
                if sum(flag==6)
                    p7(3)=p(k,3)+az;
                    if collapsed(q,s,p7,r(k)+Inter), p(k,1)=Inf; continue; end
                end
                q1=eval(sprintf('p%d',flag(1)));
                q2=eval(sprintf('p%d',flag(2)));

                if collapsed(q,s,q1,r(k)+Inter), p(k,1)=Inf; continue; end
                if collapsed(q,s,q2,r(k)+Inter), p(k,1)=Inf; continue; end

                q=[q; p7; q1; q2];
                s=[s r(k) r(k) r(k)];
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

                if collapsed(q,s,q12,r(k)), p(k,1)=Inf; continue; end
                if collapsed(q,s,q13,r(k)), p(k,1)=Inf; continue; end
                if collapsed(q,s,q23,r(k)), p(k,1)=Inf; continue; end
                if collapsed(q,s,q123,r(k)), p(k,1)=Inf; continue; end

                %disp(sprintf('\t %g \t %g \t %g \n',p(k,:),q1,q2,q3,q12,q13,q23,q123))

                q=[q; q1; q2; q3; q12; q13; q23; q123];
                s=[s r(k) r(k) r(k) r(k) r(k) r(k) r(k)];
                num=num+1;
                ind=[ind num num num num num num num];
            end

        else
            num=num+1;
        end

        
        if isempty(progress_str)
            prevLength = 0;
        else
            prevLength = numel(progress_str);
        end
        progress_str = sprintf('Number of cells positioned: %d out of %d tried (max %d)\nPorosity: %.3f\n',numel(unique(ind))- numel(unique(vasc_ind)),k,n, 1 - FracVol(a,q,s,ind));
        if verb,fprintf([repmat('\b',1,prevLength) '%s'],progress_str);end
        drawnow;
        
        %% Ecriture de la (des) nouvelle(s) sphere(s) dans les variables de sortie
        q=[q; p(k,:)];
        s=[s r(k)];
        ind=[ind num];
        %splot(q,s)
        k=k+1;
        lim=0;
        
    end
    if fextra <= (1-FracVol(a,q,s,ind)), disp('WARNING, lower porosity'),end
else
    s = radii;
    q = center;
    ind = numel(radii);
    if verb, disp('No cell positioned');end
end

compart=InCompart(repmat(a,[100 1]).*rand(100,2),q,s*1.3);

%splot(q,s,a)
[dum vect] = unique(ind,'first');

fv = FracVol(a,q,s,ind);
fv_int = FracVol(a,q,r_Intra(ind),ind);
fv_myelin = FracVol(a,q,r_Myelin(ind),ind) - fv_int;

model.Name = 'white_matter';
model.Centers = q;
model.Radius = s;
model.Indices = ind;
model.CubeSize = a;
model.FracVolIntra = fv;
model.FracVolExtra = 1-fv;
model.MainSpheres = vect;
model.Swelled = 'No';
model.Myelin = Myelin;
if strcmpi(Myelin,'yes'),
model.RadiusIntra = r_Intra(ind);
model.RadiusMyelin = r_Myelin(ind);
model.FracVolIntra = fv_int;
model.FracVolMyelin = fv_myelin;
end
% model.ParameterFile = txt;
model.Distribution.radius = Radius;
model.Distribution.centers = Centers;
model.Distribution.RadiusMean = mean(s(vect));
model.Distribution.RadiusStd = std(s(vect));
model.OverlapMatrices = OverlapCompartTable(q,s);

% splot(model)



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
% disp(sprintf('Volume intra cellaulaire: %.3g\nNombre d''itérations: %d',fn,N))

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
