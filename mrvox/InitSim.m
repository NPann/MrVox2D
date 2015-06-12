function [too, M] = InitSim(Model,G)

% Initiate the matrices and structure used in the simulation


%% Initialize variables

fov = Model.geo.vasc.R*sqrt(Model.geo.vasc.N*pi/Model.geo.vasc.Vf);
dx = fov/Model.geo.res;

[x,y]=ndgrid(-Model.geo.res/2:Model.geo.res/2-1,-Model.geo.res/2:Model.geo.res/2-1);
too.x = x;
too.y = y;

[xB,yB]=ndgrid(-fov/2+dx/2:dx:fov/2-dx/2,-fov/2+dx/2:dx:fov/2-dx/2);
too.xB = xB;
too.yB = yB;


%% Water diffusion kernel
if Model.phy.DH2O > 0
    sig = sqrt(2*Model.phy.DH2O*Model.dt);
    sig2 = sig/dx;
    if Model.Flag.ScaleSpaceKer
        too.DHker = exp(-sig2^2).*besseli(sqrt(x.^2+y.^2),sig2^2);  %     Discrete approach
    else
        too.DHker = exp(-((x*dx).^2+(y*dx).^2)/(2*sig^2));          %     Continuous approach
    end
    too.DHker = too.DHker/sum(too.DHker(:));
    too.FTDHker = fftn(fftshift(too.DHker));
    
    %% Bounce on wall.
    %     % Extra
    %     too.DHkerWextra = ifftn(fftn(full(G.vasc.P)) .* too.FTDHker);
    %     too.DHkerWextra(logical(full(G.vasc.P))) = 0;
    %     % Intra
    %     too.DHkerWintra = ifftn(fftn(full(1-G.vasc.P)) .* too.FTDHker);
    %     too.DHkerWintra(logical(full(1-G.vasc.P))) = 0;
    %
    %
    %
    % Vessel
    too.DHkerWvasc = ifftn(fftn(full(1-G.vasc.P)) .* too.FTDHker);
    too.DHkerWvasc(logical(full(1-G.vasc.P))) = 0;
    
    % Cells
    too.DHkerWcell = ifftn(fftn(full(1-G.cell.P)) .* too.FTDHker);
    too.DHkerWcell(logical(full(1-G.cell.P))) = 0;
    
    % Extra cellular space
    too.DHkerWees = ifftn(fftn(full(1-G.ees.P)) .* too.FTDHker);
    too.DHkerWees(logical(full(1-G.ees.P))) = 0;
    
    
else
    too.DHker = -1;
end


% khi
khi =   Model.phy.vasc.khi * full(G.vasc.P) +  ...
    Model.phy.cell.khi * full(G.cell.P) + ...
    Model.phy.ees.khi  * full(G.ees.P) ;

% R2
R2 =    1/Model.phy.vasc.T2 * full(G.vasc.P) + ...
    1/Model.phy.cell.T2 * full(G.cell.P) + ...
    1/Model.phy.ees.T2  * full(G.ees.P)  ;
% R1
R1 =    1/Model.phy.vasc.T1 * full(G.vasc.P) + ...
    1/Model.phy.cell.T1 * full(G.cell.P) + ...
    1/Model.phy.ees.T1  * full(G.ees.P)  ;

% Compute magnetic field perturbations dB
if Model.geo.res > 1 
    if Model.Flag.B0Ori3D == 1
        thetaB0 = [pi/2 pi/2 0];
        [b1, m1, kk1] = unique(thetaB0, 'first');
        for p=1:numel(b1)
            warning('off');
            FtB(:,:) = Model.phy.B0 * fftshift(fftn(khi)) .*(1/3 - ((too.yB * sin(thetaB0(m1(p)))).^2) ./ (too.xB.^2+too.yB.^2));
            warning('on');
            B = real(ifftn(ifftshift(FtB)));
            dB(:,:,p)=B;
        end
        dBm = mean(dB(:,:,kk1),3);
    else
        thetaB0 = Model.vox.B0theta;
        Phi = Model.vox.B0phi;
        warning('off');
        
        FtB = Model.phy.B0 * fftshift(fftn(khi)) .*(1/3 - (( (too.yB* cos(Phi) + too.xB * sin(Phi))* sin(thetaB0)).^2) ./ (too.xB.^2+too.yB.^2));
        warning('on');
        B = real(ifftn(ifftshift(FtB)));
        dBm=B;
    end
else %not a lattice
    dBm = 0;
end

too.dBm = dBm;

%% Add B0 offset
too.dBm = too.dBm + (2*pi)/Model.phy.gamma * Model.vox.B0off;

%% Add extra static gradient
if Model.geo.res > 1
    maxX = max(too.xB(:));
    dxX = diff(too.xB(1:2));
    Gx = Model.vox.B0Gr.gx*2*pi/(2*maxX+dxX)/Model.phy.gamma/Model.dt;
    % Aliasing from diffusion kernel with FFT contrains gradient to
    % produces dephasing modulo 2pi over voxel extend
    too.dBm = too.dBm + Gx * too.xB ;
end

% Compute rotation and lattice-relaxation matrices
too.rot = exp( (-1i* Model.phy.gamma *too.dBm - R2) * Model.dt);
too.R1dt = R1 * Model.dt;

%% Initialization of matrices
C = zeros(size(G.vasc.P));
M.M0 = Model.phy.M0.vasc * full(G.vasc.P) + Model.phy.M0.ees * full(G.cell.P) + Model.phy.M0.ees * full(G.ees.P);
M.par = M.M0;
M.per = zeros(size(G.vasc.P));

too.vf = Model.geo.vasc.Vf;
too.gamma = Model.phy.gamma;
too.dt = Model.dt;
too.fov = fov;
too.dx  = dx;
too.N = size(M.par,1);

too.seq.Gx = [];
too.seq.Gy = [];
too.seq.RF = [];
too.seq.RFphi = [];


too.flag.OffFlipAng = Model.Flag.OffFlipAng;

too.phac = 0;

