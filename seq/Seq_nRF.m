function [M, Seqflag,too] = Seq_nRF(tt,dt,para,M,too)

% Simulate an n RF pulse sequence
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

% NRF with read cst gradient
Tacq = para.Tacq;
TR  = para.TR;
TRF = para.RF.exc.time;
Dur = para.RF.exc.dur;
ang = para.RF.exc.ang*pi/180;
pha = (para.RF.exc.pha0)*pi/180 ;

Gxc = 0;
Gyc = 0;

RFc = 0;

%% RF pulse
if any( abs(TRF+Dur/2 - mod(tt,TR))< dt *(1-1e-9))
    
    ind = find(abs(TRF+Dur/2 - mod(tt,TR))< dt *(1-1e-9));
    
    RFc = ang(ind);
    too.phac = pha(ind);
    
    Mtmp = reshape(cat(3,real(M.per(:)),imag(M.per(:)),M.par(:)),[numel(M.per(:)) 3]);
    Mtmp = [too.RFtrans(:,1,ind) .* Mtmp(:,1) + too.RFtrans(:,2,ind) .* Mtmp(:,2) + too.RFtrans(:,3,ind) .* Mtmp(:,3), ...
                 too.RFtrans(:,4,ind) .* Mtmp(:,1) + too.RFtrans(:,5,ind) .* Mtmp(:,2) + too.RFtrans(:,6,ind) .* Mtmp(:,3), ...
                 too.RFtrans(:,7,ind) .* Mtmp(:,1) + too.RFtrans(:,8,ind) .* Mtmp(:,2) + too.RFtrans(:,9,ind) .* Mtmp(:,3)];
    Mtmp = reshape(Mtmp,[size(M.per,1) size(M.per,2) 3]);
    
    M.par = Mtmp(:,:,3);
    M.per = Mtmp(:,:,1) + 1i * Mtmp(:,:,2);
    too.inc = too.inc+1;
end

%% Test pour acquisition
if any( abs(Tacq - mod(tt,TR))< dt *(1-1e-9))
    Seqflag.Acq = 1;
else
    Seqflag.Acq = 0;
end

too.seq.Gx       = [too.seq.Gx Gxc];
too.seq.Gy       = [too.seq.Gy Gyc];
too.seq.RF      = [too.seq.RF RFc];
too.seq.RFphi   = [too.seq.RFphi mod(too.phac,2*pi)];