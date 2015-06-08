function [M, Seqflag, too] = Seq_PulseSRefoc(tt,dt,para,M,too)

% Simulate a Diffusion sequence with Narrow Pulse Single Refoc
%
% Code by N.Pannetier, Feb 13rd 2013.


Tacq    = para.Tacq;
TR      = para.TR;
TRF     = para.RF.exc.time;
ang     = para.RF.exc.ang*pi/180;

Td      = para.Diff.Td;
G       = too.Diff.G;
tg      = para.Diff.Grad.t;
pha     = para.RF.exc.pha;

Gc = 0;
RFc = 0;

%% Narrow Pulse
if abs(TRF(2) - tt - Td/2)  < tg*(1-1e-9)
    M.per = M.per.*too.DiffGradrot;
    Gc = G;
end

%% RF pulse
if any( abs(TRF - mod(tt,TR))< dt *(1-1e-9))
    
    ind = find(abs(TRF - mod(tt,TR))< dt *(1-1e-9));
    RFc = ang(ind);
    too.phac = pha(ind);
    
    Mtmp = reshape(cat(3,real(M.per(:)),imag(M.per(:)),M.par(:)),[numel(M.per(:)) 3]) * squeeze(too.RFrotation(:,:,ind));
    
    Mtmp = reshape(Mtmp,[size(M.per,1) size(M.per,2) 3]);
    
    M.par = Mtmp(:,:,3);
    M.per = Mtmp(:,:,1) + 1i * Mtmp(:,:,2);
end

%% Narrow Pulse
if abs(TRF(2) - tt + Td/2)  < tg*(1-1e-9)
    M.per = M.per.*too.DiffGradrot;
    Gc = G;
end

too.seq.Gx  = [too.seq.Gx Gc];
too.seq.Gy  = [too.seq.Gy Gc];
too.seq.RF = [too.seq.RF RFc];
too.seq.RFphi   = [too.seq.RFphi too.phac];

%% Test pour acquisition
if any( abs(Tacq - mod(tt,TR))< dt *(1-1e-9))
    Seqflag.Acq = 1;
else
    Seqflag.Acq = 0;
end


