function [M, Seqflag,too] = Seq_TrueFISP(tt,dt,para,M,too)

% Simulate a TrueFISP pulse sequence
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

% NRF with read cst gradient
Tacq = para.Tacq;
TR  = para.TR;
TRF = para.RF.exc.time;
ang = (-1)^(too.inc+1) * para.RF.exc.ang*pi/180;
dur = para.RF.exc.dur;
pha = (para.RF.exc.pha0 + para.RF.exc.phainc * (too.inc * (too.inc-1)/2))*pi/180;

Gxc = 0;
Gyc = 0;

RFc = 0;


%% RF pulse
if any( abs(TRF - mod(tt,TR))< dt *(1-1e-9))
    
    RFc = ang;
    too.phac = pha;
    ra = zeros([too.N*too.N 3]);
    rg = zeros([too.N*too.N 1]);
    RFtrans = zeros([too.N*too.N 9]);
    if too.flag.OffFlipAng
        ra(:,:) = bsxfun(@rdivide,cat(3,ang/(dur*too.gamma)*cos(pha)*ones([too.N*too.N 1]),ang/(dur*too.gamma)*sin(pha)*ones([too.N*too.N 1]),too.dBm(:)), ...
                            sqrt((ang/(dur*too.gamma)*ones([too.N*too.N 1])).^2 + (too.dBm(:)).^2));
        rg = sqrt((too.gamma*too.dBm(:)*dur).^2 + ang^2 * ones([too.N*too.N 1]));
    else
        ra(:,:) = cat(3,cos(pha) * ones([too.N*too.N 1]) , sin(pha)*ones([too.N*too.N 1]),zeros([too.N*too.N 1]));
        rg = ang * ones([too.N*too.N 1]);
    end
    RFtrans(:,1) = ra(:,1).^2 + (1 - ra(:,1).^2) .* cos (rg);
    RFtrans(:,2) = ra(:,1).*ra(:,2) .* (1 -cos(rg)) - ra(:,3) .* sin(rg);
    RFtrans(:,3) = ra(:,1).*ra(:,3) .* (1 -cos(rg)) + ra(:,2) .* sin(rg);
    RFtrans(:,4) = ra(:,1).*ra(:,2) .* (1 -cos(rg)) + ra(:,3) .* sin(rg);
    RFtrans(:,5) = ra(:,2).^2 + (1 - ra(:,2).^2) .* cos (rg);
    RFtrans(:,6) = ra(:,2).*ra(:,3) .* (1 -cos(rg)) - ra(:,1) .* sin(rg);
    RFtrans(:,7) = ra(:,1).*ra(:,3) .* (1 -cos(rg)) - ra(:,2) .* sin(rg);
    RFtrans(:,8) = ra(:,2).*ra(:,3) .* (1 -cos(rg)) + ra(:,1) .* sin(rg); 
    RFtrans(:,9) = ra(:,3).^2 + (1 - ra(:,3).^2) .* cos (rg);
    
    Mtmp = reshape(cat(3,real(M.per(:)),imag(M.per(:)),M.par(:)),[numel(M.per(:)) 3]);
    Mtmp = [RFtrans(:,1) .* Mtmp(:,1) + RFtrans(:,2) .* Mtmp(:,2) + RFtrans(:,3) .* Mtmp(:,3), ...
                 RFtrans(:,4) .* Mtmp(:,1) + RFtrans(:,5) .* Mtmp(:,2) + RFtrans(:,6) .* Mtmp(:,3), ...
                 RFtrans(:,7) .* Mtmp(:,1) + RFtrans(:,8) .* Mtmp(:,2) + RFtrans(:,9) .* Mtmp(:,3)];
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