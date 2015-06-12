function [too] = InitSeq(para,too)

% Compute RF pulses rotation operator and dynamic gradient


ang = para.RF.exc.ang*pi/180;
dur = para.RF.exc.dur;
pha = (para.RF.exc.pha0)*pi/180;
too.inc = 1;

N = numel(too.dBm(:,1));
ra = zeros([N*N 3]);
rg = zeros([N*N 1]);
too.RFtrans = zeros([N*N 9 numel(ang)]);
for ind=1:numel(para.RF.exc.ang)
    if too.flag.OffFlipAng
        ra(:,:) = bsxfun(@rdivide,cat(3,ang(ind)/(dur(ind)*too.gamma)*cos(pha(ind))*ones([N*N 1]),ang(ind)/(dur(ind)*too.gamma)*sin(pha(ind))*ones([N*N 1]),too.dBm(:)), ...
                            sqrt((ang(ind)/(dur(ind)*too.gamma)*ones([N*N 1])).^2 + (too.dBm(:)).^2));
        rg = sqrt((too.gamma*too.dBm(:)*dur(ind)).^2 + ang(ind)^2 * ones([N*N 1]));
    else
        ra(:,:) = cat(3,cos(pha(ind))*ones([N*N 1]),sin(pha(ind))*ones([N*N 1]),zeros([N*N 1]));
        rg = ang(ind) * ones([N*N 1]);
    end
    
    too.RFtrans(:,1,ind) = ra(:,1).^2 + (1 - ra(:,1).^2) .* cos (rg);
    too.RFtrans(:,2,ind) = ra(:,1).*ra(:,2) .* (1 -cos(rg)) - ra(:,3) .* sin(rg);
    too.RFtrans(:,3,ind) = ra(:,1).*ra(:,3) .* (1 -cos(rg)) + ra(:,2) .* sin(rg);
    too.RFtrans(:,4,ind) = ra(:,1).*ra(:,2) .* (1 -cos(rg)) + ra(:,3) .* sin(rg);
    too.RFtrans(:,5,ind) = ra(:,2).^2 + (1 - ra(:,2).^2) .* cos (rg);
    too.RFtrans(:,6,ind) = ra(:,2).*ra(:,3) .* (1 -cos(rg)) - ra(:,1) .* sin(rg);
    too.RFtrans(:,7,ind) = ra(:,1).*ra(:,3) .* (1 -cos(rg)) - ra(:,2) .* sin(rg);
    too.RFtrans(:,8,ind) = ra(:,2).*ra(:,3) .* (1 -cos(rg)) + ra(:,1) .* sin(rg); 
    too.RFtrans(:,9,ind) = ra(:,3).^2 + (1 - ra(:,3).^2) .* cos (rg);
end

%% Grad Spoiling
if isfield(para,'Spoil')
   maxX = max(too.xB(:));
   dxX = diff(too.xB(1:2));
   G = para.Spoil.Grad.A*2*pi/(2*maxX+dxX);
   too.SpoilGradrot = exp( - 1i* G * (too.xB + maxX));
end

%% Diffusion gradient
if isfield(para,'Diff')
   G = para.Diff.Grad.A*2*pi/(too.fov);
   too.DiffGradrot = exp( - 1i* G * (para.Diff.Grad.dir(1) * too.xB + para.Diff.Grad.dir(2) * too.yB + too.fov/2));
   too.Diff.G = para.Diff.Grad.A*2*pi/((too.fov)*too.dt*too.gamma)*1e3;
end

%% Imaging gradient
if isfield(para,'Imaging')
    too.Dummy = para.Imaging.Dummy;
    too.Seq.ph.step = -1 - 1/(para.Imaging.Ph.step/2);
    too.Spoil.phi = para.Imaging.RF.phi*pi/180;
    too.Seq.HalfTE  = para.Imaging.Ro.dt/2;
    
    too.Gro  = para.Imaging.Ro.A/too.fov/(too.gamma/(2*pi))*1e3;                                        % mT/m
    too.Grod = -para.Imaging.Ro.A*too.Seq.HalfTE/para.Imaging.Def.dt/too.fov/(too.gamma/(2*pi))*1e3;     % mT/m
    too.Gphdmx = para.Imaging.Ro.A*too.Seq.HalfTE/para.Imaging.Def.dt/too.fov/(too.gamma/(2*pi))*1e3;     % mT/m
    
    too.Gspoil = para.Imaging.Crusher.A/too.fov/(too.gamma/(2*pi))*1e3;
    
    
    too.Gphdrot = exp( 1i * too.gamma * too.Gphdmx * 1e-3 * too.dt * too.yB);
    too.Grodrot = exp( 1i * too.gamma * too.Grod   * 1e-3 * too.dt * too.xB);
    too.Grorot  = exp( 1i * too.gamma * too.Gro    * 1e-3 * too.dt * too.xB);
    too.Spoilrot = exp( 1i * too.gamma * too.Gspoil* 1e-3 * too.dt * too.xB);
    
end

