function [M, Seqflag,too] = Seq_GREImaging(tt,dt,para,M,too)

% Simulate an n RF pulse sequence
%
% Code by N.Pannetier, C.Debacker, INSERM, Grenoble, 2012

% NRF with read cst gradient
Tacq = (para.TE-too.Seq.HalfTE):(para.TE+too.Seq.HalfTE);
TR = para.TR;
TRF = para.RF.exc.time;
ang = para.RF.exc.ang*pi/180;
TE = para.TE;
ang = para.RF.exc.ang*pi/180;
RFdir = para.RF.exc.dir;
pha = para.RF.exc.pha*pi/180;

Gxc = 0;
Gyc = 0;
RFc = 0;

%% RF pulse
if any( abs(TRF - mod(tt,TR))< dt *(1-1e-9))
    
     ind = find(abs(TRF - mod(tt,TR))< dt *(1-1e-9));
     ang =ang(ind);
     pha = pha(ind);
     
     RFc = ang(ind);
     pha = mod(pha + too.Spoil.phi,2*pi);
     
     if(strcmp(RFdir,'x'))
        impulsion_x=[1 0 0;0 cos(ang) sin(ang);0 -sin(ang) cos(ang)];
        impulsion_y=[cos(pha) 0 -sin(pha);0 1 0;sin(pha) 0 cos(pha)];
        impulsion_z=eye(3);
    elseif(strcmp(RFdir,'y'))
        impulsion_x=eye(3);
        impulsion_y=[cos(ang) 0 -sin(ang);0 1 0;sin(ang) 0 cos(ang)];
        impulsion_z=[cos(pha) sin(pha) 0;-sin(pha) cos(pha) 0;0 0 1];
    elseif(strcmp(RFdir,'z'))
        impulsion_x=[1 0 0;0 cos(pha) sin(pha);0 -sin(pha) cos(pha)];
        impulsion_y=eye(3);
        impulsion_z=[cos(ang) sin(ang) 0;-sin(ang) cos(ang) 0;0 0 1];
    end
    too.RFrotation(:,:) =  impulsion_x * impulsion_y * impulsion_z;

    Mtmp = reshape(cat(3,real(M.per(:)),imag(M.per(:)),M.par(:)),[numel(M.per(:)) 3]) * squeeze(too.RFrotation(:,:));   
    Mtmp = reshape(Mtmp,[size(M.per,1) size(M.per,2) 3]);
    
    M.par = Mtmp(:,:,3);
    M.per = Mtmp(:,:,1) + 1i * Mtmp(:,:,2);
    
    too.Spoil.phi = mod(too.Spoil.phi + para.Imaging.RF.phi*pi/180,2*pi);   % increment RF spoiler
    if too.Dummy >0
        too.Dummy = too.Dummy - 1;
    else
        too.Seq.ph.step = too.Seq.ph.step + 1/(para.Imaging.Ph.step/2);         % phase encoding step
    end
    too.pha = pha;
end

%% Dephasing 
if any( abs(TE - too.Seq.HalfTE - para.Imaging.Def.dt - mod(tt,TR)) < para.Imaging.Def.dt/2 *(1-1e-9))
    M.per = M.per.*too.Grodrot;
    M.per = M.per.*(too.Gphdrot).^too.Seq.ph.step;
    Gxc = too.Grod;
    Gyc = too.Gphdmx*too.Seq.ph.step;
end

%% ADC
if any( abs(TE - mod(tt,TR)) < too.Seq.HalfTE *(1-1e-9))
    M.per = M.per.*too.Grorot;
    Gxc = too.Gro;

    Seqflag.Acq = 1;
else
     Seqflag.Acq = 0;
end

%% Grad Spoiling
if (TR - mod(tt,TR))  <= para.Imaging.Crusher.dt*(1+1e-9)
        M.per = M.per .* too.Spoilrot;
        Gxc = too.Gspoil;
end


too.seq.Gx      = [too.seq.Gx Gxc];
too.seq.Gy      = [too.seq.Gy Gyc];
too.seq.RF      = [too.seq.RF RFc];
too.seq.RFphi   = [too.seq.RFphi too.pha];