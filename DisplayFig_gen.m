function FigProp = DisplayFig_gen(G,M,too,S,t)

%% Magnetization lattices
if numel(M.par) > 1
    CastBox = contourc(full(double(imfill(too.VisuInd,'holes'))),1);
    

    f1 = figure('Color',[1 1 1],'Position',[10 500 600 500]);
    VessContour = contourc(full(double(G.vasc.P)),1);
    subplot(2,2,1),
    h1 = imagesc(abs(M.per));axis square, axis off,title('abs(Mper)'),colorbar%;colormap gray
    inds = 1;
    hold on
    subplot(2,2,1),plot(CastBox(1,2:end),CastBox(2,2:end))
    while inds < size(VessContour,2)
        subplot(2,2,1),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
        inds = inds + VessContour(2,inds)+1;
    end
    hold off
    
    subplot(2,2,2),
    h2 = imagesc(angle(M.per));axis square, axis off,title('Arg(Mper)'),colorbar
    inds = 1;
    hold on,
    subplot(2,2,2),plot(CastBox(1,2:end),CastBox(2,2:end))
    while inds < size(VessContour,2)
        subplot(2,2,2),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
        inds = inds + VessContour(2,inds)+1;
    end
    hold off
    
    subplot(2,2,3),
    tmp = too.dBm;
    % tmp(logical(G.vasc.P)) = NaN;
    h3 = imagesc(tmp);axis square, axis off,title('B (T)')
    inds = 1;
    hold on,
    subplot(2,2,3),plot(CastBox(1,2:end),CastBox(2,2:end)),colorbar
    while inds < size(VessContour,2)
        subplot(2,2,3),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
        inds = inds + VessContour(2,inds)+1;
    end
    hold off
    
    
    subplot(2,2,4),
    h4 = imagesc((M.par));axis square, axis off,title('Mpar)'),colorbar
    inds = 1;
    hold on,
    subplot(2,2,2),plot(CastBox(1,2:end),CastBox(2,2:end))
    while inds < size(VessContour,2)
        subplot(2,2,4),plot(VessContour(1,inds+1:inds+VessContour(2,inds)),VessContour(2,inds+1:inds+VessContour(2,inds)),'color',[0 0 0]); hold on,
        inds = inds + VessContour(2,inds)+1;
    end
    hold off
    
    
    FigProp.f1 = f1;
    FigProp.h1 =h1;
    FigProp.h2 =h2;
    FigProp.h3 =h3;
    FigProp.h4 =h4;
end

%% Signal + Seq
f2 = figure('Color',[1 1 1],'Position',[620 500 600 500]);
subplot(2,3,1),
p1 = plot((0:numel(too.seq.Gx)-1)*too.dt,too.seq.Gx);ylim([1.2*min(too.seq.Gx(:))-1 1.2*max(too.seq.Gx(:))+1]),ylabel('G (mT/m)'),xlabel('t (ms)'),hold on
p12 = plot((0:numel(too.seq.Gy)-1)*too.dt,too.seq.Gy,'r');ylim([1.2*min(too.seq.Gy(:))-1 1.2*max(too.seq.Gy(:))+1]),ylabel('G (mT/m)'),xlabel('t (ms)'),hold off,title('Gradients')
subplot(2,3,2),
% p2 = plot((0:numel(too.seq.Gx)-1)*too.dt,180/pi*too.seq.RF);ylabel('Angle'),xlabel('t (ms)'),title('RF'),hold on
p2 = plotyy((0:numel(too.seq.Gx)-1)*too.dt,180/pi*too.seq.RF,(0:numel(too.seq.Gx)-1)*too.dt,180/pi*too.seq.RFphi);
xlabel('t (ms)'),title('RF')%,hold on
set(get(p2(1),'Ylabel'),'String','Angle')
set(get(p2(2),'Ylabel'),'String','Phase');
set(p2(1),'YLimMode','auto')
set(p2(2),'YLimMode','auto')
set(p2(1),'YTickMode','auto')
set(p2(2),'YTickMode','auto')
set(p2(1),'XLimMode','auto')
set(p2(2),'XLimMode','auto')
set(p2(1),'XTickMode','auto')
set(p2(2),'XTickMode','auto')
% p22 = plot((0:numel(too.seq.Gx)-1)*too.dt,180/pi*too.seq.RFphi);ylabel('Angle'),xlabel('t (ms)'),title('RF'),hold off
subplot(2,3,4),
p4 = plot(t*1e3,abs(S),'.');ylabel('Abs'),xlabel('t (ms)'),title('abs(S)')
subplot(2,3,5),
p5 = plot(t*1e3,angle(S),'.');ylabel('Angle'),xlabel('t (ms)'),title('angle(S)')
subplot(2,3,6),
NbPoints = 500;
Range = 1:ceil(numel(M.per(:))/NbPoints):numel(M.per(:));
p6 = plot(real(M.per(Range)),imag(M.per(Range)),'.');ylabel('image(Mper)'),xlabel('real(Mper)'),title('Transerve Plane')
xlim([-1 1]),ylim([-1 1])

subplot(2,3,3),
p3 = plot3(real(M.per(Range)),imag(M.per(Range)),M.par(Range),'.');
xlim([-1 1]),ylim([-1 1]),zlim([-1 1])


FigProp.f2 = f2;
FigProp.p1 =p1;
FigProp.p12 =p12;
FigProp.p21 =findobj(p2(1),'Type','line');
FigProp.p22 =findobj(p2(2),'Type','line');
% FigProp.p22 =p22;
FigProp.p3 = p3;
FigProp.p4 =p4;
FigProp.p5 =p5;
FigProp.p6 =p6;
FigProp.Range = Range;