function FigProp = DisplayFig_up(M,too,S,t,FigProp)

if numel(M.par) > 1
    set(FigProp.h1,'CData',abs(M.per))
    set(FigProp.h2,'CData',angle(M.per))
    set(FigProp.h4,'CData',M.par)
end

set(FigProp.p1,'XData',(0:numel(too.seq.Gx)-1)*too.dt*1e3);
set(FigProp.p1,'YData',too.seq.Gx);
set(FigProp.p12,'XData',(0:numel(too.seq.Gy)-1)*too.dt*1e3);
set(FigProp.p12,'YData',too.seq.Gy);
set(get(FigProp.p12,'Parent'),'YLim',[1.2*min([too.seq.Gx(:); too.seq.Gy(:)])-1 1.2*max([too.seq.Gx(:); too.seq.Gy(:)])+1])
set(get(FigProp.p12,'Parent'),'XLim',[0 numel(too.seq.Gy)]*too.dt*1e3)


set(FigProp.p21,'XData',(0:numel(too.seq.Gx)-1)*too.dt*1e3);
set(FigProp.p21,'YData',180/pi*too.seq.RF);
set(FigProp.p22,'XData',(0:numel(too.seq.Gy)-1)*too.dt*1e3);
set(FigProp.p22,'YData',too.seq.RFphi*180/pi);
set(get(FigProp.p21,'Parent'),'XLim',[0 numel(too.seq.Gy)]*too.dt*1e3)
set(get(FigProp.p22,'Parent'),'XLim',[0 numel(too.seq.Gy)]*too.dt*1e3)
if 1.1*min(180/pi*too.seq.RF) <  1.1*max(180/pi*too.seq.RF)
    set(get(FigProp.p21,'Parent'),'YLim',[1.1*min(180/pi*too.seq.RF) 1.1*max(180/pi*too.seq.RF)])
end
% set(FigProp.p2(1),'YLim',[1.2*min(180/pi*too.seq.RF)-1 1.2*max(180/pi*too.seq.RF)+1])

% set(FigProp.p2(1),'XLim',[0 numel(too.seq.Gy)]*too.dt*1e3)
% set(FigProp.p2(2),'YLim',[1.2*min(180/pi*too.seq.RFphi)-1 1.2*max(180/pi*too.seq.RF)+1])
% set(FigProp.p2(2),'XLim',[0 numel(too.seq.Gy)]*too.dt*1e3)

% set(get(FigProp.p2,'Parent'),'XLim',[0 numel(too.seq.Gx)]*too.dt*1e3)

if ~isempty(FigProp.p4)
    
    set(FigProp.p4,'XData',t*1e3);
    set(FigProp.p4,'YData',abs(S));
    set(get(FigProp.p4,'Parent'),'XLim',[0 numel(too.seq.Gx)]*too.dt*1e3)
    
    set(FigProp.p5,'XData',t*1e3);
    set(FigProp.p5,'YData',angle(S));
    set(get(FigProp.p5,'Parent'),'XLim',[0 numel(too.seq.Gx)]*too.dt*1e3)
    
    set(FigProp.p6,'YData',imag(M.per(FigProp.Range)));
    set(FigProp.p6,'XData',real(M.per(FigProp.Range)));
    
    
    set(FigProp.p3,'XData',real(M.per(FigProp.Range))); 
    set(FigProp.p3,'YData',imag(M.per(FigProp.Range)));
    set(FigProp.p3,'ZData',M.par(FigProp.Range));
    
    mx = max([real(M.per(:)); imag(M.per(:)); M.par(:)]);
    mn = min([real(M.per(:)); imag(M.per(:)); M.par(:)]);
    mm = 1.2*max(abs([mn mx]));
    set(get(FigProp.p3,'Parent'),'XLim',[-mm mm],'YLim',[-mm mm],'ZLim',[-mm mm])
    
    mx = max([real(M.per(:)); imag(M.per(:));]);
    mn = min([real(M.per(:)); imag(M.per(:));]);
    mm = 1.2*max(abs([mn mx]));
    set(get(FigProp.p6,'Parent'),'XLim',[-mm mm],'YLim',[-mm mm])
else
    figure(FigProp.f2);
    if ~isempty(S)
        '';
    end
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
    hold on,plot3([-1 1],[0 0],[0 0],'k-'),plot3([0 0],[-1 1],[0 0],'k-'),plot3([0 0],[0 0],[-1 1],'k-'),hold off
    xlim([-1 1]),ylim([-1 1]),zlim([-1 1]),
    set(get(p3,'Parent'),'XGrid','On','YGrid','On','ZGrid','On')
    FigProp.p3 =p3;
    FigProp.p4 =p4;
    FigProp.p5 =p5;
    FigProp.p6 =p6;
    FigProp.Range = Range;
end
drawnow