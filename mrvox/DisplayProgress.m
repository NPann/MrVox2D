function [progress_str, tic_dt] = DisplayProgress(tt,Model,progress_str,tic_dt,ticI)

% Display algorithm progression

percentFinished = floor(tt/Model.Tmax*100);
elapsed_time = toc(ticI);

if isempty(progress_str)
    prevLength = 0;
    time_left = 0;
    rate = 0;
else
    prevLength = numel(progress_str);
    time_left     = (Model.Tmax/tt-1) * elapsed_time;
    elapsed_time_dt = toc(tic_dt);
    rate = elapsed_time_dt;
end
progress_str = sprintf('%3.0f %% simulated in %3.0f s, at %3.2f s/dt (%3.0f /%3.0f s done)\n Estimated time left: %3.0f s \n',...
    percentFinished,elapsed_time,rate,tt,Model.Tmax,time_left);
fprintf([repmat('\b',1,prevLength) '%s'],progress_str);
tic_dt=tic;