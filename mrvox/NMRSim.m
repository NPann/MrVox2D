function [M, Seqflag, too] = NMRSim(Model,too,G,M,Pulseq,Seq,tt)

% Compute the evolution of the magnetization matrices


%% Evolution
M.per = M.per .* too.rot;
M.par = (M.par - M.M0) .* exp(-too.R1dt) + M.M0; 

%% Compute diffusion of water
if Model.phy.DH2O > 0
    if Model.Flag.HinderedDiff % Hindered diffusion
        M.per = (ifftn(fftn(M.per.*G.cell.P) .* too.FTDHker) + M.per .* full(G.cell.P).*too.DHkerWcell) .* full(G.cell.P) + ...
                (ifftn(fftn(M.per.*G.vasc.P) .* too.FTDHker) + M.per .* full(G.vasc.P).*too.DHkerWvasc) .* full(G.vasc.P) + ...
                (ifftn(fftn(M.per.*G.ees.P)  .* too.FTDHker) + M.per .* full(G.ees.P) .*too.DHkerWees)  .* full(G.ees.P);
        M.par = (ifftn(fftn(M.par.*G.cell.P) .* too.FTDHker) + M.par .* full(G.cell.P).*too.DHkerWcell) .* full(G.cell.P) + ...
                (ifftn(fftn(M.par.*G.vasc.P) .* too.FTDHker) + M.par .* full(G.vasc.P).*too.DHkerWvasc) .* full(G.vasc.P) + ...
                (ifftn(fftn(M.par.*G.ees.P)  .* too.FTDHker) + M.par .* full(G.ees.P) .*too.DHkerWees)  .* full(G.ees.P);
    else
         M.per = ifftn(fftn(M.per) .* too.FTDHker);
         M.par = ifftn(fftn(M.par) .* too.FTDHker);
    end
end


%% Play MR sequence 
[M, Seqflag,too] = Pulseq(tt,Model.dt,Seq,M,too);



