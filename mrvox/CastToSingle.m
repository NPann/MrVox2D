function [G too M S] = CastToSingle(G,too,M,S)

%% Cast to single
M.M0 = single(M.M0);
M.par = single(M.par);
M.per = single(M.per);

if too.DHker ~= -1
    too.FTDHker = single(too.FTDHker);
    too.DHkerWvasc = single(too.DHkerWvasc);
    too.DHkerWcell = single(too.DHkerWcell);
    too.DHkerWees = single(too.DHkerWees);
end
too.rot = single(too.rot);
too.R1dt = single(full(too.R1dt));
too.RFtrans = single(too.RFtrans);

G.vasc.P = single(full(G.vasc.P));
G.ees.P  = single(full(G.ees.P));
G.cell.P  = single(full(G.cell.P));

if isstruct(S)
    S.per = single(S.per);
    S.par = single(S.par);
else
S = single(S);
end