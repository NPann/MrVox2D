function b = ComputebValue(FileName,LineNum)


[Model, Seq] = ReadModelAndSeq(FileName);

%% Generate Input Table
[Label, par]    = GenerateInputDico(Model);
if Model.Flag.Verbose, fprintf('Parameter Input: \n');end
if ~isempty(par)
    for a=1:numel(Label)
        eval(sprintf('%s = par(LineNum,%d);',Label{a},a));
        if Model.Flag.Verbose, fprintf('%s = %.2e\n',Label{a},eval(Label{a}));end
    end
end

%% Initiate RandStream
stream = RandStream('mt19937ar','Seed',Model.geo.vasc.Id);
RandStream.setGlobalStream(stream)

%% Other Input parameters
Seq.RF.exc.ang  = Seq.RF.exc.ang(:) * Model.vox.RFfact;
Pulseq = eval(Seq.Name);
Model.phy.gamma = 2.675*10^8;

% Geometry
G = GeoGen(Model);
% Simulation lattices
[too, M] = InitSim(Model,G);
% Sequence
too = InitSeq(Seq,too);

Gs = too.Diff.G;

fprintf('\nVoxel fov: %.1f um\n',too.fov*1e6);
fprintf('\nSeq Type: %s\n',Seq.Id);
fprintf('Grad strength: %.1f mT/m\n',Gs*1e3);
fprintf('Grad Dur: %.1f ms\n',Seq.Diff.Grad.t*1e3);
fprintf('Diff delay: %.1f ms\n\n',Seq.Diff.Td*1e3);


Model.geo.res =10;
Model.geo.N = 10;


%% Others
count = 1;
S = struct('per',[],'par',[]);
t = [];
if Model.Flag.Verbose,ticI = tic;tic_dt=tic;str = '';end
if Model.Flag.Verbose,fprintf('Voxel fov: \t%.1fum\nLattice size: \t%d\nStep time: \t%.1fms\n\n',too.fov*1e6,Model.geo.vasc.N,Model.dt*1e3);end
if Model.Flag.Display,AlreadyDraw = 0;Svisu = [];end

if isfield(Seq,'Imaging')
    Model.Tmax = Seq.Imaging.Ph.step * Seq.TR + Seq.Imaging.Dummy * Seq.TR;
else
    Model.Tmax = Seq.TR-Model.dt;
end

% Single/Double case to speed up
if Model.Flag.SingleCast,[G too M S] = CastToSingle(G,too,M,S);end

% Define summation lattice
if Model.Flag.ExtraVascOnly
    VisuInd = false(size(M.per));
    VisuInd(logical(1-full(G.vasc.P))) = 1;
else
    VisuInd = true(size(M.per));
end
if Model.Flag.CastFact < 1
    Nc = round(sqrt(Model.geo.res.^2*Model.Flag.CastFact));
    VisuCast = false(size(M.per));
    Ncc = round((Model.geo.res-Nc)/2);
    VisuCast(Ncc+1:end-Ncc,Ncc+1:end-Ncc) = 1;
    VisuInd = VisuInd & VisuCast;
end
NormFact = sum(VisuInd(:));
too.VisuInd = VisuInd;

%% TIME LOOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%
too.Diff.G = Gs;
for tt=0:Model.dt:Model.Tmax
    [M, ~,too] = Pulseq(tt,Model.dt,Seq,M,too);
end
time = 0:Model.dt:Model.Tmax;


pol = ones(size(too.seq.Gx));
PiPulse = find(too.seq.RF == pi);
for a=1:numel(PiPulse)
    pol(PiPulse(a):end) = pol(PiPulse(a):end)*(-1);
end

for a = 1:numel(time)
    tmp(a) = sum(1e-3*sqrt(too.seq.Gx(1:a).^2 + too.seq.Gy(1:a).^2).*pol(1:a) * Model.dt).^2;
end
b = (Model.phy.gamma)^2 * sum(tmp*Model.dt)/1e6; % in s/mm2
fprintf('b-value =  %.1f s/mm2\n\n',b);

