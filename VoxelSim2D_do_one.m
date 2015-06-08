function [Sa, Sphi] = VoxelSim2D_do_one(ModelFileName,SeqFileName,ModelLineNum,SeqLineNum)
%
%
% This code smulates the MR signal of a 2D plane where vessels and cells are randomly
% positioned.
%
% Input:
%   - ModelFileName:    Path to the file containing the structure of the
%   voxel parameters
%   - SeqFileName:      Path to the file containing the structure of the
%   sequence parameters
%   - LineNum:  Line Index of the set of input parameters that has to be simulated
%
%
% NB: Code can be compiled and be deployed for parallel computation using
% the Matlab RunTime package
%
% N.Pannetier, (Feb 13rd, 2013)


%% Read Input file
if nargin < 2
    fprintf('Error: Missing parameters files\n\n');
    return
elseif nargin > 2
    if isdeployed
        ModelLineNum = str2num(ModelLineNum);
        SeqLineNum   = str2num(SeqLineNum);
    end
else
    ModelLineNum = 1;
    SeqLineNum = 1;
    fprintf('Default parameters values: ModelLine %d, SeqLine %d\n',ModelLineNum,SeqLineNum);
end

Model   = ReadModel(ModelFileName);
Seq     = SelectSeq(SeqFileName,SeqLineNum);

%% Generate Input Table
[ModelLabel, ModelPar]    = GenerateModelTable(Model);
if Model.Flag.Verbose, fprintf('Parameter Input: \n');end
if ~isempty(ModelPar)
    for a=1:numel(ModelLabel)
        eval(sprintf('%s = ModelPar(ModelLineNum,%d);',ModelLabel{a},a));
        if Model.Flag.Verbose, fprintf('%s = %.2e\n',ModelLabel{a},eval(ModelLabel{a}));end
    end
end

%% Initiate RandStream
stream = RandStream('mt19937ar','Seed',Model.geo.vasc.Id);
RandStream.setGlobalStream(stream)

%% Other Input parameters
Seq.RF.exc.ang  = Seq.RF.exc.ang(:) * Model.vox.RFfact;
Pulseq = eval(Seq.Name);
Model.phy.gamma = 2.675*10^8;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%    Simu Starts here    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION
% Geometry
G = GeoGen(Model);

% For the voxel
[too, M] = InitSim(Model,G);
% figure,hist(too.dBm(:),400)

% For the MR sequence
too = InitSeq(Seq,too);

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
%     Model.Tmax = Seq.TR-Model.dt;
      if numel(Seq.TR) > 1
        Model.Tmax = sum(Seq.TR(1:end-1));  
      else
        Model.Tmax = Seq.NbRep * Seq.TR;
      end
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
for tt=0:Model.dt:Model.Tmax
    
    % NMR block
    [M, Seqflag,too] = NMRSim(Model,too,G,M,Pulseq,Seq,tt);
    
    % MR Acquisition
    if Seqflag.Acq
        t(count) = tt;
        S.per(count) = sum(M.per(VisuInd))/NormFact;
        S.par(count) = sum(M.par(VisuInd))/NormFact;
        count = count+1;
    end
    
    
    %% Display real time
    if ~isdeployed
        if Model.Flag.Verbose,[str, tic_dt] = DisplayProgress(tt,Model,str,tic_dt,ticI);end
        if Model.Flag.Display
            if ~AlreadyDraw
                FigProp = DisplayFig_gen(G,M,too,S.per,t);
                AlreadyDraw = 1;
            else
                FigProp = DisplayFig_up(M,too,S.per,t,FigProp);
            end
        end
    end
    
end

Sa   = abs(cast(S.per,'double'));
Sphi = angle(cast(S.per,'double'));

if isdeployed
    fprintf('Sa = \n');
    fprintf('%.16f\n',Sa);
    fprintf('END\n\n\n');
    fprintf('Sphi = \n');
    fprintf('%.16f\n',Sa);
    fprintf('END\n');
end

