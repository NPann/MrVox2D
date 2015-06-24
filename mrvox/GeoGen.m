function P=GeoGen(Model)

% Generate the geometry matrices
% Geometry are stored/load in Model.geo.path directory if any.
%
% Input:
%   - geo: structure containing geometry parameters
%
% Output:
%   - P: structure containing geometry matrices



%% Initialize
vasc.P = zeros(Model.geo.res,Model.geo.res);               
cell.P = zeros(Model.geo.res,Model.geo.res);
fov = Model.geo.vasc.R*sqrt(Model.geo.vasc.N*pi/Model.geo.vasc.Vf);
verbose = Model.Flag.Verbose;

%% Generate Geo filename
vasc.filename = sprintf('GeoV_v%d_res%d_R%.1f_N%d_Int%.1f_Vf%.2f', ...
    Model.geo.vasc.Id, ...
    Model.geo.res, ...
    Model.geo.vasc.R * 1e6, ...
    Model.geo.vasc.N, ...
    Model.geo.vasc.Inter * 1e6, ...
    Model.geo.vasc.Vf * 1e2);
cell.filename = sprintf('GeoC_v%d_c%d_res%d_R%.1f_N%d_Int%.1f_Vp%.2f_P%.2f_Rc%.1f_CInter%.1f', ...
    Model.geo.vasc.Id, ...
    Model.geo.cell.Id, ...
    Model.geo.res, ...
    Model.geo.vasc.R * 1e6, ...
    Model.geo.vasc.N, ...
    Model.geo.vasc.Inter * 1e6, ...
    Model.geo.vasc.Vf * 1e2, ...
    Model.geo.cell.Poro * 1e2,...
    Model.geo.cell.Rc * 1e6, ...
    Model.geo.cell.Inter * 1e6);
ind = regexp(vasc.filename,'[.]');
vasc.filename(ind) = 'p';
vasc.filename = fullfile(Model.geo.path,vasc.filename);
vasc.filename = [vasc.filename '.mat'];
ind = regexp(cell.filename,'[.]');
cell.filename(ind) = 'p';
cell.filename = fullfile(Model.geo.path,cell.filename);
cell.filename = [cell.filename '.mat'];

%% Geo already generated ?
if(exist(vasc.filename,'file')~=0)        % for the vessels ?
    if Model.Flag.Verbose,disp('Loading vessel geometry...');end
    load(vasc.filename);
    P.vasc = vasc;
    if(exist(cell.filename,'file')~=0)    % for the cells ?
        if Model.Flag.Verbose,disp('Loading cell geometry...');end
        load(cell.filename);
        P.cell = cell;
    else  % Generate the cells                               
        if Model.Flag.Verbose,disp('Generating cell geometry...');end
        model = CellGenerator(Model.geo.cell.Rc*1e6,Model.geo.cell.Inter*1e6,Model.geo.cell.Poro,vasc.radii*1e6,vasc.center*1e6,fov*1e6,vasc.indices,verbose);
        cell.radii = model.Radius*1e-6;
        cell.center = model.Centers*1e-6;
        
        % Define Lattice
        dx=fov/Model.geo.res;                  
        dy=fov/Model.geo.res;                  

        [x,y] = ndgrid(1:Model.geo.res,1:Model.geo.res);  
        x = dx * x;
        y = dy * y;
        xy = [x(:)  y(:)]';  
        Tmp2 = zeros(Model.geo.res,Model.geo.res);
        
        % Position cells
        for NbEl = 1:numel(cell.radii)
            Tmp = zeros(Model.geo.res,Model.geo.res);
            a = cell.center(NbEl,1);                
            b = cell.center(NbEl,2);
            R = cell.radii(NbEl);
            
            r=sqrt((xy(1,:) - a) .^2+(xy(2,:) -b) .^2);
            indx_inside = r<R;  
            Tmp(indx_inside)= 1;
            
            Tmp2 = Tmp + Tmp2;
        end
        
        clear x y r indx_inside
        
        cell.P = sparse(Tmp2-vasc.P);
        if Model.Flag.SaveGeo, save(cell.filename,'cell');end
        P.cell = cell;
    end
    P.ees.P = sparse(ones(size(P.vasc.P)) - full(P.vasc.P + P.cell.P));
else
    if Model.Flag.Verbose,disp('Generating vessel geometry...'); end 
    if(Model.geo.vasc.Inter ~= -1)
        model = VesselGenerator(Model.geo.vasc.Vf,Model.geo.vasc.R*1e6,Model.geo.vasc.N,Model.geo.vasc.Inter*1e6,verbose);
        vasc.radii = model.Radius*1e-6;
        vasc.center = model.Centers*1e-6;
        vasc.indices = model.Indices;
    elseif(mod(Model.geo.vasc.N^(1/2),1) == 0) % vessels on a cartesian grid
        radii = ones([1 Model.geo.vasc.N])*Model.geo.vasc.R;
        tmp = (fov/sqrt(Model.geo.vasc.N)/2:fov/sqrt(Model.geo.vasc.N):fov);
        center(:,1) = tmp;
        center = repmat(center,[sqrt(Model.geo.vasc.N) 1]);
        for aaa=1:sqrt(Model.geo.vasc.N)
            center((aaa-1)*sqrt(Model.geo.vasc.N)+1:aaa*sqrt(Model.geo.vasc.N),2) = ones([1 sqrt(Model.geo.vasc.N)])*tmp(aaa);
        end
        vasc.radii = radii;
        vasc.center = center;
        vasc.indices = (1:Model.geo.vasc.N);
    else
        if Model.Flag.Verbose,disp('Wrong Rext or Nb cylinder values'),end
        return;
    end
    
     % Define Lattice
    dx=fov/Model.geo.res;                   
    dy=fov/Model.geo.res;                   
    [x,y] = ndgrid(1:Model.geo.res,1:Model.geo.res);
    x = dx * x;
    y = dy * y;
    xy = [x(:)  y(:)]';
    
    % Position vessels
    for it=1:size(vasc.radii,2)
        a = vasc.center(it,1);
        b = vasc.center(it,2);
        Ray = vasc.radii(it);
        r=sqrt((xy(1,:) - a) .^2+(xy(2,:) -b) .^2);
        indx_inside = r<Ray;                        
        vasc.P(indx_inside)= 1;
    end
    
    % Generate cells
    if Model.Flag.Verbose,disp('Generating cell geometry...');end
    model = CellGenerator(Model.geo.cell.Rc*1e6,Model.geo.cell.Inter*1e6,Model.geo.cell.Poro,vasc.radii*1e6,vasc.center*1e6,fov*1e6,vasc.indices,verbose);
    cell.radii = model.Radius*1e-6;
    cell.center = model.Centers*1e-6;
    
    % Position cells
    Tmp2 = zeros(Model.geo.res,Model.geo.res);
    for NbEl = 1:numel(cell.radii)
        Tmp = zeros(Model.geo.res,Model.geo.res);
        a = cell.center(NbEl,1);              
        b = cell.center(NbEl,2);
        R = cell.radii(NbEl);
        
        r=sqrt((xy(1,:) - a) .^2+(xy(2,:) -b) .^2);
        indx_inside = r<R;
        Tmp(indx_inside)= 1;
        
        Tmp2 = Tmp + Tmp2;
    end
    
    clear x y r indx_inside
    
    cell.P = sparse(Tmp2-vasc.P);
    vasc.P = sparse(vasc.P);

    if Model.Flag.SaveGeo,save(vasc.filename,'vasc');end
    if Model.Flag.SaveGeo,save(cell.filename,'cell');end
        
    % Build geometry structure as sparse matrix for saving
    P.vasc = vasc;
    P.cell = cell;
    P.ees.P = sparse(ones(size(P.vasc.P)) - full(P.vasc.P + P.cell.P));
end

