% Prepares text files for Factorial Snow Model 2 configuration
%
% RELEASE NOTES
%   Version 1.0 written by Eric Smyth (esmyth@esmyth.com) Sep 2018
%   Version 2.0 edited by Eric Smyth Mar 2020 (interact with new FSM2 code
%   released Jan 2019)
%   Version 3.0 edited by Eric Smyth Apr 2020 to fix issue with mismatched
%   snow layers and layer depths (need to specify layer depths explicitly
%   if Nsmax ~= 3)
%
% SYNTAX
%   fsm_write_config(SETTINGS, GRID, DRIVE, PARAMS, SITE, OUT)
%
% INPUTS
%   SETTINGS = structure with the following fields:
%       file_path = where config file will be written
%       met_file = name of met forcing file
%       IC = choice of whether to use an initial condition file (1 = yes)
%       start_file = name of initial condition file (if IC=1)
%       OPTIONAL: ztop_file = name of DEM file for downscaling (if compiled
%       with DOWNSC=1)
%       %%% not currently working --- > dump_file = name of output state variable dump file
%   GRID = structure that needs to include "Nsmax" and any other of the following
%       Nsmax = maximum number of snow layers (default 3)
%       ------ If Nsmax does not equal 3 exactly --------------------------
%       Dzsnow = snow layer thicknesses, in meters. MUST MATCH Nsmax # of layers
%       -------------------------------------------------------------------
%       Nsoil = number of soil layers (default 4)
%             NOTE: This script only built to work for single cell,
%             Nx=Ny1=1
%       Nx = number of grid points in x direction or in sequence
%       Ny = number of grid points in y direction
%       ztop_file = DEM file name (only if DOWNSC=1 in compile)
%   DRIVE = structure with at least one of the following fields and units for data:
%       dt = timestep (s)
%       zT = temp & humidity measurement height (m)
%       zU = wind speed measurement height (m)
%       ------ If downscaling option DOWNSC = 1 ---------------------------
%       Pscl = precip adjustment scale (km^-1)
%       Tlps = temp lapse rate (K/km)
%       Tsnw = snow threshold temperature (C)
%       zaws = weather station elevation for downscaling (m)
%       ------ If SW partitioning option SWPART = 1 -----------------------
%       noon = time of solar noon
%       lat = Latitude 
%   PARAMS = structure with at least one of the following fields and units for data:
%       asmx = maximum albedo for fresh snow
%       asmn = minimum albedo for melting snow
%       gsat = surface conductance for saturated soil (m/s)
%       hfsn = snow cover fraction depth scale (m)
%       Nitr = number of iterations in energy balance calc
%       z0sn = snow surface roughness length (m)
%       rgr0 = fresh snow grain radius (m)
%       ------ If snow albedo option ALBEDO = 0 ---------------------------
%       Talb = snow albedo decay temperature threshold (C)  
%       ------ If snow albedo option ALBEDO = 1 ---------------------------
%       Salb = snowfall to refresh albedo (kg m^-2)
%       tcld = cold snow albedo decay time scale (h)
%       tmlt = melting snow albedo decay time scale (h)
%       ------ If therm conductivity option CONDCT = 0 --------------------
%       kfix = fixed thermal conductivity (W m^-1 K^-1)
%       ------ If therm conductivity option CONDCT = 1 --------------------
%       bthr = thermal conductivity exponent
%       ------ If snow density option DENSITY = 0 -------------------------
%       rho0 = fixed snow density (kg m^-3)
%       ------ If snow density option DENSITY = 1 -------------------------
%       rcld = maximum density for cold snow (kg m^-3)
%       rmlt = maximum density for melting snow (kg m^-3)
%       trho = snow compaction time scale (h)
%       ------ If snow density option DENSITY = 2 -------------------------
%       eta0 = reference snow viscosity (Pa s)
%       etab = snow viscosity parameter (m^3 kg^-1)
%       snda = snow densification parameter (s^-1)
%       ------ If snow density option DENSITY = 1 or 2 --------------------
%       rhof = fresh snow density (kg m^-3)
%       rhob = Temperature factor in fresh snow density (kg/m^3/K)
%       rhoc = Wind factor in fresh snow density (kg s^0.5/m^3.5)
%       ------ Atmospheric stability option EXCHNG = 1 --------------------
%       bstb = atmospheric stability parameter
%       ------ Bucket hydrology option HYDROL = 1 -------------------------
%       Wirr = irreducible liquid water content of snow
%       ------ Canopy parameters ------------------------------------------
%       avg0 = snow-free vegetation albedo
%       avgs = snow-covered vegetation albedo
%       cden = dense canopy turbulent canopy coefficient
%       cvai = canopy snow capacity per unit VAI (kg m^-2)
%       cveg = vegetation turbulent transfer coefficient
%       Gcn1 = Leaf angle distribution parameter
%       Gcn2 = Leaf angle distribution parameter
%       gsnf = snow-free vegetation moisture conductance (m/s)
%       kdif = diffuse radiation extinction coefficient
%       kveg = canopy cover coefficient
%       rchd = displacement height to canopy height ratio
%       rchz = roughness length to canopy height ratio
%       tcnc = canopy unloading time scale for cold snow (h)
%       tcnm = canopy unloading time scale for melting snow (h)
%   SITE = structure with at least one of the following fields and units for data:
%       alb0 = snow-free ground albedo
%       canh = canopy heat capacity (J K^-1 m^-2)
%       fcly = soil clay fraction
%       fsnd = soil sand fraction
%       fsky = sky view fraction
%       fveg = canopy cover fraction
%       hcan = canopy height (m)
%       scap = Canopy snow capacity (kg/m^2)
%       trcn = canopy transmissivity
%       VAI = vegetation area index
%       z0sf = snow-free ground roughness length
%   OUT = structure with the following fields:
%       Nave = number of timesteps in averaged outputs
%       Nsmp = timestep of sample outputs (<= Nave)
%       runid = run identifier string

function fsm_write_config(SETTINGS, GRID, DRIVE, PARAMS, SITE, OUT)

%% Build libraries of possible options
lib_grid = {'Nsmax','Nsoil','Nx','Ny'};
lib_drive = {'dt', 'zT', 'zU', 'Pscl', 'Tlps', 'Tsnw', 'zaws', 'noon', 'lat'};
lib_params = {'asmx', 'asmn', 'gsat', 'hfsn', 'Nitr', 'z0sn', 'rgr0', 'Talb', 'Salb', ...
    'tcld', 'tmlt', 'kfix', 'bthr', 'rho0', 'rcld', 'rmlt', 'rhof', 'rhob', 'rhoc', 'trho', 'eta0',  ...
    'etab', 'snda', 'bstb', 'Wirr', 'avg0', 'Gcn1', 'Gcn2', ...
    'avgs', 'cden', 'cvai', 'cveg', 'gsnf', 'kdif', 'kveg', 'rchd', 'rchz', 'tcnc', 'tcnm'};
lib_site = {'alb0', 'canh', 'fcly', 'fsnd', 'fsky', 'fveg', 'hcan', 'trcn', 'VAI', 'z0sf', 'scap'};

%% Create config file
fid = fopen(SETTINGS.file_path, 'w');

%% Grid dimensions namelist
% Check for parameters specified
FN = fieldnames(GRID);
ngrid = numel(FN); % how many were specified

fprintf(fid, '%s\n', '&gridpnts');
for hh=1:ngrid
    %%% current param specified
    xgridName = char(FN(hh));
    
    %%% check if valid
    if max(strcmp(xgridName,lib_grid))==1
        % then this is valid. write in config file
        eval(['xgrid = GRID.' xgridName ';']);
        fprintf(fid, '%s\n', ['  ' xgridName ' = ' num2str(xgrid)]);
    end
end
if isfield(SETTINGS,'ztop_file') == 1
    fprintf(fid, '%s\n', ['  ztop_file = ''' SETTINGS.ztop_file '''']);
end
fprintf(fid, '%s\n', '/');

%% Model levels namelist
fprintf(fid, '%s\n', '&gridlevs');
if isfield(GRID,'Dzsnow') == 1
    fprintf(fid, '%s', '  Dzsnow = ');
    for ii = 1:length(GRID.Dzsnow)
        if ii == length(GRID.Dzsnow)
            fprintf(fid, '%s\n', num2str(GRID.Dzsnow(ii)));
        else
            fprintf(fid, '%s', [num2str(GRID.Dzsnow(ii)) ', ']);
        end
    end
else
    % Snow layer thicknesses not specified, check if necessary
    if isfield(GRID,'Nsmax') == 1
        % Check if correct number of layers
        if GRID.Nsmax == 3
            % This is fine, can use default thicknesses
        else
            error('Must specify same number of layers as Dzsnow')
        end
    else
        % Not specifying number of layers, do nothing
    end
end
fprintf(fid, '%s\n', '/');

%% Driving data namelist
% Check for parameters specified
FN = fieldnames(DRIVE);
ndrive = numel(FN); % how many were specified

fprintf(fid, '%s\n', '&drive');
fprintf(fid, '%s\n', ['  met_file = ''' SETTINGS.met_file '''']);
for hh=1:ndrive
    %%% current param specified
    xdriveName = char(FN(hh));
    
    %%% check if valid
    if max(strcmp(xdriveName,lib_drive))==1
        % then this is valid. write in config file
        eval(['xdrive = DRIVE.' xdriveName ';']);
        fprintf(fid, '%s\n', ['  ' xdriveName ' = ' num2str(xdrive)]);
    end
end
fprintf(fid, '%s\n', '/');

%% Parameters namelist
% Check for parameters specified
FN = fieldnames(PARAMS);
nparams = numel(FN); % how many were specified

fprintf(fid, '%s\n', '&params');
for j=1:nparams
    %%% current param specified
    xparamName = char(FN(j));
    
    %%% check if valid
    if max(strcmp(xparamName,lib_params))==1
        % then this is valid. write in config file
        eval(['xparam = PARAMS.' xparamName ';']);
        fprintf(fid, '%s\n', ['  ' xparamName ' = ' num2str(xparam)]);
    end
end
fprintf(fid, '%s\n', '/');

%% Site characteristics namelist
% Check for site characteristics specified
FN = fieldnames(SITE);
nsite = numel(FN); % how many were specified

fprintf(fid, '%s\n', '&maps');
for pp=1:nsite
    %%% current param specified
    xsiteName = char(FN(pp));
    
    %%% check if valid
    if max(strcmp(xsiteName,lib_site))==1
        % then this is valid. write in config file
        eval(['xsite = SITE.' xsiteName ';']);
        fprintf(fid, '%s\n', ['  ' xsiteName ' = ' num2str(xsite)]);
    end
end
fprintf(fid, '%s\n', '/');

%% Initial conditions namelist
fprintf(fid, '%s\n', '&initial');
if SETTINGS.IC == 1
    fprintf(fid, '%s\n', ['  start_file = ''' SETTINGS.start_file '''']);
end
fprintf(fid, '%s\n', '/');

%% Output namelist
fprintf(fid, '%s\n', '&outputs');
fprintf(fid, '%s\n', ['  Nave = ' num2str(OUT.Nave)]);
fprintf(fid, '%s\n', ['  Nsmp = ' num2str(OUT.Nsmp)]);
fprintf(fid, '%s\n', ['  runid = ' num2str(OUT.runid)]);
% fprintf(fid, '%s\n', ['  dump_file = ''' SETTINGS.dump_file '''']);
fprintf(fid, '%s\n', '/');

fclose(fid);