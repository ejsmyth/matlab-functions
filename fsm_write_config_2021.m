%%% NEED TO RECOMPILE FSM2 to get etab, rhob, and rhoc parameters

% Prepares text files for Factorial Snow Model 2 configuration
%
% RELEASE NOTES
%   Version 1.0 written by Eric Smyth (esmyth@esmyth.com) Sep 2018
%   Version 2.0 edited by Eric Smyth Mar 2020 (interact with new FSM2 code
%       released Jan 2019)
%   Version 3.0 edited by Eric Smyth Apr 2020 to fix issue with mismatched
%   snow layers and layer depths (need to specify layer depths explicitly
%       if Nsmax ~= 3)
%   Version 4.0 edited by Eric Smyth Jan 2021 to interact with new FSM2
%       code released Nov 2020
%
% SYNTAX
%   fsm_write_config(SETTINGS, GRID, DRIVE, PARAMS, SITE)
%
% INPUTS
%   SETTINGS = structure with the following fields:
%       file_path = where config file will be written
%       met_file = name of met forcing file
%       IC = choice of whether to use an initial condition file (1 = yes)
%       start_file = name of initial condition file (if IC=1)
%       runid = run identifier string
%       %%% not currently working --- > dump_file = name of output state variable dump file
%   GRID = structure that needs to include "Nsmax" and any other of the following
%       Nsmax = maximum number of snow layers (default 3)
%       ------ If Nsmax does not equal 3 exactly --------------------------
%       Dzsnow = snow layer thicknesses, in meters. MUST MATCH Nsmax # of layers
%       -------------------------------------------------------------------
%       Nsoil = number of soil layers (default 4)
%             NOTE: This script only built to work for single cell,
%             Nx=Ny1=1
%       Ncols = number of columns in grid
%       Nrows = number of rows in grid
%   DRIVE = structure with at least one of the following fields and units for data:
%       dt = timestep (s)
%       zT = temp & humidity measurement height (m)
%       zU = wind speed measurement height (m)
%       zsub = subcanopy wind speed diagnostic height (m)
%       ------ If SW partitioning option SWPART = 1 -----------------------
%       noon = time of solar noon
%       lat = Latitude 
%   PARAMS = structure with at least one of the following fields and units for data:
%       asmx = maximum albedo for fresh snow
%       asmn = minimum albedo for melting snow
%       hfsn = snow cover fraction depth scale (m)
%       z0sn = snow surface roughness length (m)
%       rgr0 = fresh snow grain radius (m)
%       ------ If snow albedo option ALBEDO = 1 ---------------------------
%       Talb = snow albedo decay temperature threshold (C)  
%       ------ If snow albedo option ALBEDO = 2 ---------------------------
%       Salb = snowfall to refresh albedo (kg m^-2)
%       tcld = cold snow albedo decay time scale (h)
%       tmlt = melting snow albedo decay time scale (h)
%       ------ If therm conductivity option CONDCT = 0 --------------------
%       kfix = fixed thermal conductivity (W m^-1 K^-1)
%       ------ If snow density option DENSTY = 0 -------------------------
%       rfix = fixed snow density (kg m^-3)
%       ------ If snow density option DENSTY = 1 -------------------------
%       rcld = maximum density for cold snow (kg m^-3)
%       rhof = fresh snow density (kg m^-3)
%       rmlt = maximum density for melting snow (kg m^-3)
%       trho = snow compaction time scale (h)
%       ------ If snow density option DENSITY = 2 -------------------------
%       eta0 = reference snow viscosity (Pa s)
%       etab = snow viscosity parameter (m^3 kg^-1)%%%%%%%%%%%%%%%%%%%%%%%%%%
%       snda = snow densification parameter (s^-1)
%       ------ If snow density option DENSITY = 1 or 2 --------------------
%       rhow = wind-packed snow density kg m^-3
%       rhob = Temperature factor in fresh snow density (kg/m^3/K) %%%%%%%%%%%%
%       rhoc = Wind factor in fresh snow density (kg s^0.5/m^3.5) %%%%%%%%%%%%
%       ------ Snow hydrology option HYDROL = 1,2 ----------__-------------
%       Wirr = irreducible liquid water content of snow
%       ------ Forest canopy option CANMOD = 1 ----------------------------
%       fvg1 = fraction of vegetation in upper canopy layer
%       ------ Canopy parameters ------------------------------------------
%       cvai = vegetation heat capacity per unit VAI (J kg^-1 m^-2)
%       gsnf = snow-free vegetation moisture conductance (m/s)
%       hbas = canopy base height (m)
%       kext = canopy light extinction coefficient
%       rveg = leaf boundary resistance (s^1/2 m^-1/2)
%       svai = intercepted snow capacity per unit VAI (kg m^-2)
%       tunl = canopy snow unloading time scale (h)
%       wcan = canopy wind decay coefficient
%       ------ Forest radiation option CANRAD = 1 -------------------------
%       acn0 = snow-free dense canopy albedo
%       avgs = snow-covered dense canopy albedo
%       ------ Forest radiation option CANRAD = 2 -------------------------
%       avg0 = snow-free vegetation albedo
%       avgs = snow-covered vegetation albedo
%       
%   SITE = structure with at least one of the following fields and units for data:
%       alb0 = snow-free ground albedo
%       vegh = canopy height (m)
%       fsky = sky view fraction
%       VAI = vegetation area index


function fsm_write_config_2021(SETTINGS, GRID, DRIVE, PARAMS, SITE)

%% Build libraries of possible options
lib_gridpnts = {'Nsmax','Nsoil','Ncols','Nrows'};
lib_drive = {'dt', 'zT', 'zU', 'noon', 'lat'};
lib_params = {'asmn', 'asmx', 'eta0', 'hfsn', 'kfix', 'rcld', 'rfix', 'rgr0', 'rhof', 'rhow', ...
    'rmlt', 'Salb', 'snda', 'Talb', 'tcld', 'tmlt', 'trho', 'Wirr', 'z0sn', ...
    'acn0', 'acns', 'avg0', 'avgs', 'cvai', 'gsnf', 'hbas', 'kext', 'rveg', 'rveg', 'svai', 'tunl', ...
    'wcan', 'fcly', 'fsnd', 'gsat', 'z0sf'}; %%% Need to add etab
lib_site = {'alb0', 'vegh', 'fsky', 'VAI'};

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
    if max(strcmp(xgridName,lib_gridpnts))==1
        % then this is valid. write in config file
        eval(['xgrid = GRID.' xgridName ';']);
        fprintf(fid, '%s\n', ['  ' xgridName ' = ' num2str(xgrid)]);
    end
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
if isfield(PARAMS,'fvg1') == 1
    fprintf(fid, '%s\n', ['  fvg1 = ' num2str(PARAMS.fvg1)]);
end
if isfield(DRIVE,'zsub') == 1
    fprintf(fid, '%s\n', ['  zsub = ' num2str(DRIVE.zsub)]);
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

fprintf(fid, '%s\n', '&veg');
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
fprintf(fid, '%s\n', ['  runid = ' num2str(SETTINGS.runid)]);
% fprintf(fid, '%s\n', ['  dump_file = ''' SETTINGS.dump_file '''']);
fprintf(fid, '%s\n', '/');

fclose(fid);