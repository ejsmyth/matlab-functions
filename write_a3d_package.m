% writes package to run alpine 3d
%
% RELEASE NOTES
%   Written by Mark Raleigh (mraleig1@gmail.com), Sept 2018
%   Updated by Eric Smyth (esmyth@esmyth.com), Oct 2020
%        Added option to scale A3D precip by a gridded snow depth file (see
%        SGRIDS, A3D_SETUP, and OPTS)
%
% SYNTAX
%    write_a3d_package(METEO, SNOWFILES, SGRIDS, A3D_SETUP, OPTS)
%
% INPUTS
%   METEO = structure for met stations data/metadata with following fields
%       TIME = Lx7 time matrix
%       STA_ID = 1xS cellstr of station ID (no spaces, don't start with number)
%       STA_NAME = 1xS cellstr of station names
%       STA_LAT = 1xS array of station latitudes
%       STA_LON = 1xS array of station longitudes
%       STA_ELEV = 1xS array of station elevations (m)
%       Tair = LxS matrix of air temeperature
%       P = LxS matrix of precipitation
%       Us = LxS matrix of wind speed
%       Ud = LxS matrix of wind direction
%       RH = LxS matrix of relative humidity
%       Qsi = LxS matrix of incoming shortwave
%       Qli = LxS matrix of incoming longwave
%
%   SNOWFILES = structure for snowfiles with following fields
%       SoilAlbedo = MxN matrix or single value (if constant in space) of soil albedo (optional)
%       BareSoil_z0 = MxN matrix or single value (if constant in space) of base soil surf. roughness (optional)
%       CanopyHeight = MxN matrix of canopy height (m) (optional)
%       CanopyLeafAreaIndex = MxN matrix of canopy LAI (optional)
%       CanopyDirectThroughfall = MxN matrix of canopy direct throughfall (i.e. 1-fveg) (optional)
%
%   SGRIDS = structure for land surface data with following fields
%       XX = MxN matrix of x-coordinates (UTM)
%       YY = MxN matrix of y-coordinates (UTM)
%       ZZ = MxN matrix of elev (m) ... DEM
%       LU = MxN matrix of land use map (specify in OPTS structure the type)
%       POI = Qx2 matrix of points of interest (1st col = X-coord (UTM),
%               2nd col = Y-coord (UTM)... will snap to grid and find elev.
%               Enter [] or omit if none desired
%               Enter 1 for all points in domain (just 1x1 element)
%       grid_size = 1x1 value, grid size (m)
%       utmzone = string specifying utmzone (see deg2utm output)
%       timezone = offset from utc (e.g., -7 for MST)
%       SD = MxN matrix of snow depths (m) for ALS scaling (optional, read if OPTS.ALS_SCALING = 1)
%
%   A3D_SETUP = structure for io.ini file with following fields
%       OUTPUTS.GRIDS_DAYS_BETWEEN = number of days between output grids (1=daily,
%                           0.041666 =1h, etc.)
%       OUTPUTS.GRIDS_START = starting hour to increment grid outputs (using
%                       GRIDS_DAY_BETWEEN as interval)
%       OUTPUTS.GRIDS_PARAMETERS = cellstring array or matrix with output grids of
%               interest. possible values:
%           TA 	= Air temperature.
%           RH 	= Relative humidity.
%           VW 	= Wind velocity.
%           ISWR = Incoming short wave radiation.
%           ISWR_DIFF = Incoming short wave, diffuse.
%           ISWR_DIR = Incoming short wave, direct.
%           ILWR = Incoming long wave radiation.
%           HS = Height of snow.
%           PSUM = Water equivalent of precipitations, either solid or liquid.
%           PSUM_PH = Precipitation phase, between 0 (fully solid) and 1 (fully liquid)
%           TSG = Temperature ground surface.
%           TSS = Temperature snow surface.
%           TS0 = Temperature soil surface.
%           TSOIL = Temperature within the soil, at a given depth.
%           SWE = Snow Water Equivalent.
%           RSNO = Snow mean density.
%           TOP_ALB = Albedo from the top (ie above canopy)
%           SURF_ALB = Albedo of the surface (ie below canopy)
%           SP = sphericity
%           RB = bond radius
%           RG = grain radius
%           N3 = grain Coordination number
%           MS_SNOWPACK_RUNOFF = runoff on the surface of the soil (vitual lysimeter)
%           MS_SOIL_RUNOFF = runoff at the bottom of the snow/soil column
%           SFC_SUBL = The mass loss or gain of the top element due to snow (ice) sublimating.
%           STORE = internal usage (precipitation events that are delayed because they are too small)
%           GLACIER = mask showing the glaciated pixels
%           GLACIER_EXPOSED = mask showing the exposed glaciated pixels (ie not snow covered)
%
%               https://models.slf.ch/docserver/alpine3d/html/classSnGrids.html#aa7cc1a3714dd8ca8dc3919695f17649b
%
%       SNOWPACK.ROUGHNESS_LENGTH = aerodnyamic roughness (m)
%       SNOWPACK.HEIGHT_OF_METEO_VALUES = height of T/humidity data (m)
%       SNOWPACK.HEIGHT_OF_WIND_VALUE = height of wind data (m)
%       SNOWPACK.ATMOSPHERIC_STABILITY = string or cellstring specifying
%           atmospheric stability routine. options:
%           MO_MICHLMAYR  - this is a Monin Obukhov (MO) type
%           MO_STEARNS
%           MO_HOLTSLAG
%           MO_LOG_LINEAR
%           NEUTRAL
%           RICHARDSON
%       SNOWPACK.CANOPY = number option (turn canopy on?). values:
%           1=TRUE
%           0=FALSE
%       ALS_DATE = date of depth data used for als scaling, formatted as a 
%           string: 'YYYYMMDD' (optional)
%
%   OPTS = structure for options with following fields
%       EXP_NAME = string to identify this experiment (no spaces, start with letter)
%       LU_TYPE = what format is SGRIDS.LU? (optional)
%                   enter  0 if already in SLF/A3D format
%                   enter  1 if NLCD (default)
%                   enter  2 if LandFire 2014
%       SGRID_SUB = perform spatial subsetting? (optional). if so, this is
%           a 4 element array specifying rectangle where [xmin xmax ymin ymax]
%       TEMP_SUB = perform temporal subsetting? (optional). if so, this is
%           a 2 element array specifying serial date of start/end dates [SD_initial SD_final]
%       SNOFILE_OPTION = what approach for snowfiles / IC / canopy files?
%                   enter 0 if by unique land classes (Default)
%                   enter 1 if by X, Y location
%       FILES_WRITE = 4 element array specifying whether to write (1, default) or not write (0) certain files:
%               [METEO_FILES, SNOWFILES, SGRIDS, IO.INI FILE] (optional)
%       ALS_SCALING = have A3D scale precipitation with a provided snow depth file?
%                   enter 0 if no
%                   enter 1 if yes, also make sure to include the depth
%                   data in SGRIDS and date in A3D_SETUP
%
% OUTPUTS
%
%

function write_a3d_package(METEO, SNOWFILES, SGRIDS, A3D_SETUP, OPTS)

% error('FIX THE run_file to include --restart flag if SNOFILE_OPTION=1')
%% check options

if ispc==1
    sysSlash = '\';
else
    sysSlash = '/';
end

if isfield(SNOWFILES, 'SoilAlbedo')==1
    SoilAlbedo = SNOWFILES.SoilAlbedo;
else
    SoilAlbedo = 0.2;
end

if isfield(SNOWFILES, 'BareSoil_z0')==1
    BareSoil_z0 = SNOWFILES.BareSoil_z0;
else
    BareSoil_z0 = 0.02;
end

CanopyHeight = [];
if isfield(SNOWFILES, 'CanopyHeight')==1
    if size(SNOWFILES.CanopyHeight,1) == size(SGRIDS.XX,1) && size(SNOWFILES.CanopyHeight,2) == size(SGRIDS.XX,2)
        CanopyHeight = SNOWFILES.CanopyHeight;
    end    
end

CanopyLeafAreaIndex = [];
if isfield(SNOWFILES, 'CanopyLeafAreaIndex')==1
    if size(SNOWFILES.CanopyLeafAreaIndex,1) == size(SGRIDS.XX,1) && size(SNOWFILES.CanopyLeafAreaIndex,2) == size(SGRIDS.XX,2)
        CanopyLeafAreaIndex = SNOWFILES.CanopyLeafAreaIndex;
        
        %%% checks on LAI
        a = find(isinf(CanopyLeafAreaIndex)==1);
        if isempty(a) == 0
            error('Found infinite values of LAI!')
        end
    end    
end

CanopyDirectThroughfall = [];
if isfield(SNOWFILES, 'CanopyDirectThroughfall')==1
    if size(SNOWFILES.CanopyDirectThroughfall,1) == size(SGRIDS.XX,1) && size(SNOWFILES.CanopyDirectThroughfall,2) == size(SGRIDS.XX,2)
        CanopyDirectThroughfall = SNOWFILES.CanopyDirectThroughfall;
        CanopyDirectThroughfall(CanopyDirectThroughfall==1) = 0.99;
        CanopyDirectThroughfall(CanopyDirectThroughfall==0) = 0.01;
    end    
end


if isfield(OPTS, 'A3D_ROOT')~=1
    %%% no root path provided, assume as pwd
    OPTS.A3D_ROOT = [pwd sysSlash];
else
    OPTS.A3D_ROOT = char(OPTS.A3D_ROOT);
    
    %%% check if slash exists at end
    if OPTS.A3D_ROOT(end)~=sysSlash
        OPTS.A3D_ROOT = [OPTS.A3D_ROOT sysSlash];
    end
end

if isfield(OPTS, 'SGRID_SUB')==1
    if numel(OPTS.SGRID_SUB)~=4 && isempty(OPTS.SGRID_SUB)~=1
        error('OPTS.SGRID_SUB must be 4 element array specifying rectangular box, or set to empty')
    else
        sXmin = nanmin(OPTS.SGRID_SUB(1:2));
        sXmax = nanmax(OPTS.SGRID_SUB(1:2));
        sYmin = nanmin(OPTS.SGRID_SUB(3:4));
        sYmax = nanmax(OPTS.SGRID_SUB(3:4));
    end
else
    OPTS.SGRID_SUB = [];
end

if isfield(OPTS, 'TEMP_SUB')==1
    if numel(OPTS.TEMP_SUB)~=2 && isempty(OPTS.TEMP_SUB)~=1
        error('OPTS.TEMP_SUB must be 2 element array specifying start/end serial dates, or set to empty')
    else
        SD_i = OPTS.TEMP_SUB(1);
        SD_f = OPTS.TEMP_SUB(2);
    end
else
    SD_i = nanmin(METEO.TIME(:,7));
    SD_f = nanmax(METEO.TIME(:,7));
end

if isfield(OPTS, 'FILES_WRITE')==1
    if numel(OPTS.FILES_WRITE)~=4 && isempty(OPTS.FILES_WRITE)~=1
        error('OPTS.FILES_WRITE must be 4 element array specifying flags for files to write, or set to empty')
    end
else
    OPTS.FILES_WRITE = [1 1 1 1];
end

if isfield(OPTS, 'SNOFILE_OPTION')~=1
    OPTS.SNOFILE_OPTION = 0;
end

if isfield(A3D_SETUP.OUTPUTS, 'GRIDS_PARAMETERS')==1
    flag_gridParams = 1;
else
    flag_gridParams = 0;
end

tb = sprintf('\t');

%% setup folders

%%% input directory
PATH_INPUT = [OPTS.A3D_ROOT 'input' sysSlash];
PATH_MET = [PATH_INPUT 'meteo' sysSlash];
PATH_SNOW = [PATH_INPUT 'snowfiles' sysSlash];
PATH_GRID = [PATH_INPUT 'surface-grids' sysSlash];

if exist(PATH_INPUT, 'dir')==7
    rmdir(PATH_INPUT, 's')
end

mkdir(PATH_INPUT);
mkdir(PATH_MET);
mkdir(PATH_SNOW);
mkdir(PATH_GRID);

%%% output directory

%%% setup directory


%%

%%% repeat variables
nodata = -9999;

%%% epsg
utmzone_num = find(SGRIDS.utmzone==' ',1,'first');
utmzone_num = SGRIDS.utmzone(1:utmzone_num-1);
utmzone_num = str2double(utmzone_num);
utmzoneCompact = SGRIDS.utmzone;
utmzoneCompact(isspace(utmzoneCompact)) = [];
epsg = 32600 + utmzone_num; % 32613 is UTM Zone 13 (32600+zone for Northern Hemisphere)


%% write surface-grids

%%% perform spatial subsetting?
if isempty(OPTS.SGRID_SUB)==0
    [r1, c1] = nearest2D_mat(sXmin, sYmin, SGRIDS.XX, SGRIDS.YY, 1);
    [r2, c2] = nearest2D_mat(sXmax, sYmax, SGRIDS.XX, SGRIDS.YY, 1);
    
    if r1>r2
        r3 = r1;
        r1 = r2;
        r2 = r3;
    end
    
    if c1>c2
        c3 = c1;
        c1 = c2;
        c2 = c3;
    end
    
    [sRows, sCols] = size(SGRIDS.XX);
    
    SGRIDS.XX = SGRIDS.XX(r1:r2,c1:c2);
    SGRIDS.YY = SGRIDS.YY(r1:r2,c1:c2);
    SGRIDS.ZZ = SGRIDS.ZZ(r1:r2,c1:c2);
    SGRIDS.LU = SGRIDS.LU(r1:r2,c1:c2);
    SGRIDS.SD = SGRIDS.SD(r1:r2,c1:c2);
    
    %%% apply subsetting to other inputs if same size
    if size(SoilAlbedo,1)==sRows && size(SoilAlbedo,2)==sCols
        SoilAlbedo = SoilAlbedo(r1:r2,c1:c2);
    end
    
    if size(BareSoil_z0,1)==sRows && size(BareSoil_z0,2)==sCols
        BareSoil_z0 = BareSoil_z0(r1:r2,c1:c2);
    end
    
    if size(CanopyHeight,1)==sRows && size(CanopyHeight,2)==sCols
        CanopyHeight = CanopyHeight(r1:r2,c1:c2);
    end
    
    if size(CanopyLeafAreaIndex,1)==sRows && size(CanopyLeafAreaIndex,2)==sCols
        CanopyLeafAreaIndex = CanopyLeafAreaIndex(r1:r2,c1:c2);
    end
    
    if size(CanopyDirectThroughfall,1)==sRows && size(CanopyDirectThroughfall,2)==sCols
        CanopyDirectThroughfall = CanopyDirectThroughfall(r1:r2,c1:c2);
    end
    
    
end

%%% convert land use to SLF format?

if OPTS.LU_TYPE==1
    % NLCD: https://www.mrlc.gov/nlcd11_leg.php
    % SLF: https://models.slf.ch/docserver/alpine3d/html/inputs.html
    
    %%% SLF    NAME    NLCD CODE(S)
    % 01	water	11
    % 02	settlement	21	22	23	24
    % 03	coniferous forest	42
    % 04	decidous forest	41
    % 05	mixed forest	43
    % 06	cereals
    % 07	pasture	81
    % 08	bush	52
    % 09	undefined
    % 10	undefined
    % 11	road
    % 12	undefined
    % 13	firn
    % 14	bare ice	12
    % 15	rock	31
    % 16	undefined
    % 17	undefined
    % 18	fruit
    % 19	vegetables
    % 20	wheat
    % 21	alpine vegetation
    % 22	wetlands	90	95
    % 23	rough pasture	51
    % 24	subalpine meadow	71
    % 25	alpine meadow	72 73 74
    % 26	bare soil vegetation
    % 27	free
    % 28	corn	82
    % 29	grapes
    % 30-99	undefined
    
    LU = SGRIDS.LU.*0;
    LU(SGRIDS.LU == 11) = 1;
    LU(SGRIDS.LU == 21 | SGRIDS.LU == 22 | SGRIDS.LU == 23 | SGRIDS.LU == 24) = 2;
    LU(SGRIDS.LU == 42) = 3;
    LU(SGRIDS.LU == 41) = 4;
    LU(SGRIDS.LU == 43) = 5;
    LU(SGRIDS.LU == 81) = 7;
    LU(SGRIDS.LU == 52) = 8;
    LU(SGRIDS.LU == 12) = 14;
    LU(SGRIDS.LU == 31) = 15;
    LU(SGRIDS.LU == 90 | SGRIDS.LU == 95) = 22;
    LU(SGRIDS.LU == 51) = 23;
    LU(SGRIDS.LU == 71) = 24;
    LU(SGRIDS.LU == 72 | SGRIDS.LU == 73 | SGRIDS.LU == 74) = 25;
    LU(SGRIDS.LU == 82) = 28;
elseif OPTS.LU_TYPE==2
    % LandFire: https://www.landfire.gov/DataDictionary/evc.pdf
    % SLF: https://models.slf.ch/docserver/alpine3d/html/inputs.html
    
    
    LU = landCover_LF2SLF(SGRIDS.LU);
    %     %%% initialize
    %     LU = SGRIDS.LU.*0;
    %
    %     %%% load table to convert land fire to SLF
    %     [num,~,~] = xlsread('land_use_LF2SLF.xls');
    %
    %     for nLF = 1:size(num,1)
    %         LU(SGRIDS.LU == num(nLF,1)) = num(nLF,2);
    %     end
else
    LU = SGRIDS.LU; % save a copy
end

if OPTS.LU_TYPE>=1
    %%% find any cells with class =0 and replace with mode in the area
    a = find(LU==0 | LU>=30);
    if isempty(a)==0
        LU(a) = NaN;
        di = 10;     % +/-  cells in each direction
        for j=1:numel(a)
            [r,c]=ind2ij(LU,a(j));
            r_min = nanmax([1 r-di]);
            r_max = nanmin([size(LU,1) r+di]);
            c_min = nanmax([1 c-di]);
            c_max = nanmin([size(LU,2) c+di]);
            LU(r,c) = mode(reshape(LU(r_min:r_max,c_min:c_max),numel(LU(r_min:r_max,c_min:c_max)),1));
        end
        a = find(LU==0 | LU>=30);
        if isempty(a)==0
            error('persistent empty land classes detected!')
        end
    end
end

%%% check to make sure we have land use in 1LLDC format
% SLF land use code = 1LLDC, where LL=land use code, D is soil depth
% (between 1 and 9) and C is field capacity (between 1 and 9)
% FOR NOW: set D = 03 and C = 01
SGRIDS.LU = LU;
a = find(SGRIDS.LU < 10100);
SGRIDS.LU(a) = SGRIDS.LU(a).*0 + 10000 + (LU(a).*100) + 30 + 1;
% SGRIDS.LU = uint16(SGRIDS.LU);
a = find(SGRIDS.LU > 19999);
if isempty(a)==0
    error('invalid LU code detected')
end
a = find(SGRIDS.LU < 10100 & SGRIDS.LU > 0);
if isempty(a)==0
    error('invalid LU code detected')
end

%%% compute slope and aspect from DEM
% if OPTS.SNOFILE_OPTION==1
    disp('... computing slope and aspect from DEM')
    SGRIDS.LAT = SGRIDS.XX.*NaN;
    SGRIDS.LON = SGRIDS.XX.*NaN;
    for j=1:size(SGRIDS.XX,1)
        for k=1:size(SGRIDS.XX,2)
            [SGRIDS.LAT(j,k), SGRIDS.LON(j,k)] = utm2deg(SGRIDS.XX(j,k), SGRIDS.YY(j,k), SGRIDS.utmzone);
        end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [SGRIDS.ASP, SGRIDS.SLP,~,~] = gradientm(SGRIDS.LAT, SGRIDS.LON, SGRIDS.ZZ);
%     SGRIDS.ASP = SGRIDS.ZZ .* 0 + 104;
%     SGRIDS.SLP = SGRIDS.ZZ .* 0 + 5;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end

% SGRIDS.ASP(SGRIDS.SLP == 0) = 180; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

demfile = [OPTS.EXP_NAME '.dem'];
lusfile = [OPTS.EXP_NAME '.lus'];
slpfile = [OPTS.EXP_NAME '.slp'];
aspfile = [OPTS.EXP_NAME '.asp'];
if OPTS.ALS_SCALING == 1
    alsfile = [OPTS.EXP_NAME '_als_' A3D_SETUP.ALS_DATE '.txt'];
end

if OPTS.FILES_WRITE(3)==1
    disp('... writing surface-grids')
    
    %%% write dem file
    dem_pathfile = cellstr([PATH_GRID demfile]);
    ascii_raster_writer(SGRIDS.XX, SGRIDS.YY, SGRIDS.ZZ, SGRIDS.grid_size, nodata, dem_pathfile);
    
    %%% write land use file
    lus_pathfile = cellstr([PATH_GRID lusfile]);
    ascii_raster_writer(SGRIDS.XX, SGRIDS.YY, SGRIDS.LU, SGRIDS.grid_size, nodata, lus_pathfile);
    
    %%% write slope file (not used as model input, just for analysis)
    slp_pathfile = cellstr([PATH_GRID slpfile]);
    ascii_raster_writer(SGRIDS.XX, SGRIDS.YY, SGRIDS.SLP, SGRIDS.grid_size, nodata, slp_pathfile);
    
    %%% write aspect file (not used as model input, just for analysis)
    asp_pathfile = cellstr([PATH_GRID aspfile]);
    ascii_raster_writer(SGRIDS.XX, SGRIDS.YY, SGRIDS.ASP, SGRIDS.grid_size, nodata, asp_pathfile);
    
    if OPTS.ALS_SCALING == 1
        %%% write als scaling file
        als_pathfile = cellstr([PATH_GRID alsfile]);
        ascii_raster_writer(SGRIDS.XX, SGRIDS.YY, SGRIDS.SD, SGRIDS.grid_size, nodata, als_pathfile);
    end
    
end



%% POI

disp('... writing POI files')

%%% write poi file
poifile = [OPTS.EXP_NAME '.poi'];

fid = fopen([PATH_GRID poifile], 'w');
fprintf(fid, '%s\n', 'SMET 1.1 ASCII');
fprintf(fid, '%s\n', '[HEADER]');
fprintf(fid, '%s\n', ['station_id' tb '=' tb 'my_pts']);
fprintf(fid, '%s\n', ['comment' tb '=' tb 'POI for full stratigraphy outputs']);
fprintf(fid, '%s\n', ['epsg' tb '=' tb num2str(epsg)]);
fprintf(fid, '%s\n', ['nodata' tb '=' tb '-999']);
fprintf(fid, '%s\n', ['fields' tb '=' tb 'easting northing altitude']);

if isfield(SGRIDS, 'POI')==1
    if isempty(SGRIDS.POI)==0
        
        nPOI = size(SGRIDS.POI,1);
        allPOIflag = 0;
        if numel(SGRIDS.POI)==1
            if SGRIDS.POI==1
                % then user wants all points as POI
                nPOI = numel(SGRIDS.XX);
                allPOIflag = 1;
            end
        end
        
        fprintf(fid, '%s\n', '[DATA]');
        for j=1:nPOI
            if allPOIflag==0
                %%% snap to grid
                [poiR, poiC] = nearest2D_mat(SGRIDS.POI(j,1), SGRIDS.POI(j,2), SGRIDS.XX, SGRIDS.YY, 1);
                save dd.mat
                poiX = SGRIDS.XX(poiR,poiC);
                poiY = SGRIDS.YY(poiR,poiC);
                poiZ = SGRIDS.ZZ(poiR,poiC);
            else
                poiX = SGRIDS.XX(j);
                poiY = SGRIDS.YY(j);
                poiZ = SGRIDS.ZZ(j);
            end
            
            if j<nPOI
                fprintf(fid, '%s\n', [num2str(poiX) ' ' num2str(poiY) ' ' num2str(poiZ)]);
            else
                fprintf(fid, '%s', [num2str(poiX) ' ' num2str(poiY) ' ' num2str(poiZ)]);
            end
        end
    else
        fprintf(fid, '%s', '[DATA]');
    end
    
else
    fprintf(fid, '%s', '[DATA]');
end



fclose(fid);   



%% misc

%%% get center of domain (x/y coords and lat/lon)
cenX = nanmean(SGRIDS.XX(:));
cenY = nanmean(SGRIDS.YY(:));

[cenLat, cenLon] = utm2deg(cenX, cenY, SGRIDS.utmzone);

%%% get mean elevation of domain
cenZ = nanmean(SGRIDS.ZZ(:));

%% write snowfiles

if OPTS.FILES_WRITE(2)==1
    %%%  "defaults" for snowfile params
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  1 2    3    4    5    6    7     8   9    10  11  12  13  14 15 16  17  18  19  20   21    22  23   24    25   26    %%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    def_CanopyHeight =            [0,NaN, 10,    6,   8, NaN, 0.05, 0.8,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN, 0.025, 0.6,NaN, 0.1, 0.1, 0.01];
    def_CanopyLeafAreaIndex =     [0,NaN,  2,    1.5,   1.75, NaN,    1,   4,NaN,NaN,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,     1,   1,NaN,   1,   1,  0.5];
    def_CanopyDirectThroughfall = [1,NaN, 0.4, 0.6, 0.5, NaN,  0.8, 0.1,NaN,NaN,NaN,NaN,NaN,NaN,1,NaN,NaN,NaN,NaN,NaN,   0.8, 0.5,NaN, 0.5, 0.5,  0.9];

    %%% initialize output matrices (which will be written to ascii grids)
    SGRIDS.CanopyHeight = SGRIDS.XX.*NaN;
    SGRIDS.CanopyLeafAreaIndex = SGRIDS.XX.*NaN;
    SGRIDS.CanopyDirectThroughfall = SGRIDS.XX.*NaN;
    
    disp('... writing snowfiles')
    
    if OPTS.SNOFILE_OPTION==0
        listLC = unique(SGRIDS.LU);
    elseif OPTS.SNOFILE_OPTION==1
        listLC = SGRIDS.LU;
    end
    
    disp('found these land use classes: ')
    tempLC=(listLC(~isnan(listLC)));
    uLLC = unique(tempLC);
    for j=1:numel(uLLC)
        disp(num2str(uLLC(j))) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    for j=1:numel(listLC)
        %%% current LC (5 digit)
        LCx = listLC(j);
        LUS = floor((LCx-10000)/100);
        
        %%% filename
        if OPTS.SNOFILE_OPTION==0
            snofile = [PATH_SNOW OPTS.EXP_NAME '_' num2str(LCx) '.sno'];
        elseif OPTS.SNOFILE_OPTION==1
            %%% convention is col_row where indexing is in C convention
            %%% (starts with 0)
            %%% NOTE row counts from bottom left up... so need to flip
            [si, sj] = ind2sub(size(SGRIDS.LU),j);  % si/sj only used for naming
            [gi, gj] = ind2sub(size(SGRIDS.LU),j);  % gi/gj used for indexing of geo variables
            si = size(SGRIDS.LU,1)-si+1;
            snofile = [PATH_SNOW num2str(sj-1) '_' num2str(si-1) '_' OPTS.EXP_NAME '.sno'];
        end
        
        %%% SITE_META values
        if OPTS.SNOFILE_OPTION==0
            SITE_META.STA_ID = num2str(LCx);
            SITE_META.STA_NAME = [OPTS.EXP_NAME '_' num2str(LCx)];
            SITE_META.STA_LON = cenLon;
            SITE_META.STA_LAT = cenLat;
            SITE_META.STA_ELEV = cenZ;
        elseif OPTS.SNOFILE_OPTION==1
            %%% grid specific meta
            SITE_META.STA_ID = [num2str(sj-1) '_' num2str(si-1)];
            SITE_META.STA_NAME = [num2str(sj-1) '_' num2str(si-1) '_' OPTS.EXP_NAME];
            SITE_META.STA_LAT = SGRIDS.LAT(gi,gj);
            SITE_META.STA_LON = SGRIDS.LON(gi,gj);
%             [SITE_META.STA_LAT,SITE_META.STA_LON]=utm2deg(SGRIDS.XX(si,sj), SGRIDS.YY(si,sj), SGRIDS.utmzone);
            SITE_META.STA_ELEV = SGRIDS.ZZ(gi,gj);
        end
        
        SITE_META.STA_NODATA = nodata;
        SITE_META.STA_TIMEZONE = SGRIDS.timezone;
        SITE_META.PROFILE_SDATE = SD_i;
        
        %%% slope and aspect
        if OPTS.SNOFILE_OPTION==0
            %%% dummy values
            SITE_META.SlopeAngle = 0;
            SITE_META.SlopeAzi = 0;
        elseif OPTS.SNOFILE_OPTION==1
            SITE_META.SlopeAngle = SGRIDS.SLP(gi,gj);
            SITE_META.SlopeAzi = SGRIDS.ASP(gi,gj);
        end
        
        %%% soil
        if numel(SoilAlbedo)==1
            SITE_META.SoilAlbedo = SoilAlbedo;
        elseif numel(SoilAlbedo)==numel(SGRIDS.LU)
            if OPTS.SNOFILE_OPTION==0
                %%% take average of all cells with this lus code
                alus = find(SGRIDS.LU==LCx);
                SITE_META.SoilAlbedo = nanmean(SoilAlbedo(alus));
            elseif OPTS.SNOFILE_OPTION==1
                %%% just use value from this grid
                SITE_META.SoilAlbedo = SoilAlbedo(gi,gj);
            end
        else
            error('invalid SoilAlbedo input size')
        end
        
        if numel(BareSoil_z0)==1
            SITE_META.BareSoil_z0 = BareSoil_z0;
        elseif numel(BareSoil_z0)==numel(SGRIDS.LU)
            if OPTS.SNOFILE_OPTION==0
                %%% take average of all cells with this lus code
                alus = find(SGRIDS.LU==LCx);
                SITE_META.BareSoil_z0 = nanmean(BareSoil_z0(alus));
            elseif OPTS.SNOFILE_OPTION==1
                %%% just use value from this grid
                SITE_META.BareSoil_z0 = BareSoil_z0(gi,gj);
            end
        else
            error('invalid BareSoil_z0 input size')
        end
        
        
        %%% veg
        if isempty(CanopyHeight)==1
            % then use default value
            SITE_META.CanopyHeight = def_CanopyHeight(LUS);
        else
            % then use input value
            if OPTS.SNOFILE_OPTION==0
                %%% take average of all cells with this lus code
                alus = find(SGRIDS.LU==LCx);
                SITE_META.CanopyHeight = nanmean(CanopyHeight(alus));
            elseif OPTS.SNOFILE_OPTION==1
                %%% just use value from this grid
                SITE_META.CanopyHeight = CanopyHeight(gi,gj);
            end
        end
        
        if OPTS.SNOFILE_OPTION==0
            currLC = find(SGRIDS.LU==LCx);
            SGRIDS.CanopyHeight(currLC) = SITE_META.CanopyHeight;
        elseif OPTS.SNOFILE_OPTION==1
            SGRIDS.CanopyHeight(gi,gj) = SITE_META.CanopyHeight;
        end
        
        
        
        
        if isempty(CanopyLeafAreaIndex)==1
            % then use default value
            SITE_META.CanopyLeafAreaIndex = def_CanopyLeafAreaIndex(LUS);
        else
            % then use input value
            if OPTS.SNOFILE_OPTION==0
                %%% take average of all cells with this lus code
                alus = find(SGRIDS.LU==LCx);
                SITE_META.CanopyLeafAreaIndex = nanmean(CanopyLeafAreaIndex(alus));
            elseif OPTS.SNOFILE_OPTION==1
                %%% just use value from this grid
                SITE_META.CanopyLeafAreaIndex = CanopyLeafAreaIndex(gi,gj);
            end
        end
        
        if OPTS.SNOFILE_OPTION==0
            currLC = find(SGRIDS.LU==LCx);
            SGRIDS.CanopyLeafAreaIndex(currLC) = SITE_META.CanopyLeafAreaIndex;
        elseif OPTS.SNOFILE_OPTION==1
            SGRIDS.CanopyLeafAreaIndex(gi,gj) = SITE_META.CanopyLeafAreaIndex;
        end
        
        
        
        if isempty(CanopyDirectThroughfall)==1
            % then use default value
            SITE_META.CanopyDirectThroughfall = def_CanopyDirectThroughfall(LUS);
        else
            % then use input value
            if OPTS.SNOFILE_OPTION==0
                %%% take average of all cells with this lus code
                alus = find(SGRIDS.LU==LCx);
                SITE_META.CanopyDirectThroughfall = nanmean(CanopyDirectThroughfall(alus));
            elseif OPTS.SNOFILE_OPTION==1
                %%% just use value from this grid
                SITE_META.CanopyDirectThroughfall = CanopyDirectThroughfall(gi,gj);
            end
        end
        
        if OPTS.SNOFILE_OPTION==0
            currLC = find(SGRIDS.LU==LCx);
            SGRIDS.CanopyDirectThroughfall(currLC) = SITE_META.CanopyDirectThroughfall;
        elseif OPTS.SNOFILE_OPTION==1
            SGRIDS.CanopyDirectThroughfall(gi,gj) = SITE_META.CanopyDirectThroughfall;
        end
        
        
        %%% write the file
        write_slf_snowfile(snofile, SITE_META);
        clear SITE_META
    end
    
    %%% compute fveg from SGRIDS.CanopyDirectThroughfall
    SGRIDS.fVEG = 1-SGRIDS.CanopyDirectThroughfall;
    
    %%% write ascii grids for veg data
    CanHeightfile = [OPTS.EXP_NAME '.hcan'];
    LAIfile = [OPTS.EXP_NAME '.lai'];
    fVEGfile = [OPTS.EXP_NAME '.fveg'];

    CanopyHeight_pathfile = cellstr([PATH_GRID CanHeightfile]);
    LAI_pathfile = cellstr([PATH_GRID LAIfile]);
    fVEG_pathfile = cellstr([PATH_GRID fVEGfile]);
    
    ascii_raster_writer(SGRIDS.XX, SGRIDS.YY, SGRIDS.CanopyHeight, SGRIDS.grid_size, nodata, CanopyHeight_pathfile);
    ascii_raster_writer(SGRIDS.XX, SGRIDS.YY, SGRIDS.CanopyLeafAreaIndex, SGRIDS.grid_size, nodata, LAI_pathfile);
    ascii_raster_writer(SGRIDS.XX, SGRIDS.YY, SGRIDS.fVEG, SGRIDS.grid_size, nodata, fVEG_pathfile);
    
    
end
%% write meteo files

nsta = size(METEO.STA_ID,2);  
if OPTS.FILES_WRITE(1)==1
    disp('... writing meteo files')
    
      
    
    for j=1:nsta
        %%% make sure no spaces in METEO.STA_NAME (replace w/ _)
        sta_nameX = char(METEO.STA_NAME(j));
        sta_nameX(isspace(sta_nameX)) = '_';
        METEO.STA_NAME(j) = {sta_nameX};
        
        %%% filename
        metfile = [PATH_MET char(METEO.STA_ID(j)) '.smet'];
        
        %%% meta structure
        SITE_META.STA_ID = METEO.STA_ID(j);
        SITE_META.STA_NAME = METEO.STA_NAME(j);
        SITE_META.STA_LAT = METEO.STA_LAT(j);
        SITE_META.STA_LON = METEO.STA_LON(j);
        SITE_META.STA_ELEV = METEO.STA_ELEV(j);
        
        %%% calc easting/northing based on lat/lon
        [SITE_META.STA_EASTING,SITE_META.STA_NORTHING,~] = deg2utm(SITE_META.STA_LAT, SITE_META.STA_LON);
        
        %%% domain values (assumed to not change from sta to sta)
        SITE_META.STA_EPSG = epsg;
        SITE_META.STA_NODATA = nodata;
        SITE_META.STA_TIMEZONE = SGRIDS.timezone;
        
        %%% met structure
        SITE_METEO.TIME = METEO.TIME;
        
        % air temp (K)
        if isfield(METEO, 'Tair')==1
            if isnan(nanmean(METEO.Tair(:,j)))==0
                if nanmean(METEO.Tair(:,j)<50)
                    SITE_METEO.TA = METEO.Tair(:,j)+273.15;
                else
                    SITE_METEO.TA = METEO.Tair(:,j);
                end
            end
        end
        
        % precip
        if isfield(METEO, 'P')==1
            if isnan(nanmean(METEO.P(:,j)))==0
                SITE_METEO.PSUM = METEO.P(:,j);
            end
        end
        
        % wind speed (m/s)
        if isfield(METEO, 'Us')==1
            if isnan(nanmean(METEO.Us(:,j)))==0
                SITE_METEO.VW = METEO.Us(:,j);
            end
        end
        
        % wind direction
        if isfield(METEO, 'Ud')==1
            if isnan(nanmean(METEO.Ud(:,j)))==0
                SITE_METEO.DW = METEO.Ud(:,j);
            end
        end
        
        % RH (fractional)
        if isfield(METEO, 'RH')==1
            if isnan(nanmean(METEO.RH(:,j)))==0
                if nanmean(METEO.RH(:,j)>1)
                    SITE_METEO.RH = METEO.RH(:,j)./100;
                else
                    SITE_METEO.RH = METEO.RH(:,j);
                end
            end
        end
        
        % incoming SW (W/m2)
        if isfield(METEO, 'Qsi')==1
            if isnan(nanmean(METEO.Qsi(:,j)))==0
                SITE_METEO.ISWR = METEO.Qsi(:,j);
            end
        end
        
        % incoming LW (W/m2)
        if isfield(METEO, 'Qli')==1
            if isnan(nanmean(METEO.Qli(:,j)))==0
                SITE_METEO.ILWR = METEO.Qli(:,j);
            end
        end
        
        
        %%% write the file
        write_slf_meteo(metfile, SITE_META, SITE_METEO);
        
        clear SITE_META SITE_METEO
    end
end

%% write io.ini file

if OPTS.FILES_WRITE(4)==1
    %%% need to refine this. for now, just copying my "standard" io.ini and
    %%% allowing some changes through matlab
    
    
    
    %%% start file
    pathfile_io = [OPTS.A3D_ROOT 'io.ini'];
    fid = fopen(pathfile_io, 'w');
    
    %%% GENERAL
    fprintf(fid, '%s\n', '[GENERAL]');
    fprintf(fid, '%s\n', ['BUFF_CHUNK_SIZE' tb '=' tb '370']);
    fprintf(fid, '%s\n', ['BUFF_BEFORE' tb '=' tb '1.5']);
    fprintf(fid, '%s\n', '');
    
    %%% INPUT
    fprintf(fid, '%s\n', '[INPUT]');
    fprintf(fid, '%s\n', ['COORDSYS' tb '=' tb 'UTM']);
    fprintf(fid, '%s\n', ['COORDPARAM' tb '=' tb num2str(utmzoneCompact)]);
    fprintf(fid, '%s\n', ['TIME_ZONE' tb '=' tb num2str(SGRIDS.timezone)]);

    fprintf(fid, '%s\n', ['METEO' tb '=' tb 'SMET']);
    fprintf(fid, '%s\n', ['METEOPATH' tb '=' tb '../input/meteo']);
    
    for j=1:nsta
        fprintf(fid, '%s\n', ['STATION' num2str(j) tb '=' tb char(METEO.STA_ID(j))]);
    end

    fprintf(fid, '%s\n', ['ISWR_IS_NET' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['SNOWPATH' tb '=' tb '../input/snowfiles']);
    fprintf(fid, '%s\n', ['SNOW' tb '=' tb 'SMET']);
    
    fprintf(fid, '%s\n', ['GRID2D' tb '=' tb 'ARC']);
    fprintf(fid, '%s\n', ['GRID2DPATH' tb '=' tb '../input/surface-grids']);
    
    fprintf(fid, '%s\n', ['DEM' tb '=' tb 'ARC']);
    fprintf(fid, '%s\n', ['DEMFILE' tb '=' tb ' ../input/surface-grids/' char(demfile)]);
  
    fprintf(fid, '%s\n', ['LANDUSE' tb '=' tb 'ARC']);
    fprintf(fid, '%s\n', ['LANDUSEFILE' tb '=' tb ' ../input/surface-grids/' char(lusfile)]);
  
    fprintf(fid, '%s\n', ['POI' tb '=' tb 'SMET']);
    fprintf(fid, '%s\n', ['POIFILE' tb '=' tb ' ../input/surface-grids/' char(poifile)]);
    fprintf(fid, '%s\n', '');
 

    %%% OUTPUT
    fprintf(fid, '%s\n', '[OUTPUT]');
    fprintf(fid, '%s\n', ['COORDSYS' tb '=' tb 'UTM']);
    fprintf(fid, '%s\n', ['COORDPARAM' tb '=' tb num2str(utmzoneCompact)]);
    fprintf(fid, '%s\n', ['TIME_ZONE' tb '=' tb num2str(SGRIDS.timezone)]);

    fprintf(fid, '%s\n', ['EXPERIMENT' tb '=' tb char(OPTS.EXP_NAME)]);
    
    fprintf(fid, '%s\n', ['METEO' tb '=' tb 'SMET']);
    fprintf(fid, '%s\n', ['METEOPATH' tb '=' tb '../output']);
    
    fprintf(fid, '%s\n', ['SNOW_WRITE' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['SNOW_DAYS_BETWEEN' tb '=' tb '1']);
    fprintf(fid, '%s\n', ['SNOW' tb '=' tb 'SMET']);
    fprintf(fid, '%s\n', ['SNOWPATH' tb '=' tb '../output/snowfiles']);
    fprintf(fid, '%s\n', ['BACKUP_DAYS_BETWEEN' tb '=' tb '365.0']);
    fprintf(fid, '%s\n', ['FIRST_BACKUP' tb '=' tb '400.0']);
    
    if flag_gridParams==0
        fprintf(fid, '%s\n', ['GRIDS_WRITE' tb '=' tb 'FALSE']);
    else
        fprintf(fid, '%s\n', ['GRIDS_WRITE' tb '=' tb 'TRUE']);
        fprintf(fid, '%s\n', ['GRIDS_DAYS_BETWEEN' tb '=' tb num2str(A3D_SETUP.OUTPUTS.GRIDS_DAYS_BETWEEN)]);
        fprintf(fid, '%s\n', ['GRIDS_START' tb '=' tb num2str(A3D_SETUP.OUTPUTS.GRIDS_START)]);
        
        fprintf(fid, '%s', 'GRIDS_PARAMETERS =');
        for j=1:numel(A3D_SETUP.OUTPUTS.GRIDS_PARAMETERS)
            if j<numel(A3D_SETUP.OUTPUTS.GRIDS_PARAMETERS)
                fprintf(fid, '%s', [' ' char(A3D_SETUP.OUTPUTS.GRIDS_PARAMETERS(j))]);
            else
                fprintf(fid, '%s\n', [' ' char(A3D_SETUP.OUTPUTS.GRIDS_PARAMETERS(j))]);
            end
        end
        
        fprintf(fid, '%s\n', ['SOIL_TEMPERATURE_DEPTH' tb '=' tb '1']);
        
        fprintf(fid, '%s\n', ['GRID2D' tb '=' tb 'ARC']);
        fprintf(fid, '%s\n', ['GRID2DPATH' tb '=' tb '../output/grids']);
        fprintf(fid, '%s\n', ['A3D_VIEW' tb '=' tb 'TRUE']);
    end
    
    fprintf(fid, '%s\n', ['PROF_WRITE' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['PROF_FORMAT' tb '=' tb 'PRO']);
    fprintf(fid, '%s\n', ['PROF_START' tb '=' tb '0.0']);
    fprintf(fid, '%s\n', ['PROF_DAYS_BETWEEN' tb '=' tb '0.083333']);
    fprintf(fid, '%s\n', ['HARDNESS_IN_NEWTON' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['CLASSIFY_PROFILE' tb '=' tb 'FALSE']);

    fprintf(fid, '%s\n', ['TS_WRITE' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['TS_START' tb '=' tb '0.0']);
    fprintf(fid, '%s\n', ['TS_DAYS_BETWEEN' tb '=' tb '0.083333']);
    fprintf(fid, '%s\n', ['AVGSUM_TIME_SERIES' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['CUMSUM_MASS' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['PRECIP_RATES' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['OUT_CANOPY' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['OUT_HAZ' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['OUT_HEAT' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['OUT_T' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['OUT_LW' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['OUT_SW' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['OUT_MASS' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['OUT_METEO' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['OUT_STAB' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', '');

    %%% SNOWPACK
    fprintf(fid, '%s\n', '[SNOWPACK]');
    fprintf(fid, '%s\n', ['CALCULATION_STEP_LENGTH' tb '=' tb '15']);
    fprintf(fid, '%s\n', ['ROUGHNESS_LENGTH' tb '=' tb num2str(A3D_SETUP.SNOWPACK.ROUGHNESS_LENGTH)]);
    fprintf(fid, '%s\n', ['HEIGHT_OF_METEO_VALUES' tb '=' tb num2str(A3D_SETUP.SNOWPACK.HEIGHT_OF_METEO_VALUES)]);
    fprintf(fid, '%s\n', ['HEIGHT_OF_WIND_VALUE' tb '=' tb num2str(A3D_SETUP.SNOWPACK.HEIGHT_OF_WIND_VALUE)]);
    fprintf(fid, '%s\n', ['ENFORCE_MEASURED_SNOW_HEIGHTS' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['SW_MODE' tb '=' tb 'INCOMING']);
    fprintf(fid, '%s\n', ['ATMOSPHERIC_STABILITY' tb '=' tb char(A3D_SETUP.SNOWPACK.ATMOSPHERIC_STABILITY)]);
    fprintf(fid, '%s\n', ['CANOPY' tb '=' tb char(A3D_SETUP.SNOWPACK.CANOPY)]);
    fprintf(fid, '%s\n', ['MEAS_TSS' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['CHANGE_BC' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['THRESH_CHANGE_BC' tb '=' tb '-1.3']);
    fprintf(fid, '%s\n', ['SNP_SOIL' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['SOIL_FLUX' tb '=' tb 'FALSE']);
    fprintf(fid, '%s\n', ['GEO_HEAT' tb '=' tb '0.06']);
    fprintf(fid, '%s\n', '');

    %%% EBALANCE
    fprintf(fid, '%s\n', '[EBALANCE]');
    
    fprintf(fid, '%s\n', ['TERRAIN_RADIATION' tb '=' tb 'TRUE']);
    fprintf(fid, '%s\n', ['TERRAIN_RADIATION_METHOD' tb '=' tb 'SIMPLE']);
    fprintf(fid, '%s\n', '');
 
    %%% FILTERS
    fprintf(fid, '%s\n', '[FILTERS]');
    
    fprintf(fid, '%s\n', ['TA::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['TA::arg1' tb '=' tb '240 320']);

    fprintf(fid, '%s\n', ['RH::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['RH::arg1' tb '=' tb '0.01 1.2']);
    fprintf(fid, '%s\n', ['RH::filter2' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['RH::arg2' tb '=' tb 'soft 0.05 1.0']);

    fprintf(fid, '%s\n', ['PSUM::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['PSUM::arg1' tb '=' tb '-0.1 100.']);
    fprintf(fid, '%s\n', ['PSUM::filter2' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['PSUM::arg2' tb '=' tb 'soft 0.0 100.0']);
    fprintf(fid, '%s\n', ['PSUM::filter3' tb '=' tb 'undercatch_wmo']);
    fprintf(fid, '%s\n', ['PSUM::arg3' tb '=' tb 'Hellmannsh']);

    fprintf(fid, '%s\n', ['VW::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['VW::arg1' tb '=' tb '-2 70']);
    fprintf(fid, '%s\n', ['VW::filter2' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['VW::arg2' tb '=' tb 'soft 0.2 50.0']);

    fprintf(fid, '%s\n', ['DW::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['DW::arg1' tb '=' tb '0 360']);

    fprintf(fid, '%s\n', ['VW_MAX::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['VW_MAX::arg1' tb '=' tb '-2 70']);
    fprintf(fid, '%s\n', ['VW_MAX::filter2' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['VW_MAX::arg2' tb '=' tb 'soft 0.2 50.0']);

    fprintf(fid, '%s\n', ['ISWR::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['ISWR::arg1' tb '=' tb '-10 1500']);
    fprintf(fid, '%s\n', ['ISWR::filter2' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['ISWR::arg2' tb '=' tb 'soft 0 1500']);

    fprintf(fid, '%s\n', ['ILWR::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['ILWR::arg1' tb '=' tb '100 600']);
    fprintf(fid, '%s\n', ['ILWR::filter2' tb '=' tb 'max']);
    fprintf(fid, '%s\n', ['ILWR::arg2' tb '=' tb 'soft 500']);

    fprintf(fid, '%s\n', ['TSS::filter1' tb '=' tb 'min_max']);
    fprintf(fid, '%s\n', ['TSS::arg1' tb '=' tb '200 320']);
    fprintf(fid, '%s\n', '');
    
    %%% INTERPOLATIONS1D
    fprintf(fid, '%s\n', '[INTERPOLATIONS1D]');
    fprintf(fid, '%s\n', ['WINDOW_SIZE' tb '=' tb '86400']);

    fprintf(fid, '%s\n', ['PSUM::resample' tb '=' tb 'accumulate']);
    fprintf(fid, '%s\n', ['PSUM::accumulate' tb '=' tb '3600']);
    fprintf(fid, '%s\n', '');

    %%% INTERPOLATIONS2D
    fprintf(fid, '%s\n', '[Interpolations2D]');
    
    fprintf(fid, '%s\n', ['TA::algorithms' tb '=' tb 'ODKRIG_LAPSE IDW_LAPSE AVG_LAPSE']);
    fprintf(fid, '%s\n', ['TA::avg_lapse' tb '=' tb '-0.006']);
    fprintf(fid, '%s\n', ['TA::idw_lapse' tb '=' tb '-0.006 soft']);

    fprintf(fid, '%s\n', ['RH::algorithms' tb '=' tb 'LISTON_RH IDW_LAPSE']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OPTS.ALS_SCALING == 1
        fprintf(fid, '%s\n', ['PSUM::algorithms' tb '=' tb 'ALS_SCALING']);
        fprintf(fid, '%s\n', ['PSUM::als_scaling' tb '=' tb 'idw_lapse ' alsfile ' 275.']);
    elseif OPTS.ALS_SCALING == 0 
        fprintf(fid, '%s\n', ['PSUM::algorithms' tb '=' tb 'IDW_LAPSE AVG_LAPSE AVG CST']);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid, '%s\n', ['PSUM::idw_lapse' tb '=' tb '0.0005 frac']);
    fprintf(fid, '%s\n', ['PSUM::avg_lapse' tb '=' tb '0.0005 frac']);
    fprintf(fid, '%s\n', ['PSUM::cst' tb '=' tb '0']);

    fprintf(fid, '%s\n', ['PSUM_PH::algorithms' tb '=' tb 'PPHASE']);
    fprintf(fid, '%s\n', ['PSUM_PH::PPHASE' tb '=' tb 'THRESH 274.35']);
    
    fprintf(fid, '%s\n', ['VW::algorithms' tb '=' tb 'LISTON_WIND IDW_LAPSE AVG']);

    fprintf(fid, '%s\n', ['DW::algorithms' tb '=' tb 'LISTON_WIND']);

    fprintf(fid, '%s\n', ['VW_MAX::algorithms' tb '=' tb 'IDW_LAPSE']);

    fprintf(fid, '%s\n', ['P::algorithms' tb '=' tb 'STD_PRESS']);
    
    fprintf(fid, '%s\n', ['ILWR::algorithms' tb '=' tb 'ILWR_EPS IDW_LAPSE AVG_LAPSE IDW']);
    fprintf(fid, '%s\n', ['ILWR::ilwr_eps::rate' tb '=' tb '-1.8e-5']);
    
    fprintf(fid, '%s\n', ['ISWR::algorithms' tb '=' tb 'SWRAD']);
    fprintf(fid, '%s\n', '');

    
    %%% GENERATORS
    fprintf(fid, '%s\n', '[GENERATORS]');



    %%% close file
    fclose(fid);
    
    
    %%% the commented area below are variables I found when searching for
    %%% "cfg.get" in the *.cc within alpine3d source folder
    
%     %%% [GENERAL]
%     LOCAL_IO %%%
%     
%     %%% [INPUT]
%     COORDSYS
%     COORDPARAM
%     TIME_ZONE
%     
%     COMPUTE_IN_LOCAL_COORDS
%     CATCHMENT
%     
%     CATCHMENT_NUMBERING
%     
%     WINDFIELDS
%     
%     SNOW
%     
%     
%     
%     %%% [OUTPUT]
%     GRID2DPATH
%     TIME_ZONE
%     
%     EXPERIMENT
%     
%     GRIDS_WRITE
%     GRIDS_START
%     GRIDS_DAYS_BETWEEN
%     
%     TS_WRITE
%     TS_START
%     TS_DAYS_BETWEEN
%     
%     PROF_WRITE
%     PROF_START
%     PROF_DAYS_BETWEEN
%     
%     WRITE_RUNOFF_GRIDS
%     SNOW_WRITE
%     SNOW_DAYS_BETWEEN
%     
%     METEOPATH
%     
%     CATCHMENTS_PATH
%     RUNOFF_GRID2D
%     RUNOFF_GRID2DFILE
%     
%     MASK_DYNAMIC
%     MASK_GLACIERS
%     
%     RUNOFF_FILES_EXTRA_DATA
%     
%     GRIDS_PARAMETERS
%     
%     SOIL_TEMPERATURE_DEPTH
%     
%     %%% [ALPINE3D]
%     RUNOFF_CFG
%     
%     %%% [SNOWPACK]
%     CALCULATION_STEP_LENGTH
%     CANOPY
%     SNP_SOIL
%     
%     GLACIER_KATABATIC_FLOW
%     KATABATIC_LAYER_HEIGHT
%     KATABATIC_SCALING
%     KATABATIC_K_COEFFICIENT
%     
%     HEIGHT_OF_WIND_VALUE
%     
%     %%% [SNOWPACKADVANCED]
%     THRESH_RAIN
%     ADJUST_HEIGHT_OF_WIND_VALUE
%     ADJUST_HEIGHT_METEO_VALUE
%     
%     WATERTRANSPORTMODEL_SNOW
%     WATERTRANSPORTMODEL_SOIL
%     
%     %%% [EBALANCE]
%     Terrain_Radiation
%     Terrain_Radiation_Method
%     itEps_SW
%     itEps_LW
%     itEps1_SW
%     sw_radius
%     lw_radius
%     sub_crit
%     vf_in_ram
%     vf_file
%     tvfarea
%     
%     %%% [FILTERS]
%     
%     %%% [INTERPOLATIONS1D]
%     
%     
%     %%% [Interpolations2D]
%     TA::algorithms
%     
%     %%% [GENERATORS]
%     
%     
%     
%     
%     
%     
%     TIME_ZONE
%     CALCULATION_STEP_LENGTH
%     
%     
end

save dump.mat
