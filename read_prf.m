% Reads the 'prf' (state variables by layer) output from FSM2 (compiled Jan 2019)
%   Right now, only reads the snow layers, not soil
%
% RELEASE NOTES
%   Version 1.0 written by Eric Smyth (esmyth@esmyth.com) Mar 2020
%
% SYNTAX
%   PROFILE = read_prf(prf_file,vars,Nt,max_layers)
%
% INPUTS
%   prf_file = file name (path if necessary) of 'prf' file that you're trying to read
%   vars = cell array of strings, specifying at least one of the following
%           snow variables that you want to read.
%               MID_H = position of midpoint of snow layer above the ground
%               DS = thickness of snow layer
%               TEMP = snow layer temperature (K)
%               RADI = snow grain radii in layer
%               ICE = snow layer ice content
%               LIQ = snow layer liquid content
%   Nt = number of timesteps in FSM2 output (should be same as number of 
%               input data timesteps)
%   max_layers = maximum number of snow layers, specified in the FSM2
%               configuration before you run
% 
% OUTPUTS
%   PROFILE = structure, with fields corresponding to specified variables
%               Each field has Nt rows * max_layers columns
%               NaN values indicate that there was no snow in that given 
%                   layer/timestep


function PROFILE = read_prf(prf_file,vars,Nt,max_layers)

% prf_file = '999prf';
% vars = {'LIQ'};
% Nt = 8760;
% max_layers = 15;

%% Read prf file into matlab, clean up infinity strings
% fid = fopen(prf_file,'rt') ;
% X = fread(fid) ;
% fclose(fid) ;
% X = char(X.') ;
% Y = strrep(X, 'Infinity', '0.00000E+00');
% fid2 = fopen('888888prf','wt') ;
% fwrite(fid2,Y) ;
% fclose (fid2) ;
% MATRIX = importdata('888888prf');

find_and_replace(prf_file, 'Infinity', '0.00000E+00')
MATRIX = importdata(prf_file);

%% Set up profile structure
nvars = numel(vars);
for hh = 1:nvars
    varName = char(vars{hh});
    eval(['PROFILE. ' varName ' = nan(Nt,max_layers);'])
end

%% Variable library
var_lib = {'MID_H','DS','TEMP','RADI','ICE','LIQ'};

%% Work through rows of Matrix
prf_row_counter = 1;
for ii = 1:length(MATRIX(:,1))
    % First number indicates if this is a date row
    indicator = MATRIX(ii,1);
    if indicator > 1500
        % This is a date row, grab # of snow layers
        num_layers = MATRIX(ii,4);
        % Loop through number of layers, for each variable
        for cc = 1:nvars
            curr_var = char(vars{cc});
            curr_ind = find(strcmp(var_lib, curr_var));
            for gg = 1:num_layers
                if curr_ind <= length(MATRIX(1,:))
                    eval(['PROFILE.' curr_var '(prf_row_counter,gg) = MATRIX(ii+(gg*2)-1,curr_ind);'])
                else
                    eval(['PROFILE.' curr_var '(prf_row_counter,gg) = MATRIX(ii+(gg*2),curr_ind-length(MATRIX(1,:)));'])
                end
            end
        end
        prf_row_counter = prf_row_counter + 1;
    else
        % Not a date row, do nothing
    end
    
end


end