% Writes a variable to a pre-defined ncdf file
% 
% See 'mat2ncdf_wrapper' for an example matlab wrapper script, where you
% can simply define a file path containing any number of .mat files, and it
% will convert everything to .nc
%
% If the variable is a structure, this script will also create and write
% variables to the ncdf file, with the syntax:
%     variable.subvariable.sub2 -->  variable_subvariable_sub2
%
% All ncdf variable dimension names will be "[variable name]_rows" and
% "[variable name]_cols"
%
% NOTE: does not currently work with more than 2 dimensions
% NOTE: does not work with MATLAB versions before 2022
%
% RELEASE NOTES
%   Written by Eric Smyth (esmyth@esmyth.com), May 2022
%
% SYNTAX
%   mat2ncdf(varName,varData,ncfile)
%
% INPUTS
%   varName = string, name of variable to be written
%   varData = the variable itself, to be written
%   ncfile = name of ncdf file, already created, for the variable to be
%   written to
%
% OUTPUTS
%   Writes variable to ncdf file, no matlab outputs
%

function mat2ncdf(varName,varData,ncfile)

    if isempty(varData) == 1
        %%% If variable is empty, just skip it
    else

        if isstruct(varData) == 0
            %%% This is not a structure, just a normal variable
            varNameOut = strrep(varName,'.','_');
            varDataOut = varData;
            
            %%% Initialize variable in ncfile
            var_size = size(varDataOut);
            if iscell(varDataOut) == 1 || ischar(varDataOut) == 1 || isstring(varDataOut)
                %%% Current variable is a cell array, or character array,
                %%% or string - which needs to be written to ncdf with some
                %%% additional flags - specifically as a string (after
                %%% converting, if its a cell array, for example)
                [v,d] = version;
                d = str2num(d(end-3:end));
                if d < 2022
                    error('MATLAB releases before 2022 cant deal with writing strings to ncdf')
                end
                disp('...this variable is a cell array or character vector, converting to array of strings')
                nccreate(ncfile, char(varNameOut), 'Dimensions', {[char(varNameOut) '_rows'], ...
                    var_size(1), [char(varNameOut) '_cols'], var_size(2)}, ...
                    'Format','netcdf4','Datatype','string');
                varDataOut = string(varDataOut);

                %%% One common error stems from '<missing>' data, which
                %%% Matlab generates automatically sometimes when
                %%% converting to string arrays. Need to convert to simple
                %%% string saying 'NaN' in order to write to .nc
                if sum(ismissing(varDataOut)) > 0
                    missing_inds = ismissing(varDataOut);
                    varDataOut(missing_inds) = string('NaN');
                end

            else
                %%% Not a cell or character array - set up normally
                nccreate(ncfile, char(varNameOut), 'Dimensions', {[char(varNameOut) '_rows'], ...
                    var_size(1), [char(varNameOut) '_cols'], var_size(2)}, ...
                    'Format','netcdf4');

                %%% Check if variable is a logical, convert to double
                if islogical(varDataOut) == 1
                    varDataOut = double(varDataOut);
                end
            end
            
            %%% Write to ncdf
            ncwrite(ncfile, char(varNameOut), varDataOut);
    
        else
            %%% This variable is a structure, start recursion
    
            %%% Grab fieldnames in this layer
            FN_SUB = fieldnames(varData);
    
            %%% Loop through variables in this layer
            for gg = 1:numel(FN_SUB)
    
                %%% Save variable names and data
                varNameSave = varName;
                varDataSave = varData;
                
                %%% Rename variable
                varName = [varName '.' FN_SUB{gg}];
                disp(['Current variable: ' varName])
    
                %%% Grab variable data
                eval(['varData = varData.' FN_SUB{gg} ';'])
    
                %%% Call self
                mat2ncdf(varName,varData,ncfile);
                
                %%% Reset for next variable in the loop
                varName = varNameSave;
                varData = varDataSave;
    
            end  
        end
    end
end