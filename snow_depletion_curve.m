function [fSCA] = snow_depletion_curve(arg1,arg2,arg3,arg4,arg5,arg6,varargin)

%RELEASE NOTES
%   Written by Eric Smyth Dec 2021
%   esmyth@esmyth.com
%
% ****SYNTAX****
%   OPTION 1: Simplest curve, based on SD and a defined slope. No physical
%   meaning, though it is close to other "real" curves and requires less data
%   [fSCA] = snow_depletion_curve(curve_choice,SD,fSCA_slope)
%
%   OPTION 2: Curve based on depth, density, and other params
%   Based on Niu & Yang, 2007
%   [fSCA] = snow_depletion_curve(curve_choice,SD,RHO,RHO_new,fSCA_z0,fSCA_m)
%
%   OPTION 3: Curve based on SWE and subgrid DEM variability
%   Based on Roesch et al., 2001
%   [fSCA] = snow_depletion_curve(curve_choice,SWE,stdDEM,fSCA_limit)
%
% ****INPUTS****
% curve_choice
%   Integer value, corresponding to the options listed above:
%       1 = OPTION 1
%       2 = OPTION 2
%       3 = OPTION 3
% SD
%   Matrix or array of snow depth values, to be converted into equivalent
%   fSCA values, using the defined curve (meters)
% fSCA_slope
%   If option 1, defines how quickly the curve approaches 1 (smaller
%   numbers approach 1 more quickly). Example: 0.1
% RHO
%   If option 2, matrix or array of snow density values, which must be the
%   same dimensions as SD (kg / m^-3)
% RHO_new
%   If option 2, a parameter for the density of freshly fallen snow (e.g.,
%   100 kg / m^-3)
% fSCA_z0
%   If option 2, a parameter for surface roughness (e.g., 0.02, unitless)
% fSCA_m
%   If option 2, a calibrated parameter, example: 1.6 (unitless)
% SWE
%   If option 3, a matrix or array of snow water equivalent values to be
%   converted into fSCA, using the defined curve (mm)
% stdDEM
%   If option 3, a matrix or array (same dimensions as SWE) with values for
%   the standard deviation of elevations within the given pixel
%   corresponding to each SWE value (meters)
% fSCA_limit
%   If option 3, the maximum value that fSCA approaches with high SWE
%   (e.g., 0.95 or 1)
%
% ****OUTPUTS****
% fSCA
%   Matrix or array (same dimensions as SD or SWE depending on options)
%   with fractional snow covered area values (unitless)

%% Parse inputs
curve_choice = arg1;

if curve_choice == 1
    SD = arg2;
    fSCA_slope = arg3;
    
elseif curve_choice == 2
    SD = arg2;
    RHO = arg3;
    RHO_new = arg4;
    fSCA_z0 = arg5;
    fSCA_m = arg6;
    
elseif curve_choice == 3
    SWE = arg2;
    stdDEM = arg3;
    fSCA_limit = arg4;
    
else
    error('Please choose a valid curve_choice')
end

%% Generate curve
if curve_choice == 1
    %%% OPTION 1
    
    fSCA = nan(size(SD));

    for ii = 1:length(SD(1,:))
        for gg = 1:length(SD(:,1))
            
            curr_SD = SD(gg,ii);
            fSCA(gg,ii) = curr_SD / (curr_SD + fSCA_slope);
            
        end
    end
    
elseif curve_choice == 2
    %%% OPTION 2
    
    fSCA = nan(size(SD));

    for ii = 1:length(SD(1,:))
        for gg = 1:length(SD(:,1))

            a = 2.5 * fSCA_z0 * ((RHO(gg,ii) / RHO_new)^fSCA_m);

            curr_SD = SD(gg,ii);
            fSCA(gg,ii) = tanh(curr_SD / a);

        end
    end
    
elseif curve_choice == 3
    %%% OPTION 3
    
    fSCA = nan(size(SWE));

    for ii = 1:length(SWE(1,:))
        for gg = 1:length(SWE(:,1))

            curr_SWE = SWE(gg,ii);
            curr_std = stdDEM(gg,ii);

            fSCA(gg,ii) = fSCA_limit * tanh(0.1 * curr_SWE) * sqrt(curr_SWE/(curr_SWE + 0.0001 + (0.15 * curr_std)));

        end
    end
    
end


end






