%RELEASE NOTES
%   Written by Eric Smyth Dec 2020
%   esmyth@esmyth.com
%
%SYNTAX
%   [p1, p2, p3] = plot_filled_area(var_min, var_max, TIME, PARAMS)
%
% ****INPUTS****
% var_min 
%   An Lx1 array, the time series of data that is the bottom edge of the filled plot
%
% var_max
%   An Lx1 array, the time series of data that is the top edge of the filled plot
%
% TIME
%   An Lx1 array, the x values for the two time_series above
%
% PARAMS
%   A structure with the following fields
%       fillColor = an [R G B] triplet, the color of the filled area (default [240/255 248/255 1], a light blue]
%       edgeColor = an [R G B] triplet, the color of the bounding line (default black)
%       edgeLineWidth = a scalar, the width of the bounding line (default 1)
%       faceAlpha = a scalar between 0 (transparent) and 1 (opaque, default)
%       linestyle = a text string with the style of the bounding lines (default '-' solid line)
%
% ****OUTPUTS****
% p1, p2, and p3
%   The object handles of the top bounding line, the bottom bounding line,
%   and the fill area, respectively

function [p1, p2, p3] = plot_filled_area(var_min, var_max, TIME, PARAMS)

    %% Check everything is in right orientation
    if length(TIME(1,:)) > 1
        disp('Detecting TIME array is horizontal...changing to vertical for consistency')
        TIME = TIME';
    end
    if length(var_max(1,:)) > 1
        disp('Detecting var_max array is horizontal...changing to vertical for consistency')
        var_max = var_max';
    end
    if length(var_min(1,:)) > 1
        disp('Detecting var_min array is horizontal...changing to vertical for consistency')
        var_min = var_min';
    end
    
    %% Check for PARAMS and set defaults if necessary
    if isfield(PARAMS,'edgeColor') == 1
    else
        PARAMS.edgeColor = [0 0 0];
    end
    if isfield(PARAMS,'edgeLineWidth') == 1
    else
        PARAMS.edgeLineWidth = 1;
    end
    if isfield(PARAMS,'faceAlpha') == 1
    else
        PARAMS.faceAlpha = 1;
    end
    if isfield(PARAMS,'fillColor') == 1
    else
        PARAMS.fillColor = [240/255 248/255 1];
    end
    if isfield(PARAMS,'linestyle') == 1
    else
        PARAMS.linestyle = '-';
    end

    %% Plot outlines
    p1 = plot(TIME, var_max, PARAMS.linestyle, 'Color', PARAMS.edgeColor, ...
        'linewidth', PARAMS.edgeLineWidth);
    p2 = plot(TIME, var_min, PARAMS.linestyle, 'Color', PARAMS.edgeColor, ...
        'linewidth', PARAMS.edgeLineWidth);
    
    %% Generate fill area
    f1 = TIME;
    f2 = flipud(f1);
    f3 = var_min;
    f4 = flipud(var_max);
    f5 = [f1; f2];
    f6 = [f3; f4];
    
    %% Plot fill area
    p3 = fill(f5, f6, PARAMS.fillColor, 'FaceAlpha', PARAMS.faceAlpha, 'EdgeColor', 'none');

end