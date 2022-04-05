% Generates x and y min/max axis limits, given x and y coordinates and a
% user defined rounding value
%
% RELEASE NOTES
%   Written by Eric Smyth (esmyth@esmyth.com), Nov 2020
%
% SYNTAX
%    [x_min, x_max, y_min, y_max] = smart_axis_lims(Xround,Yround,X1,Y1)
%    [x_min, x_max, y_min, y_max] = smart_axis_lims(Xround,Yround,X1,Y1,X2,Y2)
%    [x_min, x_max, y_min, y_max] = smart_axis_lims(Xround,Yround,X1,Y1,X2,Y2,X3,Y3)
%
% INPUTS
%   Xround = a number N. The x_max value will be the lowest multiple of N that is 
%       higher than the highest number in X. E.g., if Xround is 25, and the
%       X array sits between 2 and 21, the x_max will be 25 and the x_min
%       will be 0.
%   Yround = a number N. Works the same way as Xround
%   X1 = array of numbers, what might be the x-coordinates on a 2D grid
%   Y1 = array of numbers, same length as X, that might be y-coordinates on a 2D grid
%   X2, Y2, etc. = other arrays of numbers, same lengths
%       The idea here is that sometimes you want your x and y axis limits
%       to be consistent between multiple plots. So, include the X and Y
%       data from these multiple plots. The code will find the max and min
%       of all the data, as a whole.
%
% OUTPUTS
%   x_min = a number that you can use as the lower bound of x axis
%   x_max = a number that you can use as the upper bound of x axis
%   y_min = a number that you can use as the lower bound of y axis
%   y_max = a number that you can use as the upper bound of y axis

function [x_min, x_max, y_min, y_max] = smart_axis_lims(arg01, arg02, arg03, arg04, arg05, arg06, arg07, arg08, varargin)

    %% Parse Inputs
    Xround = arg01;
    Yround = arg02;
    
    if nargin == 4
        % Just one X and Y array
        X_T = arg03;
        Y_T = arg04;
        
    elseif nargin == 6
        % Two X and Y arrays
        X_1 = arg03;
        Y_1 = arg04;
        X_2 = arg05;
        Y_2 = arg06;
        % Make sure they are in the same dimension
        try
            X_T = [X_1 X_2];
        catch
            X_T = [X_1' X_2];
        end
        try
            Y_T = [Y_1 Y_2];
        catch
            Y_T = [Y_1' Y_2];
        end
        
    elseif nargin == 8
        % Three X and Y arrays
        X_1 = arg03;
        Y_1 = arg04;
        X_2 = arg05;
        Y_2 = arg06;
        X_3 = arg07;
        Y_3 = arg08;
        % Make sure they are in the same dimension
        try
            X_T = [X_1 X_2];
        catch
            X_T = [X_1' X_2];
        end
        try
            Y_T = [Y_1 Y_2];
        catch
            Y_T = [Y_1' Y_2];
        end
        try
            X_T = [X_T X_3];
        catch
            X_T = [X_T X_3'];
        end
        try
            Y_T = [Y_T Y_3];
        catch
            Y_T = [Y_T X_3'];
        end
    end
    
    %% Find upper and lower limits
    x_upper = nanmax(X_T,[],'all');
    x_lower = nanmin(X_T,[],'all');
    y_upper = nanmax(Y_T,[],'all');
    y_lower = nanmin(Y_T,[],'all');
    
    %% Round and set min/maxs
    x_max = Xround*ceil(x_upper/Xround);
    x_min = Xround*floor(x_lower/Xround);
    y_max = Yround*ceil(y_upper/Yround);
    y_min = Yround*floor(y_lower/Yround);

end