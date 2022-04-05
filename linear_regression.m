% Generates a linear regression line, given x and y coordinates
%
% RELEASE NOTES
%   Written by Eric Smyth (esmyth@esmyth.com), Nov 2020
%
% SYNTAX
%    [Yhat,Yint,slp,Rsq] = linear_regression(X,Y)
%
% INPUTS
%   X = array of numbers, what might be the x-coordinates on a 2D grid
%   Y = array of numbers, same length as X, that might be y-coordinates on a 2D grid
%
% OUTPUTS
%   Yhat = array of numbers, same length as X, that sit along / define /
%       can be plotted as a line, that minimizes squared errors
%   Yint = y-intercept of regression line
%   slp = slope of regression line
%   Rsq = r squared of fit

%%% For testing
% X = sort(rand(20,1));
% Y = 2 + 3*X + randn(size(X));
% Y = 2 + 3*X;

function [Yhat,Yint,slp,Rsq] = linear_regression(X,Y)

    %%% Make arrays vertical, if not already
    if length(X(1,:)) > 1
        X = X';
        xflag = 1;
    else
        xflag = 0;
    end
    if length(Y(1,:)) > 1
        Y = Y';
    end

    %%% Pad the x coordinates with ones
    M = [ones(length(X),1),X];
    
    %%% Use the "normal" equations
    coef = inv(M'*M)*M'*Y; % coef contains regression estimates of the parameters
    Yhat = M*coef; % A new array, with y values defined by the regression coef

    %%% Grab slope and y-intercept
    Yint = coef(1);
    slp = coef(2);

    %%% Calculate r squared
    Rsq = 1 - sum((Y - Yhat).^2)/sum((Y - mean(Y)).^2);

    %%% Invert if necessary
    if xflag == 1
        Yhat = Yhat';
    end

end