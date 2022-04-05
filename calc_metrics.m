% Generates a linear regression line, given x and y coordinates
%
% RELEASE NOTES
%   Written by Eric Smyth (esmyth@esmyth.com), Nov 2020
%
% SYNTAX
%    [RMSE, MAE, Rsq, NSE] = calc_metrics(X,Y)
%
% INPUTS
%   X = array of numbers
%   Y = array of numbers
%       NOTE: arrays must be the same length
%       NOTE: NaN values are removed, along with corresponding value from
%       other array
%       NOTE: for the NSE metric, the two arrays are not interchangable - in that context, X
%       should be the "observed" data and Y should be the "simulated" data
%
% OUTPUTS
%   RMSE = root mean squared error between the two arrays
%   MAE = mean absolute error
%   Rsq = r squared of fit of a linear regression line, fitting the two arrays
%   NSE = Nash Sutcliffe efficiency

% %% For testing
% X = sort(rand(20,1));
% Y = 2 + 3*X + randn(size(X));
% Y = X + randn(size(X))./10;
% scatter(X,Y,'ko')
% Y = 2 + 3*X;
% Y = 2+X;

function [RMSE, MAE, Rsq, NSE] = calc_metrics(X,Y)

    %% Make both arrays vertical for consistency, if not already
    if length(X(1,:)) > 1
        X = X';
    end
    if length(Y(1,:)) > 1
        Y = Y';
    end
    
    %% Remove NaN values
    X(isnan(Y)) = [];
    Y(isnan(Y)) = [];
    Y(isnan(X)) = [];
    X(isnan(X)) = [];
    if isempty(X)
        error('X or Y array was all NaN')
    end

    
    %% RMSE
    RMSE = sqrt(nanmean((X-Y).^2));
    
    %% MAE
    MAE = nanmean(abs(X-Y));
    
    %% NSE
    NSE = nashsutcliffe(X,Y);
    
    %% Rsq
    [~,~,~,Rsq] = linear_regression(X,Y);

end