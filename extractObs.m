%ABOUT
%   This function pulls out observations and their time-indices for a
%   specified variable, within a water year (WY). The data can either be 
%   pulled at a regular interval or at user specified timestamps
%   NOTE: assumes your overall dataset is hourly
%
%RELEASE NOTES
%   Written by Eric Smyth July 2018
%   esmyth@esmyth.com
%
%SYNTAX
%   [iObs,tObs,SDObs,NObs] = extractObs(obsMethod,obsData,obsVar,WY)
%
%****INPUTS****
% obsMethod 
%   1: Observations will be pulled at a regular interval within the WY
%   2: Observations will be pulled at user specified timestamps
% obsData
%   If obsMethod=1: the interval, in WEEKS, for obs to be pulled
%   If obsMethod=2: an array of matlab serial dates for obs to be pulled
%       NOTE: make sure these dates are within the WY
%       NOTE: make sure the dates are in chronological order
% obsVar
%   The variable from which you want to sample, to create spaced-out
%   observations. Needs to be a text string, e.g. 'SWE'
%       NOTE: make sure this variable exists as a fieldname of the WY
%       structure that you pass into the function
% WY
%   A structure of data for one water year, that includes at least a TIME
%   field and the variable field that you specify with obsVar. See "getWY"
%   function.
%
% ****OUTPUTS****
% iObs
%   Array of indices for the observations, within the WY
% tObs
%   Array of TIME data for the observations, within the WY
% SDObs
%   Rename this to what you want - will be an array of the extracted
%   observations themselves
% NObs
%   Integer, the number of observations extracted

function [iObs,tObs,SDObs,NObs] = extractObs(obsMethod,obsData,obsVar,WY)

Nt = length(WY.TIME(:,1));

if obsMethod == 1
    %%% Pull out observations at regular intervals
    %Convert weeks to hours
    obsWindowHours = obsData * 7 * 24;
    %Create array of observation indices
    iObs = (obsWindowHours+1 : obsWindowHours : Nt)';
    %Calc # of obs
    NObs = length(iObs);
    %Pull out obs for specified variable
    SDObs = WY.(char(obsVar))(iObs);
    %Pull out timestamps for these indices
    tObs = WY.TIME(iObs,:);
    
elseif obsMethod == 2
    %%% Pull out observations at user specified intervals
    %Create array of observation indices
    c = ismember(WY.TIME(:,7),obsData);
    iObs = find(c);
    %Calc # of obs
    NObs = length(iObs);
    %Pull out obs for specified variable
    SDObs = WY.(char(obsVar))(iObs);
    %Pull out timestamps for these indices
    tObs = WY.TIME(iObs,:);

end