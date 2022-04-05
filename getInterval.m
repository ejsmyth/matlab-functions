function [iIntvl] = getInterval(kk,iObs,Nt,windowLength,numWindows)

%ABOUT
%   This function is designed to return indices for a particle filter or
%   particle batch smoother "window" within a water year (WY). One would use these
%   indices to pull out forcing data for a model to use within this window.
%
%RELEASE NOTES
%   Written by Eric Smyth February 2018
%   esmyth@esmyth.com
%
%SYNTAX
%   [iIntvl] = getInterval(kk,iObs,Nt,windowLength,numWindows)
%
%****INPUTS****
% kk 
%   Integer that represents the window number, up to the total number of
%   windows within the WY (i.e. the 3rd window out of 26 windows). This 
%   should be the index of a "for" loop.
% iObs
%   An array of indices within the WY that indicate when observations
%   occur. With the particle filter / PBS, windows are defined by when
%   these observations occur.
% Nt
%   The total number of timestamps within the WY.
% windowLength
%   The number of observations within each window. Ex: windowLength of 1 is
%   the normal particle filter.
% numWindows
%   The total number of windows within the modeing period
%
% ****OUTPUTS****
% iIntvl
%   Array of indices for the current window, within the WY.


if kk==1
    iIntvl=1:iObs(windowLength);
elseif kk<numWindows
    iIntvl=iObs(windowLength*kk-windowLength)+1:iObs(windowLength*kk);
elseif kk==numWindows
    iIntvl=iObs((kk-1)*windowLength)+1 : Nt;
end

return