function [WY,Nt] = getWY(waterYear,MASTER)

%RELEASE NOTES
%   Written by Eric Smyth July 2018
%   esmyth@esmyth.com
%
%SYNTAX
%   [WY,Nt] = getWY(waterYear,MASTER)
%
% ****INPUTS****
% waterYear 
%   water year integer (e.g. 2013) that is within your larger dataset
%
% MASTER
%   A structure with fields that have data for multiple years, that you
%   wish to shorten to one water year
%   Must include an Lx7 TIME array
%
% ****OUTPUTS****
% WY
%   A structure with the same fields as "MASTER" with data from selected
%   water year
% Nt
%   Length of shortened WY dataset

% Convert years to Water Years to pull out correct time period
WYcolumnLOG = MASTER.TIME(:,2) > 9;
WYcolumn = MASTER.TIME(:,1) + WYcolumnLOG;
iWY=any(ones(length(WYcolumn),1)*waterYear==WYcolumn*ones(1,length(waterYear)),2);

% Pull out the relevant inputs for the water years,
WY = [];

names = fieldnames(MASTER);
for i=1:length(names)
    eval(['WY.' names{i} '=MASTER.' names{i} '(iWY,:);'])
end

Nt = length(WY.(names{1})(:,1));

end