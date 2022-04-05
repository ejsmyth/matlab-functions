function [P_daily, TIME_daily] = precip2storms(P, TIME, time_format, waterYear, startDate, endDate, Precip_thresh)

% RELEASE NOTES
%   Written by Eric Smyth March 2019
%   esmyth@esmyth.com
%
% ****SYNTAX****
%   [P_daily, TIME_daily] = precip2storms(P, TIME, time_format, startDate, endDate, Precip_thresh)
%
% ****INPUTS****
% P
%   Array of hourly or daily precip
% TIME
%   Standard Nx7 TIME matrix (see time_builder function)
% time_format
%   1 = using hourly data (will be converted to daily)
%   2 = using daily data already
% waterYear
%   Integer of the water year that we are dealing with (e.g. 2016)
% startDate
%   matlab serial date of when you want to begin looking at storms
% endDate
%   matlab serial date of when you want to end looking at storms
% Precip_thresh
%   value of daily precip which can be considered part of a storm. In other
%   words, for a given day with precipitation above this threshold, if the
%   following/previous day also has precip above the threshold, the precips
%   will be summed into a larger storm (pay attention to units)
%
% ****OUTPUTS****
% P_daily
%   New array of daily precip, where precip is summed within storms onto
%   the last day of the storm, and other days within the storm have zero.
%   Will only be as long as your start-end dates
% TIME_daily
%   Standard TIME matrix, only be as long as your start-end dates

%% Setup
startIndex = find(TIME(:,7) == startDate,1);
endIndex = find(TIME(:,7) == endDate,1);
P_section = P(startIndex:endIndex);
TIME_section = TIME(startIndex:endIndex,:);

%% Aggregate precip to daily, for better indication of biggest snowfall events
if time_format == 1
    [P_daily, TIME_daily] = aggMET(P_section, TIME_section, 1, 24, waterYear-1, 10, 1, 0, 0, 1, 5, 1);
elseif time_format == 2
    P_daily = P_section;
    TIME_daily = TIME_section;
else
    error('Must choose valid time format option - see syntax note in function')
end

%% Combine nearby days of snow into larger snow events
group_lock = 1;
while group_lock > 0
    try
        P_adj = P_daily - Precip_thresh;
        storm_record = group_events(P_adj, 1);

        %% Adjust storm record to work better below
        stormIndices = [];
        for ii = 1:length(storm_record(:,1))
            stormIndices = [stormIndices storm_record(ii,1):storm_record(ii,2)];
        end
        
        group_lock = 0;
        
    catch
        Precip_thresh = max(Precip_thresh - 1,0);
        if group_lock == Precip_thresh
            error('Couldnt aggregate storms')
        else
            group_lock = group_lock+1;
        end
    end
end

%% Adjust P_daily to sum up snowfall in storm events
P_record = P_daily;
for ii = 1:length(P_daily)
    if find(stormIndices == ii,1) > 0
        % This day is part of larger storm, might need to adjust
        if ii == 1
            % This is the first day
            if find(stormIndices == ii+1,1) > 0
                % Next day is part of same storm, set precip to zero
                P_record(ii) = 0;
            end
        elseif ii == length(P_daily)
            % This is the last day
            if isempty(find(stormIndices == ii-1,1)) == 1
                % Previous day not part of same storm, isolated snowfall
            else
                % Must be last day of larger storm
                lock_P = 1;
                temp_P = P_daily(ii);
                while lock_P >= 1
                    if find(stormIndices == ii-lock_P,1) > 0
                        temp_P = temp_P + P_daily(ii-lock_P);
                        lock_P = lock_P + 1;
                    else
                        lock_P = 0;
                        P_record(ii) = temp_P;
                    end
                end
            end
        else
            % Not first or last day
            if find(stormIndices == ii+1,1) > 0
                % Next day is part of same storm, set precip to zero
                P_record(ii) = 0;
            elseif isempty(find(stormIndices == ii-1,1)) == 1
                % Previous day not part of same storm, isolated snowfall
            else
                % Must be last day of larger storm
                lock_P = 1;
                temp_P = P_daily(ii);
                while lock_P >= 1
                    if find(stormIndices == ii-lock_P,1) > 0
                        temp_P = temp_P + P_daily(ii-lock_P);
                        lock_P = lock_P + 1;
                    else
                        lock_P = 0;
                        P_record(ii) = temp_P;
                    end
                end 
            end
        end
    end
end

P_daily = P_record;

end