function [obsData] = precip_observations(arg1,arg2,arg3,arg4,arg5,arg6,arg7,varargin)

%RELEASE NOTES
%   Written by Eric Smyth March 2019
%   esmyth@esmyth.com
%
% ****SYNTAX****
%   OPTION 1: Choose overall # of largest precip events
%   [obsData] = precip_observations(totalObsNum,P_daily,TIME_daily)
%
%   OPTION 2: Find largest precip events at least X days apart
%   [obsData] = precip_observations(totalObsNum,P_daily,TIME_daily,minNumDays,maxSamples,sampleIncrement,max_length)
%
% ****INPUTS****
% totalObsNum
%   Integer value of the total number of observations you want to make
%   within the defined time period
% P_daily
%   Array of daily precipitation values, ideally aggregated into storms
% TIME_daily
%   Standard Nx7 TIME array corresponding to P_daily (see time_builder.m)
% minNumDays
%   If option 2, the minimum number of days between observations
% maxSamples
%   If option 2, integer value, to limit the sample size to the top X days of
%   precipitation (e.g. 100)
% sampleIncrement
%   If option 2, if the top X snowfall days don't have combinations that are far 
%   enough apart, will increase sample size to 2X and re-run all combos
%   (e.g. 10)
% max_length
%   If option 2, maximum length of a matrix that matlab can "handle." When
%   the routine needs to generate a large number of potential day
%   combinations, it will only generate X combinations at a time, evaluate
%   them, and then discard. So this number should be large (to allow for
%   many possible combinations evaluated againt each other at once) but not
%   too large (e.g. 100,000)
%
% ****OUTPUTS****
% obsData
%   Array of matlab serial dates at which observations will be taken

if nargin == 3
    %% obs_option must be 1: Choose overall # of largest precip events
    totalObsNum = arg1;
    P_daily = arg2;
    TIME_daily = arg3;
    
    [~,Pindices] = maxk(P_daily,totalObsNum);
    Pindices = sort(Pindices);
    Pindices = Pindices + 5;           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    obsData = TIME_daily(Pindices,7);
    
elseif nargin == 7
    %% Option 2: Find largest precip events at least X days apart
    totalObsNum = arg1;
    P_daily = arg2;
    TIME_daily = arg3;
    minNumDays = arg4;
    maxSamples = arg5;
    sampleIncrement = arg6;
    max_length = arg7;
    
    [Pvalues,Pindices] = maxk(P_daily,maxSamples);

    if totalObsNum == 1
        % No need to go through all this...just choose top precip day
        obsData = TIME_daily(Pindices(1),7);
    else
        lock_fsm = 1; % initialize "lock" variable for while loop
        comboStore = [];
        while lock_fsm >= 1
            A = 1:lock_fsm*sampleIncrement;
            B = totalObsNum;
            currentPrecips = Pvalues(A);
            currentIndices = Pindices(A);

            %%% Calculate length of full combination matrix
            full_length = round(factorial(length(A))/(factorial(B)*factorial(length(A)-B)));
            %%% Calculate iterations (>=1)
            num_iterations = max(ceil(full_length/max_length),1);
            %%% Initialize matrix to hold combinations (or sections if necessary)
            if num_iterations == 1
                C = nan(full_length,B);
            else
                C = nan(max_length+1,B);
            end

            %%% Generate combinations
            lock_iteration = 1;
            while lock_iteration <= num_iterations
                for ii = 1:length(C(:,1))
                    if ii == 1
                        if lock_iteration == 1
                            C(ii,:) = 1:B;
                        else
                            C(ii,:) = end_store;
                        end
                    else
                        if C(ii-1,end) < A(end)
                            C(ii,1:end-1) = C(ii-1,1:end-1);
                            C(ii,end) = C(ii-1,end)+1;
                        else
                            lock_c = 1;
                            while  lock_c < B
                                if C(ii-1,end-lock_c) < A(end-lock_c)
                                    if lock_c < (B-1)
                                        C(ii,1:end-lock_c-1) = C(ii-1,1:end-lock_c-1);
                                    end
                                    C(ii,end-lock_c) = C(ii-1,end-lock_c)+1;
                                    C(ii,end-lock_c+1:end) = C(ii,end-lock_c)+[1:lock_c];
                                    lock_c = B;
                                else
                                    lock_c = lock_c + 1;
                                end
                            end
                        end
                    end
                end

                C(~any(~isnan(C), 2),:)=[];

                %%% Calc hueristics
                hueristics = nan(1,length(C(:,1)));
                for tt = 1:length(C(:,1))
                    % Precip score
                    tempPrecips = currentPrecips(C(tt,:));
                    precipScore = nansum(tempPrecips);
                    % Date score
                    tempIndices = currentIndices(C(tt,:));
                    tempDateNums = TIME_daily(tempIndices,7);
                    tempDateNums = sort(tempDateNums);
                    minDateDiff = nanmin(diff(tempDateNums));
                    if minDateDiff >= minNumDays
                        dateScore = 1;
                    else
                        dateScore = -1;
                    end
                    % Overall score
                    hueristics(tt) = precipScore * dateScore;
                end

                if nanmax(hueristics) > 0
                    % Routine found a good collection of days
                    findBest = find(hueristics == nanmax(hueristics),1);
                    bestIndices = C(findBest,:);
                    bestSectionIndices = Pindices(bestIndices);
                    try
                        bestSectionIndices = bestSectionIndices + 5; %%%%%%%%%%%%%%%%%%%%
                        bestDateNums = TIME_daily(bestSectionIndices,7);
                        obsData = sort(bestDateNums);
                        lock_fsm = -9999;
                        lock_iteration = 99999999;
                    catch
                        end_store = C(end,:);
                        C = nan(max_length+1,B);
                        lock_iteration = lock_iteration + 1;
                    end
                else
                    end_store = C(end,:);
                    C = nan(max_length+1,B);
                    lock_iteration = lock_iteration + 1;
                end 
            end

            lock_fsm = lock_fsm + 1;                    

        end
    end

else
    error('Invalid number of input arguments. See precip_observations.m for SYNTAX')
end

end