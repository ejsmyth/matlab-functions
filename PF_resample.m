function [resample] = PF_resample(resampleMethod,weights)
%RELEASE NOTES
%   Written by Eric Smyth February 2018
%   esmyth@esmyth.com
%
%SYNTAX
%   [resample] = PF_resample(resampleMethod,weights)
%
% ****INPUTS****
% resampleMethod: integer from 1-4
%   1: Probabalistic resampling method
%   2: Residual resampling method
%   3: Stochastic universal resampling method
%   4: Metropolis-Hastings algorithm
%   See van Leeuwen (2009) pg. 4095 onwards ("Particle Filtering in
%   Geophysical systems")
%
% weights
%   An array of weights that correspond to model particles. These weights
%   will be used to resample from the particles. Can be any length
%
% ****OUTPUTS****
% resample
%   This is an array of indices, which corresponds with the array of input 
%   weights (and your array of model states at this point). For example, an
%   index of "4" means that the 4th model state (out of your ensemble of
%   model states) has been chosen. 

%% Initialize empty matrices
resample = [];
additionalResample = [];
randoms = [];

%% Generate discrete cumulative sum array of weights
cumSum = cumsum(weights);

%% Begin resampling routines
if resampleMethod == 1 % Probabalistic resampling
    for ii = 1:length(weights)
        randNum = rand;
        index = find(randNum <= cumSum,1);
        resample(ii) = index;
    end
    
elseif resampleMethod == 2 % Residual resampling
    multiples = weights .* length(weights);
    integerParts = nan(1,length(weights));
    for ii = 1:length(weights)
        step1 = multiples(ii);
        step2 = floor(step1);
        integerParts(ii) = step2;
        if step2 > 0
            copies = ii .* ones(1,step2);
            resample = [resample copies];
        end
    end
    resampleLength = length(resample);
    if resampleLength < length(weights)
        step3 = multiples - integerParts;
        remainderCDF = cumsum(step3);
        remainderSum = max(remainderCDF);
        residual = length(weights) - resampleLength;
        for gg = 1:residual
            randNum = rand * remainderSum;
            index = find(randNum <= remainderCDF,1);
            additionalResample(gg) = index;
        end
        resample = [resample additionalResample];
    end
    
elseif resampleMethod == 3 % Stochastic universal sampling
    randNum = rand * (1/length(weights));
    segments = linspace(randNum,randNum+((1/length(weights))*(length(weights)-1)),length(weights));
    for ii = 1:length(weights)
        index = find(segments(ii) <= cumSum,1);
        resample(ii) = index;
    end
    
elseif resampleMethod == 4 % Metropolis-Hastings algorithm
    firstWeight = max(weights);
    firstIndex = find(firstWeight == weights, 1);
    resample(1) = firstIndex;
    if firstIndex > 1
        if firstIndex < length(weights)
            remainingWeights = [weights(1:(firstIndex-1)) weights((firstIndex+1):end)];
            remainingIndices = [1:(firstIndex-1) (firstIndex+1):length(weights)];
        else
            remainingWeights = weights(1:(end-1));
            remainingIndices = 1:(length(weights)-1);
        end
    else
        remainingWeights = weights(2:end);
        remainingIndices = 2:length(weights);
    end
    for ii = 1:length(remainingWeights)
        currWeight = remainingWeights(ii);
        if ii == 1
            compWeight = firstWeight;
        end
        if currWeight > compWeight
            resample(ii+1) = remainingIndices(ii);
            compWeight = currWeight;
        else
            ratio = currWeight/compWeight;
            randNum = rand;
            if randNum < ratio
                resample(ii+1) = remainingIndices(ii);
                compWeight = currWeight;
            else
                resample(ii+1) = resample(ii);
            end
        end
    end
else
    error('Must choose resampleMethod of 1-4')
end

end