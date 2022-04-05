function [C] = generate_combinations(A,B,max_length)

%RELEASE NOTES
%   Written by Eric Smyth March 2019
%   esmyth@esmyth.com
%
% ****SYNTAX****
%   [C] = uniqueRows(A,B)
%
% ****INPUTS****
% A
%   Array from 1:N, e.g. 1:20
% B
%   Integer less than or equal to N, e.g. 10
% max_length
%   Integer value of the maximum length of a matrix that matlab can
%   "handle." For example, 30c10 produces a ~30 million row matrix that
%   matlab can't really handle. So max_length could be specified as
%   100,000. This is really only useful if you "do" something with these
%   100k long intervals while the code is running. Right now, this function
%   just combines all the intervals together at the end into one long
%   matrix
%
% ****OUTPUTS****
% C
%   This will be a matrix of B combinations within A, i.e. all the different
%   combinations that represent "A choose B"
%
% ****NOTE****
%   This function only works in finding the combinations of 1:N, i.e. a continuous
%   array from 1 to N. If you want to generate combinations of a different
%   array, like [4 6 2 9], just generate combinations of [1 2 3 4] and use
%   the values as indices to the [4 6 2 9] array.

%% Calculate length of full matrix
full_length = round(factorial(length(A))/(factorial(B)*factorial(length(A)-B)));

%% Calculate iterations
num_iterations = max(ceil(full_length/max_length),1);

%% Initialize matrix to hold combinations (or sections if necessary)
if num_iterations == 1
    C = nan(full_length,B);
else
    C = nan(max_length+1,B);
    C_store = [];
end

%% Generate combinations
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
    
    if num_iterations > 1
        end_store = C(end,:);
        C_store = [C_store; C];
        if lock_iteration < num_iterations
            C = nan(max_length+1,B);
        end
    end
    
    lock_iteration = lock_iteration + 1;
    
end

if num_iterations > 1
    C = C_store;
end

C(~any(~isnan(C), 2),:)=[];

end
