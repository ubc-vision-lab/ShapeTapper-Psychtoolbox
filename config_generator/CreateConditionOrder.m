%HELP for function CREATECONDITIONORDER
%
%Author: Rob Whitwell
%
%CREATECONDITIONORDER creates a user specified number of randomly permuted
%orderings of a user-specified number of conditions.
%
%@PARAMS:
%           'conditions' - the total number of unique conditions to be permuted.
%           'n' - the number of times each unique condition occurs in the trial order. 'n'
%               should be divisable by 'r'.
%           'rep' - the maximum number of back-to-back repetitions of any given unique
%               condition that is permitted to occur in the condition order.
%           'N' - the number of 
%           'fNm' - the file name to store the order as.

function [A] = CreateConditionOrder(conditions,n,rep,N,fNm)

A = NaN;

switch nargin
    case 0
        disp(['Error!!! CREATECONDITIONORDER: The parameters ''conditions'', ''n'', '...
            'and ''rep'' are required. Returned NaN.'])
        return;
    case 1
        disp(['Error!!! CREATECONDITIONORDER: The parameters ''n'' and '...
            '''rep'' are required. Returned NaN.'])
        return;
    case 2
        disp(['Error!!! CREATECONDITIONORDER: The parameter ''rep'' is required. Returned'...
            ' NaN.'])
        return;
    case 3
        N=1;
        fNm=[];
    case 4
        fNm=[];
end

%if 'rep' equals 0, set it to 1.
if rep == 0
	disp(['Note CREATECONDITIONORDER: The parameter ''r'' must be greater than zero. '...
        '''r'' is now set to 1.']);
        rep=1;
end

%check the parameter 'fNm'
if ~ischar(fNm) && ~isempty(fNm)
    disp(['Error!!! CREATECONDITIONORDER: The parameter ''fNm'' must be either a string'...
        ' or empty. Returned NaN.'])
    return;
end

%pre-allocate the condition order array
A = zeros([conditions*n N]);

%determine the total number of trials
totTrials = conditions*n;

%for each desired trial order
for i=1:N
    
    %row index for the temporary matrix
    row = 1;

    %first, determine the length of the first condition permutation piece. Use an
    %adjustment to the repetition limit in case the length of the piece is too long.
    adj = 0;

    %determine if the length of the new piece is too long and find the adjustment if
    %necessary
    while row + (conditions*(rep-adj)) - 1 > totTrials && adj < rep
        adj = adj+1;
    end

    %temporary placeholder for the random permutations
    newPerm = nan(rep-adj,conditions);

    %create the next series of random permutations to complete the next piece
    for j = 1:rep-adj
        newPerm(j,:) = randperm(conditions,conditions);
    end

    %convert the new permutation piece into a column vector
    newPerm = reshape(newPerm,[],1);

    %store length of the next piece of permuted condition order
    pieceLngth = length(newPerm);

    %produce new pieces of permutations provided there is room for them
    while row + pieceLngth - 1 <= totTrials

        %if this is not the first new piece, check to make sure the repetitions rule is
        %not violated at the joints between pieces
        if row > 1

            %create a matrix of probe sequences where rows reflect each trial in the probe
            %sequence. The test uses the trial-to-trial difference between condition numbers,
            %and so 0s would reflect no difference (i.e. repeats).
            probe = zeros(rep+1,1);

            %create the test piece
            testPiece = [0; diff([A(row-rep:row-1,i);newPerm(1:rep)])];

            %assume the new piece doesn't violate the repetition limit at the joint, but
            %update it if a violation is discovered
            violated = false;

            %testPiece row index
            indx = 1;
            while ~violated && indx+rep-1 < length(testPiece)

                %if the sum of the logical comparison of the violation probe against
                %the test piece is equal to the repeat limit + 1, then this
                %constitutes a violation
                if sum(probe(:,1)==testPiece(indx:indx+rep))==rep+1
                    violated = true;
                end

                %increment the 'testPiece' row index
                indx=indx+1;
            end

            %if the repeat sequence limit was not violated, add the test piece to factor
            %trial order and increment the row index.
            if ~violated
                A(row:row+pieceLngth-1,i)=newPerm;
                row = row+pieceLngth;  
            end

        %otherwise, this is the first piece, so add it to 'A', update 'row', and determine the
        %size of the next piece
        else        
            A(row:row+pieceLngth-1,i)=newPerm;
            row = row+pieceLngth; 
        end

        %find the length of the next piece provided there is still room in the condition
        %order.
        if row < totTrials
            %Use an adjustment to the repetition limit in case the default length of the next
            %piece is too long.
            adj = 0;

            %determine if the length of the new piece is too long and find the adjustment if
            %necessary
            while row + (conditions*(rep-adj)) - 1 > totTrials && adj < rep
                adj = adj+1;
            end

            %temporary placeholder for the random permutations
            newPerm = nan(conditions, rep-adj);

            %create the next series of random permutations to complete the next piece
            for j = 1:rep-adj
                newPerm(:,j) = randperm(conditions,conditions);
            end

            %convert the new permutation piece into a column vector
            newPerm = reshape(newPerm,[],1);

            %store length of the next piece of permuted condition order
            pieceLngth = length(newPerm);
        end
    end

    %if the user wishes to store the results to disk, then do so.
    if ~isempty(fNm)
        fid = fopen([fNm '_' num2str(i) '.txt'],'w');
        fclose(fid); 
        dlmwrite([fNm '_' num2str(i) '.txt'],[transpose(1:trialTot) A(:,i)],'-append','delimiter','\t','newline','pc');   
    end
end
return;
end

