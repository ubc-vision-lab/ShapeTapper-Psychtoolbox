%HELP for function CREATECONDITIONORDER
%
%CREATECONDITIONORDER creates a user specified number of randomly permuted
%orderings of a user-specified number of conditions.
%
%@PARAMS:
%           'conditions' - the total number of unique conditions to be permuted.
%           'n' - the number of times each unique condition occurs in the trial order. 'n'
%               should be divisable by 'r'.
%           'r' - the maximum number of back-to-back repetitions of any given unique
%               condition that is permitted to occur in the condition order. r=1 would
%               mean no repeats permitted.
%           'fNm' - the file name to store the order as.

function [A] = CreateConditionOrder(conditions,n,r,fNm)

%pre-allocate the condition order array
A = zeros(conditions*n,1);

switch nargin
    case 0
        disp(['Error!!! CREATECONDITIONORDER: The parameters ''conditions'', ''n'', '...
            'and ''r'' are required. Returned NaNs.'])
        return;
    case 1
        disp(['Error!!! CREATECONDITIONORDER: The parameters ''n'' and '...
            '''r'' are required. Returned NaNs.'])
        return;
    case 2
        disp(['Error!!! CREATECONDITIONORDER: The parameter ''r'' is required. Returned'...
            ' NaNs.'])
        return;
    case 3
        fNm = [];     
end

%check to make sure 'n' >= 'r'
if n < r
	disp(['Error!!! CREATECONDITIONORDER: The parameter ''n'' must be greater than or'...
            ' equal to the parameter ''r''. Returned NaNs.']);
        return;
end

%if 'r' does not equal zero, then check to make sure 'n' is divisable by 'r'.
if r == 0
	disp(['Note CREATECONDITIONORDER: The parameter ''r'' must be greater than zero. '...
        '''r'' is now set to 1.']);
        r=1;
elseif rem(n,r) ~=0
	disp(['Error!!! CREATECONDITIONORDER: The parameter ''n'' must be evenly divisable'...
        ' by the parameter ''r''. Returned NaNs.']);
        return;
end

%check the parameter 'fNm'
if ~ischar(fNm) && ~isempty(fNm)
    disp(['Error!!! CREATECONDITIONORDER: The parameter ''fNm'' must be either a string'...
        ' or empty. Returned NaNs.'])
    return;
end

%populate the condition order with 'n/r' iterations of unique condition orders
maxNumSegments = n/r;
currSegNum = 1;
while currSegNum <= maxNumSegments
    
    newSeg = [];
    
    %generate as many new segments of the condition order as repeats will permit. This
    %will result in a segment that will not violate the repeat max in and of itself when
    %randomized below.
    for i=1:r
        newSeg = [newSeg randperm(conditions)];
    end
    
    %randomize this segment of the condition order
    newSeg = transpose([1:length(newSeg); newSeg]);
    newSeg = sortrows(newSeg,[1 2]);
    
    %drop the randomizer column
    newSeg(:,1) = []; 
    
    %test this new segment for violations of the repeat limit at the joint between
    %segments. Do this by creating an array of the violating sequence to test against.
    %Start after the first segment has been generated.
    if currSegNum > 1
        
        %define the probe to be r+1 run of condition values that correspond to the first
        %condition value in the new segment
        probe = newSeg(1)*ones(r+1,1);
        
        repeatLimitExceeded = false;
        
        %test the probe against a test section
        for i=1:r
            
            %define the test section    
            testSec = [oldSeg(end-(r-1):end);newSeg(1:r)];
            testSec = testSec(i:i+r);
            
            %if the sum of the logical test of the probe against the test segment is equal
            %to the length of the probe, then the repeat limit has not been violated
            %and we can add the new segment to existing segments
            if sum(probe==testSec)==r+1
                repeatLimitExceeded=true;
            end
        end
        
        %if there is no violation, then add the new segment to the condition order
        if ~repeatLimitExceeded
        	A((conditions*r*(currSegNum-1))+1:conditions*r*currSegNum,1) = newSeg;
            
            %increment the current segment number
            currSegNum=currSegNum+1;
                
            %store 'newSeg' as 'oldSeg'
            oldSeg = newSeg;
        end
    else
        %add 'newSeg' to the condition order, store 'newSeg' as 'oldSeg', and increment
        %the current segment number
        A(1:conditions*r*currSegNum,1)=newSeg; 
        oldSeg = newSeg;
        currSegNum=currSegNum+1;      
    end
end

%create the trial order
%A = reshape(A,prod(size(A)),1);

%if the user wishes to store the results to disk, then do so.
if ~isempty(fNm)
    fid = fopen([fNm '.txt'],'w');
    fclose(fid); 
    dlmwrite([fNm '.txt'],A,'-append','delimiter','\t','newline','pc');   
end
return;
end

