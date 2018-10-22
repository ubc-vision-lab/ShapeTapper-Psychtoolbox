%HELP for function CREATETRIALORDER
%
%CREATETRIALORDER creates a tx(f+1) (trials x factors + the trial number column) matrix
%in which the first column denotes the trial number and the additional columns denotes the
%particular level of a given factor administered for that trial.
%
%@PARAMS:
%           'factors' - the number of elements corresponds to the number of factors and 
%               the value in each element denotes the number of levels in that factor.
%           'n' - the number of times each unique condition occurs in the trial order
%           'r' - the maximum number of back-to-back repetitions of any given unique
%               condition that is permitted to occur in the condition order. r=1 would
%               mean no repeats permitted.
%           'fNm' - the file name. If empty then the trial order is not saved to disk.
%           'fExt' - the file extension of the file to be saved to disk.
%           'delim' - the desired delimiter.
%
%@OUTPUT:
%           'A' - a tx(f+1)  matrix in which the first column reflects the trial number
%               and the additional columns denotes the particular levels of the factors
%               to be administered on that trial.
%

function [A] = CreateTrialOrder(factors,n,r,fNm,fExt,delim)

switch nargin
    case 0
        disp(['Error!!! CREATECONDITIONORDER: The parameters ''factors'', '...
            '''n'', and ''r'' are required. Returned NaNs.'])
        return;
    case 1
        disp(['Error!!! CREATECONDITIONORDER: The parameters ''n'' and, '...
            '''r'' are required. Returned NaNs.'])
        return;
    case 2
        disp(['Error!!! CREATECONDITIONORDER: The parameters ''r'' is '...
            'required. Returned NaNs.'])
        return;
    case 3
        fNm = [];
        fExt = [];
        delim = [];
        disp(['Note! CREATECONDITIONORDER: The parameters ''fNm'', ''fExt'', and ''delim'''...
            ' are missing. No file will be saved to disk.'])
    case 4
        fExt = '.txt';
        delim = '\t';
        disp(['Note! CREATECONDITIONORDER: The parameters ''fExt'', and ''delim'' are '...
            'missing. Defaulting to ''.txt'' and a tab delimiter.'])
    case 5
        delim = '\t';
        disp(['Note! CREATECONDITIONORDER: The parameter ''delim'' is missing. Defaulting '...
            'to a tab delimiter.'])
end

%reject condition totals that are more than 10
if sum(factors>10) > 0
    disp(['Note! CREATECONDITIONORDER: The parameter ''factors'' should '...
        'contain elements 10 or less. Returned NaNs.'])
    return;
end

%check the parameter 'fNm'
if ~ischar(fNm) && ~isempty(fNm)
	disp(['Error: CREATECONDITIONORDER!!! The parameter ''fNm'' '...
            ' must be either a string or empty. Returned NaNs.'])
return;
end

%determine some critical variables.
totConds = prod(factors);
totTrials = totConds*n;

%create the trial order but add a column for the trial number
A = nan(totTrials,length(factors)+1);
A(:,1) = 1:totTrials;

%factors = sort(factors,'descend');

%loop through each factor (the number of factors corresponds to the length of 'factors'
for i = 1:length(factors)
    
    %temporary matrix to hold the order for the current factor
    temp = [];

    %generate a range of integers reflecting the levels of the current
    %factor
    range = 1:factors(i);

    %generate all possible permutations (each unique permutation is a row)
    permPop = perms(range);
    
    %the number of permuted sequences to be sampled
    N = totTrials/factors(i);
    
    count = 1;
    while count <= N
        
        %select one random permutation from the population of permutations.
        %Note that this process is occuring with replacement
        permSmpl = transpose(permPop(randsample(length(permPop(:,1)),1),:));
        
        %check for repeats at the joints (the code needs to be updated)
        if r < 2 && count > 1 && temp(end) == permSmpl(1)
    
        else
            temp = [temp; permSmpl];
            count = count+1;
        end
    end
    
    A(:,i+1) = temp;
    A = sortrows(A,i+1);
end

A = sortrows(A);

%if the user wishes to store the results to disk, then do so.
if ~isempty(fNm)
    fid = fopen([fNm fExt],'w');
    fclose(fid); 
    dlmwrite([fNm fExt],A,'-append','delimiter',delim,'newline','pc');   
end
return;
end

