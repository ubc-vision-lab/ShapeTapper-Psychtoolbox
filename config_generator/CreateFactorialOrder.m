%HELP for function CREATEFACTORIALORDER
%
%Author: Rob Whitwell
%
%CREATETRIALORDER creates a tx(f+1) (trials x factors + the trial number column) matrix
%in which the first column denotes the trial number and the additional columns denotes the
%particular level of a given factor administered for that trial.
%
%@PARAMS:
%           'factors' - the number of elements corresponds to the number of factors and 
%               the value in each element denotes the number of levels in that factor.
%           'n' - the number of times each unique condition occurs in the trial order
%           'rep' - the maximum number of back-to-back repetitions of any given unique
%               condition that is permitted to occur in the condition order. rep=1 would
%               mean a given pair of trials t and t+1 can contain the same level of a
%               factor. In short, level back-to-back-repeats are limited to 1.
%           'N' - the desired number of factorial trial orders 
%           'fNm' - the file name. If empty then the trial order is not saved to disk.
%           'fExt' - the file extension of the file to be saved to disk. Default is .txt.
%           'delim' - the desired delimiter. Default is tab delimited.
%
%@OUTPUT:
%           'A' - a trial (rows) x factors (columns) x number of orders (pages) 3D matrix
%               in which rows of all pages reflect trials; the first column of all pages
%               reflects the trial number and the remaining columns reflect the factors
%               in the order in which they appear in the parameter 'factors'. The factor
%               columns in A are integers denoting the levels in that factor.

function [A] = CreateFactorialOrder(factors,n,rep,N,fNm,fExt,delim)

MAXINVALIDORDERS = 100000;

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
        N=1;
        fNm = [];
        fExt = [];
        delim = [];
        disp(['Note! CREATECONDITIONORDER: The parameters ''fNm'', ''fExt'', and ''delim'''...
            ' are missing. No file will be saved to disk.'])
        N=1;
    case 4
        fNm = [];
        fExt = '.txt';
        delim = '\t';
        disp(['Note! CREATECONDITIONORDER: The parameters ''fExt'', and ''delim'' are '...
            'missing. Defaulting to ''.txt'' and a tab delimiter.'])        
    case 5
        fExt = '.txt';
        delim = '\t';
        disp(['Note! CREATECONDITIONORDER: The parameter ''delim'' is missing. Defaulting '...
            'to a tab delimiter.'])
    case 6
        delim = '\t';
        disp(['Note! CREATECONDITIONORDER: The parameter ''delim'' is missing. Defaulting '...
            'to a tab delimiter.'])
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

%create the trial order
A = zeros(totTrials,length(factors)+1,N);

%fill up the trial column
A(:,1,:)=repmat(1:totTrials,[1 1 N]);

%use the first factor to generate a random trial order
A(:,2,:)=CreateConditionOrder(factors(1),totTrials/factors(1),rep,N);

%create the user-defined number of factorial orders
for i=1:N
    
    invalidOrders = zeros(length(factors)-1); %tally the number of invalid Orders

    %loop through the remaining factors, filling up the trial order
    for j = 2:length(factors)

        finish = false; 

        %loop through the process until finished.
        while ~finish
            
            %determine the sorting order
            srt=[];  
            for k=1:j-1
                srt(k)=k+1;
            end
            srt=[srt 1];

            %sort the rows according to the previous factor
            A(:,:,i)=sortrows(A(:,:,i),srt);

            %determine the number of trials per unique conditions produced by the current
            %factor and the all previous ones
            nUniqueCnds = prod(factors(1:j-1));
            trialsPerUniqueCnds = totTrials/nUniqueCnds;

            %loop through the unique combinations of levels of the previous factors
            for k = 1:nUniqueCnds

                %create and populate an column array pool of levels 
                levelPool=zeros(trialsPerUniqueCnds,1);
                for m = 1:trialsPerUniqueCnds/factors(j)
                    levelPool(1+(factors(j)*(m-1)):m*factors(j),1) = randperm(factors(j),factors(j));
                end

                %add the permutations for this combination of previous factor level and
                %current factor to the trial order
                A(1+trialsPerUniqueCnds*(k-1):k*trialsPerUniqueCnds,j+1,i) = levelPool;

            end

            %sort the trial order according to the trial number
            A(:,:,i)=sortrows(A(:,:,i),1);

            %counter for the number of valid trials
            validTrials = rep;

            %trial counter
            curTrl = rep+1;

            %check each trial for a violation of the repetition limit
            while curTrl <= totTrials

                %the level in the current trial
                curLvl = A(curTrl,j+1,i);

                %construct the probe and the test piece
                probe = curLvl*ones(rep+1,1);
                testPiece = [A(curTrl-rep:curTrl-1,j+1,i); curLvl];

                %test the test piece. Testpiece fails if the next level completes a repeat
                %sequence that exceeds the repetition limit.
                if sum(testPiece==probe)==rep+1
                    curTrl=totTrials+1; %should exit the loop
                    invalidOrders(j-1) = invalidOrders(j-1)+1;
                else
                    validTrials=validTrials+1;
                    curTrl=curTrl+1;
                end
            end

            %if the trial order works, then move on to the next factor
            if validTrials == totTrials || invalidOrders(j-1) > MAXINVALIDORDERS
                finish = true; 
            end
        end
    end

    %if the user wishes to store the results to disk, then do so.
    if ~isempty(fNm)
        fid = fopen([fNm '_' num2str(i) fExt],'w');
        fclose(fid); 
        dlmwrite([fNm '_' num2str(i) fExt],uint32(A(:,:,i)),'-append','delimiter',delim,'newline','pc');   
    end
end
return;
end

