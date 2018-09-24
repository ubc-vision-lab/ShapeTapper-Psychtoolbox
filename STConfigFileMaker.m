%HELP for function STCONGIFFILEMAKER
%
%STCONGIFFILEMAKER creates a configuration file for Matlab-based ShapeTapper. It requires
%STCONFIGFILEMAKERASSIST and CREATECONDITIONORDER to operate. It relies on
%CREATECONDITIONORDER to generate a condition order which requires the user to specify the
%total number of unique conditions, the number of times each unique condition is to be
%administered, and the trial-span over which unique trials must occur.
%
%@PARAMS:
%           'conditions' - the total number of unique conditions to be permuted.
%           'n' - the number of times each unique condition occurs in the trial order per
%               block. 'n' should be divisable by 'r'.
%           'r' - the trial-span over which unique trials must occur. e.g., '1' would mean
%               that unique conditions must always occur from one trial to the next,
%               whereas '2' would mean that unique conditions can occur after 2 trials,
%               and permits a condition to repeat from one trial to the next.
%           'blocks' - the number of blocks in the file.
%           'eventsPerTrial' - the number of events per trial.
%           'fNm' - the name of the file to be stored to disc
%
%@OUTPUT:
%           'A' - the config file. Rows reflect fields and and the first column is the
%               header. Additional columns reflect trials and trial events.
%               Headers. Column headers to be saved to the file:
%               1: block_num
%               2: is_practice_block
%               3: block_has_feedback
%               4: block_feedback_thresh
%               5: trial_num
%               6: trial_max_time
%               7: trial_feedback
%               8: stim_img_name
%               9: stim_onset
%               10: stim_duration
%               11: stim_cent_x
%               12: stim_cent_y
%               13: stim_size_x
%               14: stim_size_y
%               15: stim_rotation
%               16: stim_is_target
%               17: mask_onset
%               18: mask_duration
%               19: mask_size
%               20: mask_color
%               21: mask_rotation
%               22: mask_fit
%               23: mask_margin
%               24: mask_touchable

function [A] = STConfigFileMaker(conditions,n,r,blocks,eventsPerTrial,fNm)

IMG = 'IMG_'; %part of the dummy image label 
EVNT = '_E_'; %image label denoting the event. Event number TBD.

[headers, formatSpec] = STConfigFileMakerAssist; %obtain the headers
numFields = length(headers);

%the number of trials per block
trialsPerBlock = conditions*n;

%total number of trials (will differ from the number of rows in the config file)
trialTot = trialsPerBlock*blocks;

%total number of rows in the config file, including the header row
rowsTotCnfgFl = trialTot*eventsPerTrial+1;

%we'll use FPRINT, so rows in 'A' will end up as columns (and cols in 'A' end up rows)
%when 'A' is saved to disk.
A = cell(numFields,rowsTotCnfgFl);

%the condition order of 1-event trials only. Cols 1: condition order; 2: block number; 3:
%trial number within the block
conditionOrder = zeros(trialTot,3);

%loop through each block
for i = 1:blocks
    
    %create a base-case condition order
    conditionOrder(1+(i-1)*trialsPerBlock:i*trialsPerBlock,1) = CreateConditionOrder(conditions,n,r);
    conditionOrder(1+(i-1)*trialsPerBlock:i*trialsPerBlock,2) = i;
    conditionOrder(1+(i-1)*trialsPerBlock:i*trialsPerBlock,3) = 1:trialsPerBlock;
end

%column index for A
ACol = 2;

%populate 'A' for FPRINTF
for i = 1:trialTot
    for j=1:eventsPerTrial
  
        A{1,ACol} = uint32(conditionOrder(i,2));    %1: block number
        A{2,ACol} = uint32(0);                      %2: is_practice_block?
        A{3,ACol} = uint32(0);                      %3: block_has_feedback?
        A{4,ACol} = 0;                              %4: block_feedback_thresh
        A{5,ACol} = uint32(conditionOrder(i,3));    %5: trial_num
        A{6,ACol} = 0;                              %6: trial_max_time
        A{7,ACol} = uint32(0);                      %7: trial_feedback
        
        %determine the image dummy-name
        if conditionOrder(i,1) < 10
            imgNmCndStr = ['00' num2str(conditionOrder(i,1))];
        elseif conditionOrder(i,1) > 9 && conditionOrder(i,1) < 100
            imgNmCndStr = ['0' num2str(conditionOrder(i,1))];
        else
            imgNmCndStr = num2str(conditionOrder(i,1)); 
        end
        
        if j < 10
            imgNmEvntStr = ['00' num2str(j)];
        elseif j > 9 && j < 100
            imgNmEvntStr = ['0' num2str(j)];
        else
           imgNmEvntStr = num2str(j); 
        end
        A{8,ACol} = [IMG imgNmCndStr EVNT imgNmEvntStr];   %8: stim_img_name
       
        A{9,ACol} = 0;                                  %9: stim_onset
        A{10,ACol} = 0;                                 %10: stim_duration
        A{11,ACol} = 0;                                 %11: stim_cent_x
        A{12,ACol} = 0;                                 %12: stim_cent_y
        A{13,ACol} = 0;                                 %13: stim_size_x
        A{14,ACol} = 0;                                 %14: stim_size_y
        A{15,ACol} = 0;                                 %15: stim_rotation
        A{16,ACol} = uint32(0);                         %16: stim_is_target
        A{17,ACol} = 0;                                 %17: mask_onset
        A{18,ACol} = 0;                                 %18: mask_duration
        A{19,ACol} = 0;                                 %19: mask_size
        A{20,ACol} = 0;                                 %20: mask_color
        A{21,ACol} = 0;                                 %21: mask_rotation
        A{22,ACol} = 0;                                 %22: mask_fit
        A{23,ACol} = 0;                                 %23: mask_margin
        A{24,ACol} = uint32(0);                         %24: mask_touchable
        
        %increment ACol
        ACol=ACol+1;
    end
end

%acquire the standard set of column headers
[A(:,1),headerLine,formatSpec] = STConfigFileMakerAssist;

%pick a directory to store them in, navigate to that directory
oldDir = pwd;
targetDir = uigetdir();
cd(targetDir);

%save the header line to file, store 'A' to file, then close the file.
fid = fopen([fNm '.csv'],'w');
fprintf(fid,headerLine);
fprintf(fid,formatSpec,A{:,2:end});
fclose(fid);

%return to the original directory
cd(oldDir);
return;



