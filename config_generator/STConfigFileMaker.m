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
%               whereas '2' would mean that unique conditions are permitted to occur every
%                2 trials which allows a unique condition to repeat.
%           'blocks' - the number of blocks in the file.
%           'eventsPerTrial' - the number of events per trial.
%           'fNm' - the name of the file to be stored to disc
%
%@OUTPUT:
%           'A' - the config file. Rows reflect fields and and the first column is the
%               header. Additional columns reflect trials and trial events.
%               Headers. Column headers to be saved to the file:
%               1: block_num %u
%               2: is_practice_block (0="no",1="yes") %u
%               3: block_has_feedback (0="no",1="yes") %u
%               4: block_feedback_thresh (%) %f
%               5: trial_num
%               6: trial_max_time
%               7: trial_feedback
%               8: background colour (a string "white" or "black") %s
%               9: stim_img_name
%               10: stim_onset
%               11: stim_duration
%               12: stim_cent_x (pixels)
%               13: stim_cent_y (pixels)
%               14: stim_size_x (cm)
%               15: stim_size_y (cm)
%               16: stim_rotation (degrees)
%               17: stim_is_touchable (0=not touchable; 1=touchable)
%               18: stim_is_target
%               19: subj_fixation_type
%               20: subj_fixation_onset (seconds)
%               21: subj_fixation_duration (seconds)
%               22: mask_onset (seconds)
%               23: mask_duration (seconds)
%               24: mask_size (cm)
%               25: mask_color (0=block;1=white;0-1 in RGB 3-element vector)
%               26: mask_rotation (degrees)
%               27: mask_fit (0=mask is bounding rectangle; 1=fitted along shape borders
%                       of stim image)
%               28: mask_margin (margin of mask points around stim image)

function [A] = STConfigFileMaker(conditions,n,r,blocks,eventsPerTrial,fNm)

IMG = 'IMG_'; %part of the dummy image label 
EVNT = '_E_'; %image label denoting the event. Event number TBD.

[FIELDNAMES,HEADERLINE,SAVEFORMATSPEC,~] = STConfigFileMakerAssist; %obtain the headers
numFields = length(FIELDNAMES);

%the number of trials per block
trialsPerBlock = conditions*n;

%total number of trials (NOTE: this value might differ from the number of rows in the config file)
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
        A{8,ACol} = BCKGRNDCLRNM;                   %8: background colour name
        
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
        A{9,ACol} = [IMG imgNmCndStr EVNT imgNmEvntStr];   %9: stim_img_name
       
        A{10,ACol} = 0;                                 %10: stim_onset
        A{11,ACol} = 0;                                 %11: stim_duration
        A{12,ACol} = 0;                                 %12: stim_cent_x (pixels)
        A{13,ACol} = 0;                                 %13: stim_cent_y (pixels)
        A{14,ACol} = 0;                                 %14: stim_size_x (cm)
        A{15,ACol} = 0;                                 %15: stim_size_y (cm)
        A{16,ACol} = 0;                                 %16: stim_rotation (degrees)
        A{17,ACol} = uint32(0);                         %17: stim_is_target
        A{18,ACol} = uint32(0);                         %18: stim_is_touchable
        A{19,ACol} = uint32(0);                         %19: subj_fixation_type
        A{20,ACol} = 0;                                 %20: subj_fixation_onset (seconds)
        A{21,ACol} = 0;                                 %21: subj_fixation_duration (seconds)
        A{22,ACol} = 0;                                 %22: mask_onset (seconds)
        A{23,ACol} = 0;                                 %23: mask_duration (seconds)
        A{24,ACol} = 0;                                 %24: mask_size
        A{25,ACol} = 0;                                 %25: mask_color
        A{26,ACol} = 0;                                 %26: mask_rotation (degrees)
        A{27,ACol} = uint32(0);                         %27: mask_fit
        A{28,ACol} = 0;                                 %28: mask_margin
        
        %increment ACol
        ACol=ACol+1;
    end
end

%add the filed name header to 'A'
A(:,1) = transpose(FIELDNAMES);

%pick a directory to store them in, navigate to that directory
oldDir = pwd;
targetDir = uigetdir();
cd(targetDir);

%save the header line to file, store 'A' to file, then close the file.
fid = fopen([fNm '.csv'],'w');
fprintf(fid,HEADERLINE);
fprintf(fid,SAVEFORMATSPEC,A{:,2:end});
fclose(fid);

%return to the original directory
cd(oldDir);
return;



