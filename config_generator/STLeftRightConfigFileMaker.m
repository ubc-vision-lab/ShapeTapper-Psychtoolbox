%HELP for function STLEFTRIGHTCONFIGMAKER. This function makes a 2-image per trial
%left/right (same/different) 2AFC configuration file for the Matlab-based ShapeTapper
%program and stores it to disk.
%
%Author: Rob Whitwell
%
%@PARAMS
%           'factors' - an array in which each element is a factor the value of which
%               specifies the number of levels in that given factor. Assumes the following
%               column structure 1:number of shapes; 2:number of unique orientations;
%               3:the quadrant of first image of the pair (1-4); 4:the same (1) vs.
%               different (2) nature of the pair of shapes.
%           'n' - the number of times each unique condition occurs in the trial order
%           'r' - the maximum number of back-to-back repetitions of any given unique
%               condition that is permitted to occur in the condition order. r=1 would
%               mean no repeats permitted.
%           'b' - the number of blocks of trials in the trial order
%           'scrnXPxls - the width of the display screen (in pixels)
%           'scrnYPxls - the height of the display screen (in pixels)
%           'dpi' - the dots per inch of the display screen
%           'maxDispImgSizeCm' - the maximum size of the image on the display screen (cm)
%           'cnfFNm' - the configuration file name. If empty then the trial order is not
%               saved to disk.

function [A] = STLeftRightConfigFileMaker(factors,ncp,repLim,blocks,scrnXPxls,scrnYPxls,dpi,maxDispImgSizeCm,cnfFNm)

%the number of events (images to present) per trial
EVENTSPERTRIAL = 2;

%image file extension
IMGEXT = '*.png';

%the background colour
BCKGRNDCLRNM = 'black';
STIMDUR = 10; %(1/30)*10;
TRIALDUR = 10;
STIMON = [.3 .7]; %stimulus onset varies.
MSKCLR = 0; %mask colour.
MSKON = 0; %mask onset
MSKDUR = 0; %mask duration

%margin (in cm) around each quadrant of the screen.
MRGN = 1;

%cm per inches
CMPERINCH=2.540000076;

%convert dpi to dots per cm.
dotsPerCm = dpi/CMPERINCH;

%acquire some useful info.
[FIELDNAMES,HEADERLINE,SAVEFORMATSPEC,~] = STConfigFileMakerAssist;
numFields = length(FIELDNAMES);

%determine the number of trials per block and the total number of trials.
trialTot = prod(factors)*ncp;

%check that the number of trials is evenly divisble by 'blocks'
if mod(trialTot,blocks)~=0
    disp(['Error! STODDBALLCONFIGFILEMAKER: check that trial total will be evenly '...
        'by ''blocks'''])
    return;
    
%otherwise, the total number of trials is divisible by the number of requested blocks.
%Determine the number of trials in each block.
else
    trialsPerBlock = trialTot/blocks;
end

%generate the factor trial order which provides info about the target image presentation
%on each trial
factorTrialOrder = CreateFactorialOrder(factors,ncp,repLim);

%generate a list of the desired image rotations (in degrees)
rotationsLst = linspace(0,360,factors(2)+1);
rotationsLst(end)=[];

%total number of rows in the config file, including the header row
rowsTotCnfgFl = trialTot*EVENTSPERTRIAL+1;

%determine the block number and trial number within the block list
blockNmLst = nan(rowsTotCnfgFl,1);
trialNmWthnBlckLst = nan(rowsTotCnfgFl,1);
for i = 1:blocks
    blockNmLst(1+(i-1)*trialsPerBlock:i*trialsPerBlock,1)=i*ones(trialsPerBlock,1);
    trialNmWthnBlckLst(1+(i-1)*trialsPerBlock:i*trialsPerBlock,1)=1:trialsPerBlock;
end

%acquire the images that will be used to populate the generated config file
[imgFNm, imgFDir] = uigetfile(IMGEXT,'Multiselect','On');

%image file name will be a cell if more than one image file was selected
if ~iscell(imgFNm)
    numImages=1;
    imgFNm = {imgFNm};
else  
    numImages = length(imgFNm);
end

%acquire the pixel size of each image. Determine the image display X and Y size in cm.
%Also determine the limits of the centre X and Y positions for each image. Note that the
%display is split up into 4 quadrants within which a single image can appear. These are
%stored in a rx2x4 (r: image number x min then max x 'p' quadrant (1: upper left, 2: upper
%right; 3: bottom left; 4: bottom right) 3D matrix.
minMaxImgCntrXLmtsCm=zeros(numImages,2,4);
maxMaxImgCntrYLmtsCm=zeros(numImages,2,4);

scrnXCm = scrnXPxls/dotsPerCm; %screen width in cm
scrnYCm = scrnYPxls/dotsPerCm; %screen height in cm

for i = 1:numImages
    [imgSzYPxls(i),imgSzXPxls(i),~] = size(imread([imgFDir imgFNm{i}]));

    %scale picture so its max dimension = maxImgSize
    scaleFactor = (maxDispImgSizeCm*dotsPerCm)/max([imgSzYPxls(i) imgSzXPxls(i)]);
    dispImgSzXCm(i) = (imgSzXPxls(i)*scaleFactor)/dotsPerCm;
    dispImgSzYCm(i) = (imgSzYPxls(i)*scaleFactor)/dotsPerCm;
   
    %use the display image diagonal to help determine the centre position limits.
    dispImgDgnlCm(i) = norm([dispImgSzXCm(i) dispImgSzYCm(i)]);
   
    %X min limits (each quadrant)
    minMaxImgCntrXLmtsCm(i,1,1) = MRGN+(dispImgDgnlCm(i)/2); %left
    minMaxImgCntrXLmtsCm(i,1,2) = (scrnXCm/2)+MRGN+(dispImgDgnlCm(i)/2); %right
    minMaxImgCntrXLmtsCm(i,1,3) = (scrnXCm/2)+MRGN+(dispImgDgnlCm(i)/2); %right
    minMaxImgCntrXLmtsCm(i,1,4) = MRGN+(dispImgDgnlCm(i)/2); %left
    
    %X max limits (each quadrant)
    minMaxImgCntrXLmtsCm(i,2,1) = (scrnXCm/2)-MRGN-(dispImgDgnlCm(i)/2); %left
    minMaxImgCntrXLmtsCm(i,2,2) = scrnXCm-MRGN-(dispImgDgnlCm(i)/2); %right
    minMaxImgCntrXLmtsCm(i,2,3) = scrnXCm-MRGN-(dispImgDgnlCm(i)/2); %right
    minMaxImgCntrXLmtsCm(i,2,4) = (scrnXCm/2)-MRGN-(dispImgDgnlCm(i)/2); %left
    
    %Y min limits (each quadrant)
    minMaxImgCntrYLmtsCm(i,1,1) = (scrnYCm/2)+MRGN+(dispImgDgnlCm(i)/2); %bottom
    minMaxImgCntrYLmtsCm(i,1,2) = (scrnYCm/2)+MRGN+(dispImgDgnlCm(i)/2); %bottom
    minMaxImgCntrYLmtsCm(i,1,3) = MRGN+(dispImgDgnlCm(i)/2); %top
    minMaxImgCntrYLmtsCm(i,1,4) = MRGN+(dispImgDgnlCm(i)/2); %top
    
    %Y max limits (each quadrant)
    minMaxImgCntrYLmtsCm(i,2,1) = scrnYCm-MRGN-(dispImgDgnlCm(i)/2); %bottom
    minMaxImgCntrYLmtsCm(i,2,2) = scrnYCm-MRGN-(dispImgDgnlCm(i)/2); %bottom
    minMaxImgCntrYLmtsCm(i,2,3) = (scrnYCm/2)-MRGN-(dispImgDgnlCm(i)/2); %top
    minMaxImgCntrYLmtsCm(i,2,4) = (scrnYCm/2)-MRGN-(dispImgDgnlCm(i)/2); %top
end

%we'll use FPRINT, so rows in 'A' will end up as columns (and cols in 'A' end up rows)
%when 'A' is saved to disk.
A = cell(numFields,rowsTotCnfgFl);

%column index for A
ACol = 2;

%populate 'A' for FPRINTF
for i = 1:trialTot
    
    %temporary list of the image names from which to select the distractor image
    tempImgLst = imgFNm;

    %temporary list of rotations 
    tempRtns = rotationsLst;
    
    %loop through each event in the trial
    for j=1:EVENTSPERTRIAL
  
        A{1,ACol} = uint32(blockNmLst(i));          %1: block number
        A{2,ACol} = uint32(0);                      %2: is_practice_block?
        A{3,ACol} = uint32(0);                      %3: block_has_feedback?
        A{4,ACol} = 0;                              %4: block_feedback_thresh
        A{5,ACol} = uint32(trialNmWthnBlckLst(i));  %5: trial_num within the block
        A{7,ACol} = uint32(0);                      %7: trial_feedback
        A{8,ACol} = BCKGRNDCLRNM;                   %8: background colour name
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   DETERMINE THE IMAGE NAME    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        %the image number
        imgNumIndx = factorTrialOrder(i,2);

        %the quadrant it will appear in
        qdrntNumImgIndx = factorTrialOrder(i,4);

        %the first image
        if j==1

            %9: stim_img_name
            A{9,ACol} = imgFNm{factorTrialOrder(i,2)}(1:end-4);                  

            %6: trial_max_time and 10: stim_onset
            if length(STIMON) == 2 %if range, select a random time within the range
                rndStimOn = (STIMON(2)-STIMON(1))*rand();
                A{6,ACol} = TRIALDUR-rndStimOn;
                A{10,ACol} = rndStimOn;
            else
                A{6,ACol} = TRIALDUR-STIMON; %otherwise, use a constant trial duration.
                A{10,ACol} = STIMON; %otherwise, use a constant onset.
            end

            %12: stim_cent_x (proportion of the screen width)
            A{12,ACol} = (minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumImgIndx)+...
                ((minMaxImgCntrXLmtsCm(imgNumIndx,2,qdrntNumImgIndx)-...
                minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumImgIndx))*rand()))/scrnXCm;

            %13: stim_cent_y (proportion of the screen height)
            A{13,ACol} = (minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumImgIndx)+...
                ((minMaxImgCntrYLmtsCm(imgNumIndx,2,qdrntNumImgIndx)-...
                minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumImgIndx))*rand()))/scrnYCm;

            %store the target image X and Y-size and rotation
            A{14,ACol} = dispImgSzXCm(imgNumIndx);                %14: stim_size_x (cm)
            A{15,ACol} = dispImgSzYCm(imgNumIndx);                %15: stim_size_y (cm)   
            A{16,ACol} = rotationsLst(factorTrialOrder(i,3));   %16: stim_rotation (degrees)
            A{18,ACol} = uint32(2);                             %18: stim_is_target

        %second image    
        elseif j==2
            
            %determine the quadrant for the second image
            switch qdrntNumImgIndx
                case 1 %bottom left
                    qdrntNumImgIndx = 2;
                case 2 %bottom right
                    qdrntNumImgIndx = 1;
                case 3 %top right
                    qdrntNumImgIndx = 4;
                case 4 %top left
                    qdrntNumImgIndx = 3;
            end 

            %'same' shape condition
            if factorTrialOrder(i,end)==1

                %9: stim_img_name
                A{9,ACol} = imgFNm{factorTrialOrder(i,2)}(1:end-4);

            %'different' shape condition
            elseif factorTrialOrder(i,end)==2

                %remove the target image name from the temporary file name list from which we
                %will draw a random distracter image name
                tempImgLst(imgNumIndx)=[];

                %update 'imgNumIndx' to reflect a different image
                imgNumIndx = randi(length(tempImgLst));

                %9: stim_img_name
                A{9,ACol} = tempImgLst{imgNumIndx}(1:end-4);      
            end

            %6: trial_max_time and 10: stim_onset
            if length(STIMON) == 2
                A{6,ACol} = TRIALDUR-rndStimOn;
                A{10,ACol} = rndStimOn;
            else
                A{6,ACol} = TRIALDUR-STIMON; %otherwise, use a constant trial duration.
                A{10,ACol} = STIMON; %otherwise, use a constant onset.
            end

            %12: stim_cent_x (proportion of the screen width)
            A{12,ACol} = (minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumImgIndx)+...
                ((minMaxImgCntrXLmtsCm(imgNumIndx,2,qdrntNumImgIndx)-...
                minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumImgIndx))*rand()))/scrnXCm;

            %13: stim_cent_y (proportion of the screen height)
            A{13,ACol} = (minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumImgIndx)+...
                ((minMaxImgCntrYLmtsCm(imgNumIndx,2,qdrntNumImgIndx)-...
                minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumImgIndx))*rand()))/scrnYCm;

            %store the target image X and Y-size and rotation
            A{14,ACol} = dispImgSzXCm(imgNumIndx);                %14: stim_size_x (cm)
            A{15,ACol} = dispImgSzYCm(imgNumIndx);                %15: stim_size_y (cm)

            %generate a random index to select a rotation at random from the temporary
            %list of rotations
            rtnIndx = randi(length(tempRtns));
            A{16,ACol} = rotationsLst(rtnIndx);                 %16: stim_rotation (degrees)
            A{18,ACol} = uint32(2);                             %18: stim_is_target
        end
            
        %11: stim duration and %23: mask_duration (seconds)      
        if length(STIMON) == 2 %if range, use the random stim onset.
            
            %if stim onset plus the stim duration last longer than the trial duration,
            %truncate the stimulus duration
            if (STIMDUR + rndStimOn) > TRIALDUR
                A{11,ACol} = TRIALDUR - rndStimOn;
                
                %if the mask duration is empty, make it equivalent to the stimulus
                %duration
                if isempty(MSKDUR)
                    A{23,ACol} = TRIALDUR - rndStimOn;
                else
                    A{23,ACol} = MSKDUR;
                end
            
            %otherwise, use 'STIMDUR' as-is
            else
                A{11,ACol} = STIMDUR;
                
                %if 'MSKDUR' is empty, use STIMDUR, otherwise, use 'MSKDUR'
                if isempty(MSKDUR)
                    A{23,ACol} = STIMDUR;
                else
                    A{23,ACol} = MSKDUR;
                end                  
            end
        
        %no range used. Thus, 'STIMON' is fixed.
        else
            
            %if stim onset plus the stim duration last longer than the trial duration,
            %truncate the stimulus duration
            if (STIMDUR + STIMON) > TRIALDUR
                A{11,ACol} = TRIALDUR - STIMON;
                
                %if the mask duration is empty, make it equivalent to the stimulus
                %duration
                if isempty(MSKDUR)
                    A{23,ACol} = TRIALDUR - rndStimOn;
                else
                    A{23,ACol} = MSKDUR;
                end
            
            %otherwise, use 'STIMDUR' as-is
            else
                A{11,ACol} = STIMDUR;
                
                %if the mask duration is empty, make it equivalent to the stimulus
                %duration
                if isempty(MSKDUR)
                    A{23,ACol} = TRIALDUR - rndStimOn;
                else
                    A{23,ACol} = MSKDUR;
                end
            end
        end

        A{17,ACol} = uint32(1);                         %17: stim_is_touchable
        A{19,ACol} = uint32(0);                         %19: subj_fixation_type
        A{20,ACol} = 0;                                 %20: subj_fixation_onset (seconds)
        A{21,ACol} = 0;                                 %21: subj_fixation_duration (seconds)
        
        %22: mask_onset (seconds)
        
        %if 'MSKON' is empty, make it coincident with the stimulus onset
        if isempty(MSKON)
            if length(STIMON) == 2 %use the same onset as the targets and distractors
                A{22,ACol} = rndStimOn;
            else
                A{22,ACol} = STIMON; %otherwise, use a constant onset.
            end
            
        %otherwise, 'MSKON' is not empty, use that value.    
        else
            A{22,ACol} = MSKON; %otherwise, use a constant onset.
        end      
                             
        A{24,ACol} = 0;                                 %24: mask_size
        A{25,ACol} = uint32(MSKCLR);                    %25: mask_color
        A{26,ACol} = 0;                                 %26: mask_rotation (degrees)
        A{27,ACol} = uint32(0);                         %27: mask_fit
        A{28,ACol} = 0;                                 %28: mask_margin
     
        %increment ACol
        ACol=ACol+1;
    end
end

%save the header line to file, store 'A' to file, then close the file.
if ~isempty(cnfFNm)
    fid = fopen([cnfFNm '.csv'],'w');
    fprintf(fid,HEADERLINE);
    fprintf(fid,SAVEFORMATSPEC,A{:,2:end});
    fclose(fid);
end
return;