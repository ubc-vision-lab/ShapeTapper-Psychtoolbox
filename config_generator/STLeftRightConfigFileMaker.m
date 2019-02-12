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
%           'numTrialOrders' - the number of configuration files to create.
%           'cnfFNm' - the configuration file name. If empty then the trial order is not
%               saved to disk.

function [A] = STLeftRightConfigFileMaker(factors,ncp,repLim,blocks,scrnXPxls,scrnYPxls,dpi,maxDispImgSizeCm,numTrialOrders,cnfFNm)

%the number of events (images to present) per trial
EVENTSPERTRIAL = 2;

%image file extension
IMGEXT = '*.png';

%the background colour
BCKGRNDCLRNM = '0 0 0';
TXTCOLOUR = '1 1 1';
%STIMDUR = single(10);
STIMDUR = single((1/30)*10);
TRIALDUR = single(10);
STIMON = single([.3 .7]); %stimulus onset varies.
MSKCLR = single(0); %mask colour.
MSKON = single(0); %mask onset
MSKDUR = single(0); %mask duration

%margin (in cm) around each quadrant of the screen.
MRGN = single(1);

%cm per inches
CMPERINCH=2.540000076;

%convert dpi to dots per cm.
dotsPerCm = dpi/CMPERINCH;

%acquire some useful info.
[FIELDNAMES,HEADERLINE,SAVEFORMATSPEC,~,CNDFIELDNAMES,CNDHEADERLINE,CNDSAVEFORMATSPEC,~] = STLeftRightConfigFileMakerAssist;
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

%acquire the images that will be used to populate the generated config file
[imgFNm, imgFDir] = uigetfile(IMGEXT,'Multiselect','On');

%image file name will be a cell if more than one image file was selected
if ~iscell(imgFNm)
    numImages=1;
    imgFNm = {imgFNm};
else  
    numImages = length(imgFNm);
end

%create each trial order...
for i=1:numTrialOrders
    
    %'B' stores the condition information to save to disk as a seperate file
    B=cell(13,trialTot);

    %overall trial number
    B(1,:)=num2cell(uint8(1:trialTot));
    
    %store the block #
    B(2,:) = num2cell(uint8(reshape(repmat(1:blocks,[trialsPerBlock 1]),[1 trialTot])));

    %store the trial # w/n block
    B(3,:)=num2cell(uint8(repmat(1:trialsPerBlock,[1 blocks])));

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
    for j = 1:blocks
        blockNmLst(1+(j-1)*trialsPerBlock:j*trialsPerBlock,1)=j*ones(trialsPerBlock,1);
        trialNmWthnBlckLst(1+(j-1)*trialsPerBlock:j*trialsPerBlock,1)=1:trialsPerBlock;
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

    for j = 1:numImages
        [imgSzYPxls(j),imgSzXPxls(j),~] = size(imread([imgFDir imgFNm{j}]));

        %scale picture so its max dimension = maxImgSize
        scaleFactor = (maxDispImgSizeCm*dotsPerCm)/max([imgSzYPxls(j) imgSzXPxls(j)]);
        dispImgSzXCm(j) = (imgSzXPxls(j)*scaleFactor)/dotsPerCm;
        dispImgSzYCm(j) = (imgSzYPxls(j)*scaleFactor)/dotsPerCm;

        %use the display image diagonal to help determine the centre position limits.
        dispImgDgnlCm(j) = norm([dispImgSzXCm(j) dispImgSzYCm(j)]);

        %X min limits (each quadrant)
        minMaxImgCntrXLmtsCm(j,1,1) = MRGN+(dispImgDgnlCm(j)/2); %left
        minMaxImgCntrXLmtsCm(j,1,2) = (scrnXCm/2)+MRGN+(dispImgDgnlCm(j)/2); %right
        minMaxImgCntrXLmtsCm(j,1,3) = (scrnXCm/2)+MRGN+(dispImgDgnlCm(j)/2); %right
        minMaxImgCntrXLmtsCm(j,1,4) = MRGN+(dispImgDgnlCm(j)/2); %left

        %X max limits (each quadrant)
        minMaxImgCntrXLmtsCm(j,2,1) = (scrnXCm/2)-MRGN-(dispImgDgnlCm(j)/2); %left
        minMaxImgCntrXLmtsCm(j,2,2) = scrnXCm-MRGN-(dispImgDgnlCm(j)/2); %right
        minMaxImgCntrXLmtsCm(j,2,3) = scrnXCm-MRGN-(dispImgDgnlCm(j)/2); %right
        minMaxImgCntrXLmtsCm(j,2,4) = (scrnXCm/2)-MRGN-(dispImgDgnlCm(j)/2); %left

        %Y min limits (each quadrant)
        minMaxImgCntrYLmtsCm(j,1,1) = (scrnYCm/2)+MRGN+(dispImgDgnlCm(j)/2); %bottom
        minMaxImgCntrYLmtsCm(j,1,2) = (scrnYCm/2)+MRGN+(dispImgDgnlCm(j)/2); %bottom
        minMaxImgCntrYLmtsCm(j,1,3) = MRGN+(dispImgDgnlCm(j)/2); %top
        minMaxImgCntrYLmtsCm(j,1,4) = MRGN+(dispImgDgnlCm(j)/2); %top

        %Y max limits (each quadrant)
        minMaxImgCntrYLmtsCm(j,2,1) = scrnYCm-MRGN-(dispImgDgnlCm(j)/2); %bottom
        minMaxImgCntrYLmtsCm(j,2,2) = scrnYCm-MRGN-(dispImgDgnlCm(j)/2); %bottom
        minMaxImgCntrYLmtsCm(j,2,3) = (scrnYCm/2)-MRGN-(dispImgDgnlCm(j)/2); %top
        minMaxImgCntrYLmtsCm(j,2,4) = (scrnYCm/2)-MRGN-(dispImgDgnlCm(j)/2); %top
    end

    %we'll use FPRINT, so rows in 'A' will end up as columns (and cols in 'A' end up rows)
    %when 'A' is saved to disk.
    A = cell(numFields,rowsTotCnfgFl);

    %column index for A
    ACol = 2;

    %populate 'A' for FPRINTF
    for j = 1:trialTot

        %temporary list of the image names from which to select the distractor image
        tempImgLst = imgFNm;

        %temporary list of rotations 
        tempRtns = rotationsLst;

        %loop through each event in the trial
        for k=1:EVENTSPERTRIAL

            A{1,ACol} = uint8(blockNmLst(j));           %1: block number
            A{2,ACol} = uint8(0);                       %2: is_practice_block?
            A{3,ACol} = uint8(0);                       %3: block_has_feedback?
            A{4,ACol} = single(0);                      %4: block_feedback_thresh
            A{5,ACol} = uint16(trialNmWthnBlckLst(j));  %5: trial_num within the block
            A{7,ACol} = uint8(0);                       %7: trial_timeout_msg
            A{8,ACol} = uint8(0);                       %8: trial_kb_resp?
            A{9,ACol} = 'None';                         %9: correct_kb_resp
            A{10,ACol} = uint8(0);                      %10: trial_feedback_type
            A{11,ACol} = 'None';                        %11: trial_feedback
            A{12,ACol} = BCKGRNDCLRNM;                  %12: background_colour RGB
            A{13,ACol} = TXTCOLOUR;                     %13: text_colour

            A{22,ACol} = uint8(1);                      %22: stim_is_touchable
            A{24,ACol} = uint8(0);                      %24: subj_fixation_type
            A{25,ACol} = uint8(0);                      %25: subj_fixation_pause?
            A{26,ACol} = single(0);                     %26: subj_fixation_onset (seconds)
            A{27,ACol} = single(0);                     %27: subj_fixation_duration (seconds)

            A{30,ACol} = single(0);                     %30: mask_size
            A{31,ACol} = MSKCLR;                        %31: mask_color
            A{32,ACol} = single(0);                     %32: mask_rotation (degrees)
            A{33,ACol} = uint8(0);                      %33: mask_fit
            A{34,ACol} = single(0);                     %34: mask_margin

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %   DETERMINE THE IMAGE NAME    %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %the image number
            imgNumIndx = factorTrialOrder(j,2);

            %the quadrant it will appear in
            qdrntNumImgIndx = factorTrialOrder(j,4);

            %the first image
            if k==1

                %14: stim_img_name
                A{14,ACol} = imgFNm{factorTrialOrder(j,2)}(1:end-4);
                
                %store the target img name
                B{4,j}=A{14,ACol};

                %6: trial_max_time and 15: stim_onset
                if length(STIMON) == 2 %if range, select a random time within the range
                    rndStimOn = (STIMON(2)-STIMON(1))*single(rand());
                    A{6,ACol} = TRIALDUR-rndStimOn;
                    A{15,ACol} = rndStimOn;
                else
                    A{6,ACol} = TRIALDUR-STIMON; %otherwise, use a constant trial duration.
                    A{15,ACol} = STIMON; %otherwise, use a constant onset.
                end

                %17: stim_cent_x (proportion of the screen width)
                A{17,ACol} = (minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumImgIndx)+...
                    ((minMaxImgCntrXLmtsCm(imgNumIndx,2,qdrntNumImgIndx)-...
                    minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumImgIndx))*rand()))/scrnXCm;

                %store the target img x-position
                B{5,j}=A{17,ACol};
                
                %18: stim_cent_y (proportion of the screen height)
                A{18,ACol} = (minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumImgIndx)+...
                    ((minMaxImgCntrYLmtsCm(imgNumIndx,2,qdrntNumImgIndx)-...
                    minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumImgIndx))*rand()))/scrnYCm;
                
                %store the target img y-position
                B{6,j}=A{18,ACol};
                
                %store the target image X and Y-size and rotation
                A{19,ACol} = single(dispImgSzXCm(imgNumIndx));              %19: stim_size_x (cm)
                A{20,ACol} = single(dispImgSzYCm(imgNumIndx));              %20: stim_size_y (cm)   
                A{21,ACol} = single(rotationsLst(factorTrialOrder(j,3)));   %21: stim_rotation (degrees)
                A{23,ACol} = uint8(1);                                      %23: stim_is_target
                
                %store the target rotation
                B{7,j}=A{21,ACol};

            %second image    
            elseif k==2

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
                if factorTrialOrder(j,end)==1
                    
                    %14: stim_img_name
                    A{14,ACol} = imgFNm{factorTrialOrder(j,2)}(1:end-4);
                    
                    %store the distractor img1 name
                    B{8,j}=A{14,ACol};

                %'different' shape condition
                elseif factorTrialOrder(j,end)==2

                    %remove the target image name from the temporary file name list from which we
                    %will draw a random distracter image name
                    tempImgLst(imgNumIndx)=[];

                    %update 'imgNumIndx' to reflect a different image
                    imgNumIndx = randi(length(tempImgLst));

                    %14: stim_img_name
                    A{14,ACol} = tempImgLst{imgNumIndx}(1:end-4); 
                    
                    %store the distractor img1 name
                    B{8,j}=A{14,ACol};
                end

                %6: trial_max_time and 15: stim_onset
                if length(STIMON) == 2
                    A{6,ACol} = TRIALDUR-rndStimOn;
                    A{15,ACol} = rndStimOn;
                else
                    A{6,ACol} = TRIALDUR-STIMON; %otherwise, use a constant trial duration.
                    A{15,ACol} = STIMON; %otherwise, use a constant onset.
                end

                %17: stim_cent_x (proportion of the screen width)
                A{17,ACol} = (minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumImgIndx)+...
                    ((minMaxImgCntrXLmtsCm(imgNumIndx,2,qdrntNumImgIndx)-...
                    minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumImgIndx))*rand()))/scrnXCm;
                
                %store the distractor img1 x-position
                B{9,j}=A{17,ACol};

                %18: stim_cent_y (proportion of the screen height)
                A{18,ACol} = (minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumImgIndx)+...
                    ((minMaxImgCntrYLmtsCm(imgNumIndx,2,qdrntNumImgIndx)-...
                    minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumImgIndx))*rand()))/scrnYCm;
                
                %store the distractor img1 y-position
                B{10,j}=A{18,ACol};

                %store the target image X and Y-size and rotation
                A{19,ACol} = dispImgSzXCm(imgNumIndx);                %19: stim_size_x (cm)
                A{20,ACol} = dispImgSzYCm(imgNumIndx);                %20: stim_size_y (cm)

                %generate a random index to select a rotation at random from the temporary
                %list of rotations
                rtnIndx = randi(length(tempRtns));
                A{21,ACol} = rotationsLst(rtnIndx);                 %21: stim_rotation (degrees)
                A{23,ACol} = uint8(1);                              %23: stim_is_target
                
                %store the distractor img1 rotation
                B{11,j}=A{21,ACol};
            end

            %16: stim duration and %29: mask_duration (seconds)      
            if length(STIMON) == 2 %if range, use the random stim onset.

                %if stim onset plus the stim duration last longer than the trial duration,
                %truncate the stimulus duration
                if (STIMDUR + rndStimOn) > TRIALDUR
                    A{16,ACol} = TRIALDUR - rndStimOn;

                    %if the mask duration is empty, make it equivalent to the stimulus
                    %duration
                    if isempty(MSKDUR)
                        A{29,ACol} = TRIALDUR - rndStimOn;
                    else
                        A{29,ACol} = MSKDUR;
                    end

                %otherwise, use 'STIMDUR' as-is
                else
                    A{16,ACol} = STIMDUR;

                    %if 'MSKDUR' is empty, use STIMDUR, otherwise, use 'MSKDUR'
                    if isempty(MSKDUR)
                        A{29,ACol} = STIMDUR;
                    else
                        A{29,ACol} = MSKDUR;
                    end                  
                end

            %no range used. Thus, 'STIMON' is fixed.
            else

                %if stim onset plus the stim duration last longer than the trial duration,
                %truncate the stimulus duration
                if (STIMDUR + STIMON) > TRIALDUR
                    A{16,ACol} = TRIALDUR - STIMON;

                    %if the mask duration is empty, make it equivalent to the stimulus
                    %duration
                    if isempty(MSKDUR)
                        A{29,ACol} = TRIALDUR - rndStimOn;
                    else
                        A{29,ACol} = MSKDUR;
                    end

                %otherwise, use 'STIMDUR' as-is
                else
                    A{16,ACol} = STIMDUR;

                    %if the mask duration is empty, make it equivalent to the stimulus
                    %duration
                    if isempty(MSKDUR)
                        A{29,ACol} = TRIALDUR - rndStimOn;
                    else
                        A{29,ACol} = MSKDUR;
                    end
                end
            end

            %28: mask_onset (seconds)

            %if 'MSKON' is empty, make it coincident with the stimulus onset
            if isempty(MSKON)
                if length(STIMON) == 2 %use the same onset as the targets and distractors
                    A{28,ACol} = rndStimOn;
                else
                    A{28,ACol} = STIMON; %otherwise, use a constant onset.
                end

            %otherwise, 'MSKON' is not empty, use that value.    
            else
                A{28,ACol} = MSKON; %otherwise, use a constant onset.
            end
            
            %store the img onset
            B{12,j}=A{15,ACol};
                
            %store the img duration
            B{13,j}=A{16,ACol};

            %increment ACol
            ACol=ACol+1;
        end
    end

    %save the header line to file, store 'A' to file, then close the file.
    if ~isempty(cnfFNm)
        
        %add the trial order number to the file name
        sffxFNm=[];
        if i<10
            sffxFNm=[sffxFNm '_0' num2str(i)];
        elseif i>9 && i<100
            sffxFNm=[sffxFNm num2str(i)];
        end
          
        %save the configuration file to disk    
        fid = fopen([cnfFNm sffxFNm '.csv'],'w');
        fprintf(fid,HEADERLINE);
        fprintf(fid,SAVEFORMATSPEC,A{:,2:end});
        fclose(fid);
        
        %save a csv file of the condition order and other relevant information.
        fid = fopen([cnfFNm sffxFNm '_cond_file.csv'],'w');
        fprintf(fid,CNDHEADERLINE);
        fprintf(fid,CNDSAVEFORMATSPEC,B{:,:});
        fclose(fid);
        B=[CNDFIELDNAMES B];
        B=transpose(B);
    end

%end of the trial order loop
end
return;