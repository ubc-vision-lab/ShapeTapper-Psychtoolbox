%HELP for function STODDBALLCONFIGMAKER. This function generates a 3 image odd-ball
%comma-separated configuration file for the Matlab-based ShapeTapper program and stores
%it to disk.
%
%@PARAMS
%           'factors' - a 3-element array in which each element is a factor the value of which
%               specifies the number of levels in that given factor. Assumes the following
%               structure: (1) number of shapes; (2) number of unique orientations; (3)
%               the quadrant in which the target (oddball) shape should occur.
%           'tpc' - the number of trials per unique condition in the trial order
%           'repLim' - the maximum number of back-to-back repetitions of any given unique
%               condition that is permitted to occur in the condition order.
%           'blocks' - the number of blocks of trials in the trial order
%           'scrnXPxls - the width of the display screen (in pixels)
%           'scrnYPxls - the height of the display screen (in pixels)
%           'dpi' - the dots per inch of the display screen
%           'maxDispImgSizeCm' - the maximum size of the image on the display screen (cm)
%           'numTrialOrders' - the number of configuration files to create.
%           'cnfFNm' - the configuration file name. If empty then the trial order is not
%               saved to disk.
%
%e.g., >>STOddBallConfigFileMaker([6 10 4],1,4,5,1600,900,96,8,'ob_long_5_48tpb_6s');

function [A,B] = STOddBallConfigFileMaker(factors,tpc,repLim,blocks,scrnXPxls,scrnYPxls,dpi,maxDispImgSizeCm,numTrialOrders,cnfFNm)

%the number of events (images to present) per trial
EVENTSPERTRIAL = 3;

%image file extension
IMGEXT = '*.png';

%the background colour
BCKGRNDCLRNM = '0 0 0';
TXTCOLOUR = '1 1 1';
STIMDUR = single(10);
%STIMDUR = single((1/30)*10); 
TRIALDUR = single(10);
STIMON = single([.3 .7]); %stimulus onset varies.
MSKCLR = single(0); %mask colour.
MSKON = single(0); %mask onset.
MSKDUR = single(0); %mask duration.

%margin (in cm) around each quadrant of the screen.
MRGN = single(1);

%cm per inches
CMPERINCH=2.540000076;

%convert dpi to dots per cm.
dotsPerCm = dpi/CMPERINCH;

%acquire some useful info.
[FIELDNAMES,HEADERLINE,SAVEFORMATSPEC,~,CNDFIELDNAMES,CNDHEADERLINE,CNDSAVEFORMATSPEC,~] = STOddBallConfigFileMakerAssist;
numFields = length(FIELDNAMES);

%determine the number of trials per block and the total number of trials.
trialTot = prod(factors)*tpc;

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

%acquire the pixel size of each image. Determine the image display X and Y size in cm.
%Also determine the limits of the centre X and Y positions for each image. Note that the
%display is split up into 4 quadrants within which a single image can appear. These are
%stored in a rx2x4 (r: image number x min then max x 'p' quadrant (1: bottom left, 2:
%bottom right; 3: top right; 4: top left) 3D matrix.
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

%for each desired trial order...
for i=1:numTrialOrders
    
    %'B' stores the condition information to save to disk as a seperate file
    B=cell(17,trialTot);

    %overall trial number
    B(1,:)=num2cell(uint8(1:trialTot));
    
    %store the block #
    B(2,:) = num2cell(uint8(reshape(repmat(1:blocks,[trialsPerBlock 1]),[1 trialTot])));

    %store the trial # w/n block
    B(3,:)=num2cell(uint8(repmat(1:trialsPerBlock,[1 blocks])));

    %generate the factor trial order which provides info about the target image presentation
    %on each trial
    factorTrialOrder = CreateFactorialOrder(factors,tpc,repLim);

    %generate a list of the desired image rotations (in degrees)
    rotationsLst = single(linspace(0,360,factors(2)+1));
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

    %we'll use FPRINT, so rows in 'A' will end up as columns (and cols in 'A' end up rows)
    %when 'A' is saved to disk.
    A = cell(numFields,rowsTotCnfgFl);

    %column index for A
    ACol = 2;

    %populate 'A' for FPRINTF
    for j = 1:trialTot

        %temporary list of the image names from which to select the distractor image
        tempImgLst = imgFNm;

        %temporary list of quadrants the image can be positioned at
        tempQudrntLst=1:4;

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
            A{11,ACol} = 'None';                        %11: trial_feedback_img
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

            %target image
            if k==1

                %14: stim_img_name
                A{14,ACol} = imgFNm{factorTrialOrder(j,2)}(1:end-4);                  
                
                %store the target img name
                B{4,j}=A{14,ACol};
                
                %the image number
                imgNumIndx = factorTrialOrder(j,2);
                qdrntNumIndx = factorTrialOrder(j,4);

                %6: trial_max_time and 15: stim_onset
                if length(STIMON) == 2 %if range, select a random time within the range
                    rndStimOn = single((STIMON(2)-STIMON(1))*rand());
                    A{6,ACol} = TRIALDUR-rndStimOn;
                    A{15,ACol} = rndStimOn;
                else
                    A{6,ACol} = TRIALDUR-STIMON; %otherwise, use a constant trial duration.
                    A{15,ACol} = STIMON; %otherwise, use a constant onset.
                end
                
                %17: stim_cent_x (proportion of the screen width)
                A{17,ACol} = (minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumIndx)+...
                    ((minMaxImgCntrXLmtsCm(imgNumIndx,2,qdrntNumIndx)-...
                    minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumIndx))*rand()))/scrnXCm;
                
                %store the target img x-position
                B{5,j}=A{17,ACol};

                %18: stim_cent_y (proportion of the screen height)
                A{18,ACol} = (minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumIndx)+...
                    ((minMaxImgCntrYLmtsCm(imgNumIndx,2,qdrntNumIndx)-...
                    minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumIndx))*rand()))/scrnYCm;
                
                %store the target img y-position
                B{6,j}=A{18,ACol};

                %store the target image X and Y-size and rotation
                A{19,ACol} = dispImgSzXCm(imgNumIndx);             %19: stim_size_x (cm)
                A{20,ACol} = dispImgSzYCm(imgNumIndx);             %20: stim_size_y (cm)   
                A{21,ACol} = rotationsLst(factorTrialOrder(j,3));  %21: stim_rotation (degrees)
                A{23,ACol} = uint8(1);                             %23: stim_is_target
                
                %store the target rotation
                B{7,j}=A{21,ACol};
                
                %remove the target image name from the temporary file name list that we
                %use to select a random distracter image name
                tempImgLst(imgNumIndx)=[];
                rndmImgNmIndx = randi(length(tempImgLst));

                %remove the target image quadrant position from the temporary list of
                %quadrants from which we will draw a random quadrant from
                tempQudrntLst(qdrntNumIndx)=[];
                rndmQdrntNumIndx = randi(3);

                %generate a random index to select a rotation at random from the temporary
                %list of rotations
                rdnmRtnIndx = randi(length(tempRtns));

            %distractor image 1
            elseif k==2

                %select a random image file from the list of remaining images 
                A{14,ACol} = tempImgLst{rndmImgNmIndx}(1:end-4);
                
                %store the distractor img1 name
                B{8,j}=A{14,ACol};

                %6: trial_max_time and 15: stim_onset
                if length(STIMON) == 2 %if range, select a random time within the range
                    A{6,ACol} = TRIALDUR-rndStimOn;
                    A{15,ACol} = rndStimOn;
                else
                    A{6,ACol} = TRIALDUR-STIMON; %otherwise, use a constant trial duration.
                    A{15,ACol} = STIMON; %otherwise, use a constant onset.
                end

                %17: stim_cent_x (proportion of the screen width)
                A{17,ACol} = (minMaxImgCntrXLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx))+...
                    ((minMaxImgCntrXLmtsCm(rndmImgNmIndx,2,tempQudrntLst(rndmQdrntNumIndx))-...
                    minMaxImgCntrXLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx)))*rand()))/scrnXCm;
                
                %store the distractor img1 x-position
                B{9,j}=A{17,ACol};

                %18: stim_cent_y (proportion of the screen height)
                A{18,ACol} = (minMaxImgCntrYLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx))+...
                    ((minMaxImgCntrYLmtsCm(rndmImgNmIndx,2,tempQudrntLst(rndmQdrntNumIndx))-...
                    minMaxImgCntrYLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx)))*rand()))/scrnYCm;
                
                %store the distractor img1 y-position
                B{10,j}=A{18,ACol};

                %store the computed display size for the distractor image and its rotation
                A{19,ACol} = dispImgSzXCm(rndmImgNmIndx);           %19: stim_size_x (cm)
                A{20,ACol} = dispImgSzYCm(rndmImgNmIndx);           %20: stim_size_y (cm)
                A{21,ACol} = tempRtns(rdnmRtnIndx);                 %21: stim_rotation  
                A{23,ACol} = uint8(0);                              %23: stim_is_target
                
                %store the distractor img1 rotation
                B{11,j}=A{21,ACol};

                %remove the target image quadrant position from the temporary list of
                %quadrants from which we will draw a random quadrant from
                tempQudrntLst(rndmQdrntNumIndx)=[];
                rndmQdrntNumIndx = randi(2);

                %remove the used rotation from the list of rotations and generate a new index
                %at random from those that remain
                tempRtns(rdnmRtnIndx) = [];
                rdnmRtnIndx = randi(length(tempRtns));

            %distractor image 2
            elseif k==3

                %select the same random image file from the list of remaining images 
                A{14,ACol} = tempImgLst{rndmImgNmIndx}(1:end-4);
                
                %store the distractor img2 name
                B{12,j}=A{14,ACol};
                
                %if range, select a random time within the range.
                if length(STIMON) == 2 
                    A{6,ACol} = TRIALDUR-rndStimOn;             %6: trial_max_time
                    A{15,ACol} = rndStimOn;                     %15: stim_onset
                
                %otherwise, use the constant trial duration and onset
                else
                    A{6,ACol} = TRIALDUR-STIMON;                %6: trial_max_time
                    A{15,ACol} = STIMON;                        %15: stim_onset
                end

                %17: stim_cent_x (proportion of the screen width)
                A{17,ACol} = (minMaxImgCntrXLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx))+...
                    ((minMaxImgCntrXLmtsCm(rndmImgNmIndx,2,tempQudrntLst(rndmQdrntNumIndx))-...
                    minMaxImgCntrXLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx)))*rand()))/scrnXCm;
                
                %store the distractor img2 x-position
                B{13,j}=A{17,ACol};

                %18: stim_cent_y (proportion of the screen height)
                A{18,ACol} = (minMaxImgCntrYLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx))+...
                    ((minMaxImgCntrYLmtsCm(rndmImgNmIndx,2,tempQudrntLst(rndmQdrntNumIndx))-...
                    minMaxImgCntrYLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx)))*rand()))/scrnYCm;
                
                %store the distractor img2 y-position
                B{14,j}=A{18,ACol};

                %store the computed image size for the distractors
                A{19,ACol} = dispImgSzXCm(rndmImgNmIndx);           %19: stim_size_x (cm)
                A{20,ACol} = dispImgSzYCm(rndmImgNmIndx);           %20: stim_size_y (cm)
                A{21,ACol} = tempRtns(rdnmRtnIndx);                 %21: stim_rotation
                A{23,ACol} = uint8(0);                              %23: stim_is_target
                
                %store the distractor img2 rotation
                B{15,j}=A{21,ACol};
            end

            %16: stim duration and %28: mask_duration (seconds)      
            if length(STIMON) == 2 %if range, use the random stim onset.

                %if stim onset plus the stim duration last longer than the trial duration,
                %truncate the stimulus duration
                if (STIMDUR + rndStimOn) > TRIALDUR
                    A{16,ACol} = TRIALDUR - rndStimOn;          %16: stim_duration

                    %if the mask duration is empty, make it equivalent to the stimulus
                    %duration
                    if isempty(MSKDUR)
                        A{29,ACol} = TRIALDUR - rndStimOn;      %29: mask_duration
                    else
                        A{29,ACol} = MSKDUR;                    %29: mask_duration
                    end

                %otherwise, use 'STIMDUR' as-is
                else
                    A{16,ACol} = STIMDUR;                       %16: stim_duration

                    %if 'MSKDUR' is empty, use STIMDUR, otherwise, use 'MSKDUR'
                    if isempty(MSKDUR)
                        A{29,ACol} = STIMDUR;                   %29: mask_duration
                    else
                        A{29,ACol} = MSKDUR;                    %29: mask_duration
                    end                  
                end

            %no range used. Thus, 'STIMON' is fixed.
            else

                %if stim onset plus the stim duration last longer than the trial duration,
                %truncate the stimulus duration
                if (STIMDUR + STIMON) > TRIALDUR
                    A{16,ACol} = TRIALDUR - STIMON;             %16: stim_duration

                    %if the mask duration is empty, make it equivalent to the stimulus
                    %duration
                    if isempty(MSKDUR)
                        A{29,ACol} = TRIALDUR - rndStimOn;      %29: mask_duration
                    else
                        A{29,ACol} = MSKDUR;                    %29: mask_duration
                    end

                %otherwise, use 'STIMDUR' as-is
                else
                    A{16,ACol} = STIMDUR;                       %16: stim_duration

                    %if the mask duration is empty, make it equivalent to the stimulus
                    %duration
                    if isempty(MSKDUR)
                        A{29,ACol} = TRIALDUR - rndStimOn;      %29: mask_duration
                    else
                        A{29,ACol} = MSKDUR;                    %29: mask_duration
                    end
                end
            end

            %if 'MSKON' is empty, make it coincident with the stimulus onset
            if isempty(MSKON)
                if length(STIMON) == 2 %use the same onset as the targets and distractors
                    A{28,ACol} = rndStimOn;                     %28: mask_onset
                else
                    A{28,ACol} = STIMON;                        %28: mask_onset
                end

            %otherwise, use 'MSKON'  
            else
                A{28,ACol} = MSKON;                             %28: mask_onset
            end
            
            %store the img onset
            B{16,j}=A{15,ACol};
                
            %store the img duration
            B{17,j}=A{16,ACol};

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
end
return;