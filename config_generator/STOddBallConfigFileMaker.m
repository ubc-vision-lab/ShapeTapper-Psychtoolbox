function [A] = STOddBallConfigFileMaker(factors,n,r,blocks,scrnXPxls,scrnYPxls,dpi,maxImgSizeCm,cnfFNm)

%the number of events (images to present) per trial
EVENTSPERTRIAL = 3;

%the background colour
BCKGRNDCLRNM = 'black';
STIMDUR = 10;
TRIALDUR = 10;

%margin (in cm) around each quadrant of the screen.
MRGN = 1;

%cm per inches
CMPERINCH=2.540000076;

%convert dpi to dots per cm.
dotsPerCm = dpi/CMPERINCH;

%acquire some useful info.
[FIELDNAMES,HEADERLINE,SAVEFORMATSPEC,~] = STConfigFileMakerAssist;
numFields = length(FIELDNAMES);

%generate the factor trial order which provides info about the target image presentation
%on each trial
factorTrialOrder = CreateTrialOrder(factors,n,r,[],[],[]);

%generate a list of the desired image rotations (in degrees)
rotationsLst = linspace(0,360,factors(2)+1);
rotationsLst(end)=[];

%total number of trials (NOTE: this value will differ from the number of rows in the
%config file if there is more than one event per trial)
[trialTot, ~] = size(factorTrialOrder);

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
[imgFNm, imgFDir] = uigetfile('*.png','Multiselect','On');

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
    scaleFactor = (maxImgSizeCm*dotsPerCm)/max([imgSzYPxls(i) imgSzXPxls(i)]);
    dispImgSzXCm(i) = (imgSzXPxls(i)*scaleFactor)/dotsPerCm;
    dispImgSzYCm(i) = (imgSzYPxls(i)*scaleFactor)/dotsPerCm;
   
    %use the display image diagonal to help determine the centre position limits.
    dispImgDgnlCm(i) = norm([dispImgSzXCm(i) dispImgSzYCm(i)]);
   
    %X min limits (each quadrant)
    minMaxImgCntrXLmtsCm(i,1,1) = MRGN+(dispImgDgnlCm(i)/2);
    minMaxImgCntrXLmtsCm(i,1,2) = (scrnXCm/2)+MRGN+(dispImgDgnlCm(i)/2);
    minMaxImgCntrXLmtsCm(i,1,3) = (scrnXCm/2)+MRGN+(dispImgDgnlCm(i)/2);
    minMaxImgCntrXLmtsCm(i,1,4) = MRGN+(dispImgDgnlCm(i)/2);
    
    %X max limits (each quadrant)
    minMaxImgCntrXLmtsCm(i,2,1) = (scrnXCm/2)-MRGN-(dispImgDgnlCm(i)/2);
    minMaxImgCntrXLmtsCm(i,2,2) = scrnXCm-MRGN-(dispImgDgnlCm(i)/2);
    minMaxImgCntrXLmtsCm(i,2,3) = scrnXCm-MRGN-(dispImgDgnlCm(i)/2);
    minMaxImgCntrXLmtsCm(i,2,4) = (scrnXCm/2)-MRGN-(dispImgDgnlCm(i)/2);
    
    %Y min limits (each quadrant)
    minMaxImgCntrYLmtsCm(i,1,1) = (scrnYCm/2)+MRGN+(dispImgDgnlCm(i)/2);
    minMaxImgCntrYLmtsCm(i,1,2) = (scrnYCm/2)+MRGN+(dispImgDgnlCm(i)/2);
    minMaxImgCntrYLmtsCm(i,1,3) = MRGN+(dispImgDgnlCm(i)/2);
    minMaxImgCntrYLmtsCm(i,1,4) = MRGN+(dispImgDgnlCm(i)/2);
    
    %Y max limits (each quadrant)
    minMaxImgCntrYLmtsCm(i,2,1) = scrnYCm-MRGN-(dispImgDgnlCm(i)/2);
    minMaxImgCntrYLmtsCm(i,2,2) = scrnYCm-MRGN-(dispImgDgnlCm(i)/2);
    minMaxImgCntrYLmtsCm(i,2,3) = (scrnYCm/2)-MRGN-(dispImgDgnlCm(i)/2);
    minMaxImgCntrYLmtsCm(i,2,4) = (scrnYCm/2)-MRGN-(dispImgDgnlCm(i)/2); 
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
 
    %temporary list of quadrants the image can be positioned at
    tempQudrntLst=1:4;

    %temporary list of rotations
    tempRtns = rotationsLst;
    
    %loop through each event in the trial
    for j=1:EVENTSPERTRIAL
  
        A{1,ACol} = uint32(blockNmLst(i));          %1: block number
        A{2,ACol} = uint32(0);                      %2: is_practice_block?
        A{3,ACol} = uint32(0);                      %3: block_has_feedback?
        A{4,ACol} = 0;                              %4: block_feedback_thresh
        A{5,ACol} = uint32(trialNmWthnBlckLst(i));  %5: trial_num within the block
        A{6,ACol} = TRIALDUR;                       %6: trial_max_time
        A{7,ACol} = uint32(0);                      %7: trial_feedback
        A{8,ACol} = BCKGRNDCLRNM;                   %8: background colour name
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   DETERMINE THE IMAGE NAME    %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
        %target image
        if j==1
            
            %9: stim_img_name
            A{9,ACol} = imgFNm{factorTrialOrder(i,2)}(1:end-4);                  
       
            %the image number
            imgNumIndx = factorTrialOrder(i,2);
            qdrntNumIndx = factorTrialOrder(i,4);
            
            %12: stim_cent_x (proportion of the screen width)
            A{12,ACol} = (minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumIndx)+...
                ((minMaxImgCntrXLmtsCm(imgNumIndx,2,qdrntNumIndx)-...
                minMaxImgCntrXLmtsCm(imgNumIndx,1,qdrntNumIndx))*rand()))/scrnXCm;
            
            %13: stim_cent_y (proportion of the screen height)
            A{13,ACol} = (minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumIndx)+...
                ((minMaxImgCntrYLmtsCm(imgNumIndx,2,qdrntNumIndx)-...
                minMaxImgCntrYLmtsCm(imgNumIndx,1,qdrntNumIndx))*rand()))/scrnYCm;
            
            %store the target image X and Y-size and rotation
            A{14,ACol} = dispImgSzXCm(imgNumIndx);                %14: stim_size_x (cm)
            A{15,ACol} = dispImgSzYCm(imgNumIndx);                %15: stim_size_y (cm)   
            A{16,ACol} = rotationsLst(factorTrialOrder(i,3));   %16: stim_rotation (degrees)
            A{18,ACol} = uint32(2);                             %18: stim_is_target
            
            %remove the target image name from the temporary file name list from which we
            %will draw a random distracter image name
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
        elseif j==2

            %select the a random image file from the list of remaining images 
            A{9,ACol} = tempImgLst{rndmImgNmIndx}(1:end-4);
                   
            %12: stim_cent_x (proportion of the screen width)
            A{12,ACol} = (minMaxImgCntrXLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx))+...
                ((minMaxImgCntrXLmtsCm(rndmImgNmIndx,2,tempQudrntLst(rndmQdrntNumIndx))-...
                minMaxImgCntrXLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx)))*rand()))/scrnXCm;
            
            %13: stim_cent_y (proportion of the screen height)
            A{13,ACol} = (minMaxImgCntrYLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx))+...
                ((minMaxImgCntrYLmtsCm(rndmImgNmIndx,2,tempQudrntLst(rndmQdrntNumIndx))-...
                minMaxImgCntrYLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx)))*rand()))/scrnYCm;
            
            %store the computed display size for the distractor image and its rotation
            A{14,ACol} = dispImgSzXCm(rndmImgNmIndx);             %14: stim_size_x (cm)
            A{15,ACol} = dispImgSzYCm(rndmImgNmIndx);             %15: stim_size_y (cm)
            A{16,ACol} = tempRtns(rdnmRtnIndx);    
            A{18,ACol} = uint32(0);                             %18: stim_is_target
            
            %remove the target image quadrant position from the temporary list of
            %quadrants from which we will draw a random quadrant from
            tempQudrntLst(rndmQdrntNumIndx)=[];
            rndmQdrntNumIndx = randi(2);
            
            %remove the used rotation from the list of rotations and generate a new index
            %at random from those that remain
            tempRtns(rdnmRtnIndx) = [];
            rdnmRtnIndx = randi(length(tempRtns));
        
        %distractor image 2
        elseif j==3
            
            %select the same random image file from the list of remaining images 
            A{9,ACol} = tempImgLst{rndmImgNmIndx}(1:end-4);
            
            %12: stim_cent_x (proportion of the screen width)
            A{12,ACol} = (minMaxImgCntrXLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx))+...
                ((minMaxImgCntrXLmtsCm(rndmImgNmIndx,2,tempQudrntLst(rndmQdrntNumIndx))-...
                minMaxImgCntrXLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx)))*rand()))/scrnXCm;
            
            %13: stim_cent_y (proportion of the screen height)
            A{13,ACol} = (minMaxImgCntrYLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx))+...
                ((minMaxImgCntrYLmtsCm(rndmImgNmIndx,2,tempQudrntLst(rndmQdrntNumIndx))-...
                minMaxImgCntrYLmtsCm(rndmImgNmIndx,1,tempQudrntLst(rndmQdrntNumIndx)))*rand()))/scrnYCm;
            
            %store the computed image size for the distractors
            A{14,ACol} = dispImgSzXCm(rndmImgNmIndx);             %14: stim_size_x (cm)
            A{15,ACol} = dispImgSzYCm(rndmImgNmIndx);             %15: stim_size_y (cm)
            A{16,ACol} = tempRtns(rdnmRtnIndx);
            A{18,ACol} = uint32(0);                             %17: stim_is_target
        end       
            
        A{10,ACol} = 0;                                 %10: stim_onset
        A{11,ACol} = STIMDUR;                           %11: stim_duration   
        A{17,ACol} = uint32(1);                         %17: stim_is_touchable
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

%save the header line to file, store 'A' to file, then close the file.
if ~isempty(cnfFNm)
    fid = fopen([cnfFNm '.csv'],'w');
    fprintf(fid,HEADERLINE);
    fprintf(fid,SAVEFORMATSPEC,A{:,2:end});
    fclose(fid);
end
return;

