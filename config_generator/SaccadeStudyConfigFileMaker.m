%HELP for function SPIDERSTUDYCONFIGFILEMAKER.
%
%Author: Rob Whitwell
%
%This function generates trial orders for the shape spider study. It will prompt the user
%to select a 'top directory' within which subfolders for each image type (negative,
%neutral negative, neutral positive, and positive) will be sampled from randomly for each
%trial order. The subfolder names are "Negative", "Neutral_Negative", "Neutral Positive",
%and "Positive" but can be changed below in the relevant fixed variable(s). This function
%assumes 5 core conditions (3 'no-jump' positions and 2 'jump' conditions in which the
%image starts in the same middle positiona and jumps either left or right).
%
%@PARAMS:
%           'tpc' - the number of trials per unique condition
%           'repLim' - the maximum number of back-to-back repetitions of any given unique
%               condition that is permitted to occur in the condition order.
%           'blocks' - the number of blocks of trials in the trial order
%           'dualTask' - whether or not the config file is a dual task file. The dual task
%               is to report the letter in the alphanumeric stream at the end of the
%               trial.
%           'numTrialOrders' - the number of configuration files to create.
%           'cnfFNm' - the configuration file name. If empty then the trial order is not
%               saved to disk.
%
%@OUTPUT:
%           'A' - A cell matrix denoting the configuration file. Events organized by rows.
%               First row is the header-line.
%           'B' - A cell matrix denoting a summary of the conditions of 'A'.
%               First row is the header-line.
%
%e.g., create two dual task trial orders with 3 blocks of trials and use 'dual' in the
%       file name:
%           >> SpiderStudyConfigFileMaker(12,3,3,1,2,'dual');

function A = SaccadeStudyConfigFileMaker(tpc,repLim,blocks,numTrialOrders,cnfFNm)

%cm per inches
CMPERINCH=2.540000076;
SCREEN_WIDTH_INCHES = 37; %47.5;
SCREEN_PPI = 52.45; %80.69;
VIEWER_DIST_INCHES = 20; %35;

%FIXED VARIABLES
FIXIMGDISPSZXCM=single(2); %the desired width of the target image
FIXIMGDISPSZYCM=single(2); %the desired height of the target image
TARGETIMGDISPSZXCM=single(3); %the desired width of the target image
TARGETIMGDISPSZYCM=single(3); %the desired height of the target image
% ALPHAS='knvxyz'; %the target characters
% POSIMGDIR='Positive';
% NEGIMGDIR='Negative';
% NEUTPOSIMGDIR='Neutral_Positive';
% NEUTNEGIMGDIR='Neutral_Negative';

%display info: 32" Elo (698.4mm x 392.9mm)
FIXTNIMGFNM='fix_dot'; %fixation image name.
FIXTNON = single(0); %fixation onset (seconds)
FIXTNDUR_INITL = single(3); %fixation duration (seconds). Should match trial duration.
FIXTNDUR_TRIAL = single(0.5); %fixation duration (seconds). Should match trial duration.
FIXTNPROPXSCRN = single(.5); %proportion of the screen width the centre of the fixation dot
FIXTNPROPYSCRN = single(.5); %proportion of the screen height the centre of the fixation dot
% FIXTNSZ=single(3); %width and height dimension (in cm) for the fixation dot.

% HOMEIMGFLNM='finger_home';
% HOMEPROPXSCRN = FIXTNPROPXSCRN; %finger-home button x-centre
% HOMEPROPYSCRN = single(.5585); %finger-home button y-centre
% HOMESZ=single(.6);
% HOMEON = single(0);
% HOMEDUR = single(.5);

IMGEXT = '.png'; %image file extension
IMGSMPLSZ = 1; %the number of unique images to present from the pool of images.

TARGETOFFSETX15 = (VIEWER_DIST_INCHES * tand(15)) / SCREEN_WIDTH_INCHES;
TARGETOFFSETX30 = (VIEWER_DIST_INCHES * tand(30)) / SCREEN_WIDTH_INCHES;

TARGETIMGPROPXL15 = single(.5 - TARGETOFFSETX15);
TARGETIMGPROPXL30 = single(.5 - TARGETOFFSETX30);
TARGETIMGPROPXR15 = single(.5 + TARGETOFFSETX15);
TARGETIMGPROPXR30 = single(.5 + TARGETOFFSETX30);
TARGETIMGPROPYALL = FIXTNPROPYSCRN;

TRLFDBCKIMG = 'Red_X';

SNDFILES = {'A_snd','B_snd','C_snd','D_snd'};

%STIMDUR = 10; %(1/30)*10; 
TRIALDUR = single(2);

BCKGRNDCLRNM = '0 0 0'; %the background colour
TXTCOLOUR = '1 1 1'; %the text font colour
STIMON = single(0); %stimulus onset varies.
MSKCLR = single(0); %mask colour.
MSKON = single(0); %mask onset.
MSKDUR = single(0); %mask duration.

MRGN = single(1); %margin (in cm) around each quadrant of the screen.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ACQUIRE THE POOL OF IMAGES FROM WHICH TO DRAW SAMPLES FROM   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% imgRootDirPath = uigetdir();
% imgName = ls(fullfile(imgRootDirPath));
% 
% 
% % % imgNmsNeg = ls(fullfile(imgRootDirPath,NEGIMGDIR));
% % % imgNmsNtrlNeg = ls(fullfile(imgRootDirPath,NEUTNEGIMGDIR));
% % % imgNmsNtrlPos = ls(fullfile(imgRootDirPath,NEUTPOSIMGDIR));
% % % imgNmsPos = ls(fullfile(imgRootDirPath,POSIMGDIR));
% % 
% % imgNms(1:2,:)=[]; %delete the '.' and '..' entries
% % % imgNmsNeg(1:2,:)=[]; %delete the '.' and '..' entries
% % % imgNmsNtrlNeg(1:2,:)=[]; %delete the '.' and '..' entries
% % % imgNmsNtrlPos(1:2,:)=[]; %delete the '.' and '..' entries
% % % imgNmsPos(1:2,:)=[]; %delete the '.' and '..' entries
% % 
% % %convert the lists to cell arrays (b/c the file names differ in length) and remove the
% % %file extension (no need for the file extension in the config file).
% % [imgNegPopltnSz,c]=size(imgNmsNeg);
% % imgNmsNeg = strtrim(mat2cell(imgNmsNeg,ones(imgNegPopltnSz,1),c));
% % 
% % [imgNtrlNegPopltnSz,c]=size(imgNmsNtrlNeg);
% % imgNmsNtrlNeg = strtrim(mat2cell(imgNmsNtrlNeg,ones(imgNtrlNegPopltnSz,1),c));
% % 
% % [imgNtrlPosPopltnSz,c]=size(imgNmsNtrlPos); 
% % imgNmsNtrlPos = strtrim(mat2cell(imgNmsNtrlPos,ones(imgNtrlPosPopltnSz,1),c));
% % 
% % [imgPosPopltnSz,c]=size(imgNmsPos);
% % imgNmsPos = strtrim(mat2cell(imgNmsPos,ones(imgPosPopltnSz,1),c));
% % 
% % %check to make sure there are the same number of images in each category. 
% % if imgNtrlNegPopltnSz ~= imgNtrlPosPopltnSz || imgNegPopltnSz ~= imgPosPopltnSz || ...
% %         imgNtrlNegPopltnSz ~= imgNegPopltnSz || imgNtrlNegPopltnSz ~= imgPosPopltnSz
% %     disp('Error! SPIDERSTUDYCONFIGFILEMAKER: Number of Images in each subfolder must be equal!');
% %     return;
% % end
% % 
% % %remove the image file extension for the lists of image file names...
% % for j=1:imgNtrlNegPopltnSz
% %     imgNmsNeg{j}=cell2mat(strsplit(imgNmsNeg{j},IMGEXT));
% %     imgNmsNtrlNeg{j}=cell2mat(strsplit(imgNmsNtrlNeg{j},IMGEXT));
% %     imgNmsNtrlPos{j}=cell2mat(strsplit(imgNmsNtrlPos{j},IMGEXT));
% %     imgNmsPos{j}=cell2mat(strsplit(imgNmsPos{j},IMGEXT));
% % end

%for each desired trial order...
for i=1:numTrialOrders

    %create the condition order. 1-4 = stim dots A,B,C,D
    cndnOrder=CreateFactorialOrder([5 4],tpc,repLim);

    %store the trial total
    [trialTot,~] = size(cndnOrder);

%     %'B' stores the condition information to save to disk as a seperate file
% %     B=cell(2,trialTot);
%     B=cell(12,trialTot);
% 
%     %overall trial number
%     B(1,:)=num2cell(transpose(uint8(cndnOrder(:,1))));
% 
%     %stim dot number
%     B(4,:)=stim_dot_letters(cndnOrder(:,2));
%     
%     %jump vs. no jump
%     B(4,:)=num2cell(transpose(cndnOrder(:,2)>3));
% 
%     %valence switch
%     B(7,:)=num2cell(transpose(and(cndnOrder(:,2)>3,or(cndnOrder(:,3)==2,cndnOrder(:,3)==3))));

%     %initial image position
%     temp = cndnOrder(:,2);
%     temp(temp==4)=2; %convert the jump to position 1 condition to an initial position 2
%     temp(temp==5)=2; %convert the jump to position 3 condition to an initial position 2
%     B(8,:)=num2cell(transpose(uint8(temp)));
% 
%     %final image position
%     temp = cndnOrder(:,2);
%     temp(temp==4)=1; %convert the jump to position 1 condition to a final position 1
%     temp(temp==5)=3; %convert the jump to position 3 condition to a final position 3
%     B(9,:)=num2cell(transpose(uint8(temp)));

    %remove the trial number column since it is no longer required.
    cndnOrder(:,1)=[];

    %convert dpi to dots per cm.
    %dotsPerCm = dpi/CMPERINCH;

    %acquire some useful info.
    [FIELDNAMES,HEADERLINE,SAVEFORMATSPEC,~,CNDFIELDNAMES,CNDHEADERLINE,CNDSAVEFORMATSPEC,~] = ...
        SaccadeStudyConfigFileMakerAssist;
    numFields = length(FIELDNAMES);

    %check that the number of trials is evenly divisble by 'blocks'
    if mod(trialTot,blocks)~=0
        disp(['Error! SPIDERSTUDYCONFIGFILEMAKER: check that trial total will be evenly '...
            'by ''blocks'''])
        return;

    %otherwise, the total number of trials is divisible by the number of requested blocks.
    %Determine the number of trials in each block.
    else
        trialsPerBlock = trialTot/blocks;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  GENERATE THE RAPID SERIAL VISUAL PRESENTATION OF ALPHANUMERIC CHARACTERS  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     %determine the number of alphanumeric characters for each series (20-25 characters)
%     numChrctrsPerTrial=randi(6,[1 trialTot])+19;
%     alphasPerTrial=randi(6,[1 trialTot]);
%     alphaNumericSequences = cell(25,trialTot);

    %the number of events (images to present) per trial. fixation + home-buttonx2 + max of 25
    %alphanumeric characters + fixation again + initial target image + final target image.
    %Thus, 31 events per trial is the max.
    numEventsPerTrial=5;
    
    %determine the block number list
    blockRange = [1:blocks]';
    blRngRep = blockRange(:,ones(trialsPerBlock*numEventsPerTrial,1))';
    blockNmLst = blRngRep(:)';

    %the list of trial numbers within the block
    trialNmWthnBlckLst=mod(1:trialTot,trialsPerBlock);
    trialNmWthnBlckLst(trialNmWthnBlckLst == 0) = trialsPerBlock;
    trialNmWthnBlckLst = transpose(trialNmWthnBlckLst);
    trNmsRep = trialNmWthnBlckLst(:,ones(numEventsPerTrial,1))';
    trialNmWthnBlckLst = trNmsRep(:)';
%     trialNmWthnBlckLst=[];
    onsetOfEventsPerTrial=cell(trialTot,1);
    durationOfEventsPerTrial=cell(trialTot,1);

%     %determine the alphanumeric sequence itself. Final series is always 14 digits long.
%     for j=1:trialTot
% 
%         %number of digits to appear first depends on characters per trial minus 15 (the target
%         %and the 14 digits that always follow)
%         numFirstDigitsPerTrial(j) = numChrctrsPerTrial(j)-15;
%         firstDigits = randi(9,[numFirstDigitsPerTrial(j) 1]);
%         while sum(diff(firstDigits)==0)~=0
%             firstDigits = randi(9,numFirstDigitsPerTrial(j));
%         end
% 
%         %always 14 digits presented after the target...
%         endDigits = randi(9,[14 1]);
%         while sum(diff(endDigits)==0)~=0
%             endDigits = randi(9,[14 1]);
%         end
% 
%         %put the character sequences together
%         %store the first digits
%         r=1;
%         for k=1:numFirstDigitsPerTrial(j)
%             alphaNumericSequences{k,j}=num2str(firstDigits(k));
%             r=r+1;
%         end
%         targetLetters(j) = ALPHAS(randi(6,1));
%         alphaNumericSequences{r,j}=targetLetters(j); %add the target letter
%         r=r+1;
%         
%         %store the end digits (always 14 of them)
%         for k=1:14
%             alphaNumericSequences{r,j}=num2str(endDigits(k)); %add the remaining digits
%             r=r+1;
%         end
% 
%         %number of events per trial includes fixation, home button (x2), alphanumeric 
%         %character events, return to the fixation-image, the initial target image, and the 
%         %final target image 
%         numEventsPerTrial(j)=1+2+numChrctrsPerTrial(j)+1+1+1;
%         trialNmWthnBlckLst=[trialNmWthnBlckLst; j*ones(numEventsPerTrial(j),1)];
%         initialTargetImgOnset = (numFirstDigitsPerTrial(j)+1+randi(3))*.1;
% 
%         %-1 for final image onset dependence on the finger-release
%         onsetOfEventsPerTrial{j}=single([0 0 .0167 (.1:.1:(numChrctrsPerTrial(j)+1)*.1) ...
%             initialTargetImgOnset -1]);
% 
%         %-1 for initial image duration dependence on the finger-release and -3 for final image
%         %dependence on max trial time.
%         durationOfEventsPerTrial{j}=single([0 0 initialTargetImgOnset+.0835 ...
%             repmat(.03,[1 numChrctrsPerTrial(j)]) TRIALDUR-(numChrctrsPerTrial(j)*.1) -1 -3]); 
%     end
    
%     %store the trial # w/n block
%     B(3,:)=num2cell(uint8(repmat(1:trialsPerBlock,[1 blocks])));
%     
%     %store the target letters.
%     B(12,:)=transpose(cellstr(transpose(targetLetters))); 

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %   SAMPLE FROM THE POOL OF IMAGES   %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     %pear down the list of image file-names to a sample of images from each pool, then 
%     %create a n (image sample size) x 5 conditions (3 no-jump positions + 2 jump trial 
%     %positions = 5) x 2 (page 1: initial image; page 2: final image). The levels of 
%     %factor 1 (col 1 in 'cndnOrder') determine the columns in the image lists used to 
%     %populate the event. The levels of factor 2 ('col 2 in 'cndnOrder') determine the
%     %image lists used to populate the event (1: neg; 2: neut neg; 3: neut pos;
%     %4:pos). This setup requires indices to keep track of which images have already been
%     %used each to march down the rows of a given column/image list.
%     
%     %select the same random subset of images to stock the negative and neutral negative
%     %samples
%     selection=randperm(imgNtrlNegPopltnSz,IMGSMPLSZ);  
%     imgNms = repmat(imgNmsNeg(selection,:),[1 5]); %negative image names (page 1)
%     imgNms(:,:,2) = repmat(imgNmsNtrlNeg(selection,:),[1 5]); %neutral neg image names (page 2)
% 
%     %select the same random subset of images to stock the positive and neutral positive
%     %samples
%     selection=randperm(imgNtrlNegPopltnSz,IMGSMPLSZ);
%     imgNms(:,:,3) = repmat(imgNmsNtrlPos(selection,:),[1 5]); %positive image names (page 3)
%     imgNms(:,:,4) = repmat(imgNmsPos(selection,:),[1 5]); %neutral pos image names (page 4)
% 
%     %no-jump trials ('j' increments through image sets for positions 1, 2, and 3)
%     for j=1:3
%         for k=1:4
%             imgNms1(:,j,k) = imgNms(randperm(IMGSMPLSZ,IMGSMPLSZ),j,k);
%         end
%     end
%     
%     imgNms2 = imgNms1;
%     
%     %create the jump trials for the intial and final image sets
%     for j=4:5
%         selection=randperm(IMGSMPLSZ,IMGSMPLSZ);
%         imgNms1(:,j,1) = imgNms(selection,j,1); %neg initial
%         imgNms2(:,j,1) = imgNms(selection,j,1); %neg final
%         
%         selection=randperm(IMGSMPLSZ,IMGSMPLSZ);
%         imgNms1(:,j,2) = imgNms(selection,j,2); %neutr neg initial 
%         imgNms2(:,j,2) = imgNms(selection,j,1); %neg final
%         
%         selection=randperm(IMGSMPLSZ,IMGSMPLSZ);
%         imgNms1(:,j,3) = imgNms(selection,j,3); %neutr pos initial 
%         imgNms2(:,j,3) = imgNms(selection,j,4); %pos final
%         
%         selection=randperm(IMGSMPLSZ,IMGSMPLSZ);
%         imgNms1(:,j,4) = imgNms(selection,j,4); %pos initial 
%         imgNms2(:,j,4) = imgNms(selection,j,4); %pos final 
%     end
%     
%     %indices to denote the current image for each of the 5 conditions (cols). Each
%     %category of image has its own matrix.
%     img1Indcs=ones(1,5,4); %first image
%     img2Indcs=ones(1,5,4); %second image
% 
    %total number of rows in the config file, including the header row
    rowsTotCnfgFl = numEventsPerTrial*trialTot;

    %determine the block number list
    blockRange = [1:blocks]';
    blRngRep = blockRange(:,ones(trialsPerBlock*numEventsPerTrial,1))';
    blockNmLst = blRngRep(:)';
%     blockNmLst = nan(rowsTotCnfgFl,1);
%     index=1;
%     for j = 1:blocks
%         numEventsInCurBlck = sum(numEventsPerTrial(1+(j-1)*trialsPerBlock:j*trialsPerBlock));
%         blockNmLst(index:index-1+numEventsInCurBlck)=j*ones(numEventsInCurBlck,1);
%         index = index+numEventsInCurBlck;
%     end

%     %store the block #
%     B(2,:) = num2cell(uint8(reshape(repmat(1:blocks,[trialsPerBlock 1]),[1 trialTot]))); 
%    
    %we'll use FPRINT, so rows in 'A' will end up as columns (and cols in 'A' end up rows)
    %when 'A' is saved to disk.
    A = cell(numFields,rowsTotCnfgFl+1);
   
    %column index for A
    ACol = 2;
    
    %populate 'A' for FPRINTF
    for j = 1:trialTot 

        %loop through each event in the current trial
        for k=1:numEventsPerTrial

            A{1,ACol} = uint8(blockNmLst(ACol-1));          %1: block number
            A{2,ACol} = uint8(0);                           %2: is_practice_block?
            A{3,ACol} = uint8(0);                           %3: block_has_feedback?
            A{4,ACol} = single(0);                          %4: block_feedback_thresh
            A{5,ACol} = uint16(trialNmWthnBlckLst(ACol-1)); %5: trial_num within the block
            A{6,ACol} = TRIALDUR;                           %6: trial_max_time
            A{7,ACol} = uint8(1);                           %7: trial_timeout_msg
            A{8,ACol} = uint8(0);                       %8: trial_kb_resp?
            A{9,ACol} = 'None';                         %9: correct_kb_resp
            A{10,ACol} = uint8(1);                          %10: trial_feedback_type
            A{11,ACol} = TRLFDBCKIMG;                       %11: trial_feedback_img
            A{12,ACol} = BCKGRNDCLRNM;                      %12: background_colour RGB
            A{13,ACol} = TXTCOLOUR;                         %13: text_colour

            A{29,ACol} = single(0);                         %28: mask_onset (seconds)
            A{30,ACol} = single(0);                         %29: mask_duration (seconds)
            A{31,ACol} = single(0);                         %30: mask_size
            A{32,ACol} = uint16(MSKCLR);                    %31: mask_color
            A{33,ACol} = single(0);                         %32: mask_rotation (degrees)
            A{34,ACol} = uint8(0);                          %33: mask_fit
            A{35,ACol} = single(0);                         %34: mask_margin
            
            %fixation image
            if k==1            
                A{14,ACol} = FIXTNIMGFNM;                       %14: stim_img_name
                A{15,ACol} = 'None';                            %15: stim_sound_file
                A{16,ACol} = FIXTNON;                           %16: stim_onset (s)
                A{17,ACol} = single(0);                         %17: stim_duration (s)
                A{18,ACol} = FIXTNPROPXSCRN; %18: stim_cent_x (proportion of the screen width)
                A{19,ACol} = FIXTNPROPYSCRN; %19: stim_cent_y (proportion of the screen height)
                A{20,ACol} = FIXIMGDISPSZXCM;                   %20: stim_size_x (cm)
                A{21,ACol} = FIXIMGDISPSZYCM;                   %21: stim_size_y (cm)
                A{22,ACol} = single(0);                         %22: stim rotation
                A{23,ACol} = uint8(0);                          %23: stim_is_touchable
                A{24,ACol} = uint8(0);                          %24: stim_is_target
                A{25,ACol} = uint8(2);                          %25: subj_fixation_type
                A{26,ACol} = uint8(1);                          %26: subj_fixation_pause
                A{27,ACol} = FIXTNON;                           %27: subj_fixation_onset (s)
                if trialNmWthnBlckLst(ACol-1) == 1
                    A{28,ACol} = FIXTNDUR_INITL;                %28: subj_fixation_duration (s)
                else
                    A{28,ACol} = FIXTNDUR_TRIAL;                %28: subj_fixation_duration (s)
                end
                
            %target dots
            else
                A{14,ACol} = FIXTNIMGFNM;                       %14: stim_img_name
                
                %15: stim_sound_file
                if (k-1) == cndnOrder(j,2) % stim_is_target  
                    if (k-1)==1                         %position A
                        A{15,ACol} = SNDFILES{1};                         
                    elseif (k-1)==2                     %position B
                        A{15,ACol} = SNDFILES{2};           
                    elseif (k-1)==3                     %position C
                        A{15,ACol} = SNDFILES{3}; 
                    else                                %position D
                        A{15,ACol} = SNDFILES{4}; 
                    end
                else
                    A{15,ACol} = 'None';
                end
                
                A{16,ACol} = single(0);                         %16: stim_onset (s)
                A{17,ACol} = single(-3);                        %17: stim_duration (s)

                %18: stim_cent_x (proportion of the screen width)    
                if (k-1)==1                         %position A
                    A{18,ACol} = TARGETIMGPROPXL30;                         
                elseif (k-1)==2                     %position B
                    A{18,ACol} = TARGETIMGPROPXL15;           
                elseif (k-1)==3                     %position C
                    A{18,ACol} = TARGETIMGPROPXR15; 
                else                                %position D
                    A{18,ACol} = TARGETIMGPROPXR30; 
                end

                %19: stim_cent_y (proportion of the screen height)
                A{19,ACol} = TARGETIMGPROPYALL;
                A{20,ACol} = TARGETIMGDISPSZXCM;                %20: stim_size_x (cm)
                A{21,ACol} = TARGETIMGDISPSZYCM;                %21: stim_size_y (cm)
                A{22,ACol} = single(0);                         %22: stim rotation
                A{23,ACol} = uint8(2);                          %23: stim_is_touchable
                
                %24: stim_is_target   
                if (k-1) == cndnOrder(j,2)
                    A{24,ACol} = uint8(2);
                else
                    A{24,ACol} = uint8(0);
                end
                
                A{25,ACol} = uint8(0);                          %25: subj_fixation_type
                A{26,ACol} = uint8(0);                          %26: subj_fixation_pause
                A{27,ACol} = single(0);                         %27: subj_fixation_onset (s)
                A{28,ACol} = single(0);                         %28: subj_fixation_duration (s)
                
%                 %store the final image name
%                 B{11,j}=A{14,ACol};
            end

            %increment ACol
            ACol=ACol+1;
        end
    end

    %save the header line to file, store 'A' to file, then close the file.
    if ~isempty(cnfFNm)
        
%         %add single vs. dual task designation to the file name
%         if dualTask==0
%             sfxFlNm='_sing';
%         else
%             sfxFlNm='_dual';
%         end
%         
        %add the trial order number to the file name
        if i<10
            sfxFlNm=['_0' num2str(i)];
        elseif i>9 && i<100
            sfxFlNm=[num2str(i)];
        end
          
        %save the file to disk    
        fid = fopen([cnfFNm sfxFlNm '.csv'],'w');
        fprintf(fid,HEADERLINE);
        fprintf(fid,SAVEFORMATSPEC,A{:,2:end});
        fclose(fid);
        A=transpose(A);
        
%         %save a csv file of the condition order and other relevant information.
%         fid = fopen([cnfFNm sfxFlNm '_cond_file.csv'],'w');
%         fprintf(fid,CNDHEADERLINE);
%         fprintf(fid,CNDSAVEFORMATSPEC,B{:,:});
%         fclose(fid);
%         B=[CNDFIELDNAMES B];
%         B=transpose(B);
    end
end
return;