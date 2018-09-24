%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shape Tapper Experiment Framework
%
% Displays image stimuli and dot masks at user-specified times.
% Allows for mouse, touch and SMI EyeLink feedback while recording motion
% capture data from the NDI Optotrack
%
% Config file format designed for flexibility in defining new behaviours.
% See ShapeTapper_generateConfig for use of config files.
%
% -Jamie Dunkle, UBC Vision Lab, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ShapeTapper_test(config_fname)

% Clear the workspace and the screen
sca;
close all;
clearvars;

% Set messages to participant
welcomeMsg = ['Welcome to the experiment!\n\n\n',...
              'When you see a cross (+) in the center of the screen\n',...
              'press the spacebar to start each trial\n\n\n',...
              'Press Any Key To Begin\n',...
              'Esc to exit at any time'];
          
blockEndMsg = 'Block Finished\n\n\nPress Any Key To Continue';

expEndMsg = 'Experiment Finished\n\n\nPress Any Key To Exit';

% Acceptable stimulus image formats, must be compatible with imread()
img_formats = {'.png', '.jpg'};


% Subject name
subj = 'test';

% Specify directory containing stimulus images
stim_dir = 'Image_Files\';

% Specify directory containing config files
switch nargin
    case 0
        [config_fname,config_path] = uigetfile('./*.csv',...
                                       'Select an experiment config file');
        config_fname = [config_path config_fname];
    case 1
        if isempty(config_fname) || isnan(config_fname)
            config_fname = 'Config_Files\test_exp1.csv';
        end
end


%----------------------------------------------------------------------
%                       PyschToolbox Setup
%----------------------------------------------------------------------
% Run function to init PsychToolbox screen and store relevant vars for
% later rendering and timing
ptb = ptb_initscreen();

% Get window handle for PTB rendering
window = ptb.window;

% Get monitors Inter-Frame Interval and waittime (# of frames b/t renders)
ifi = ptb.ifi;
waitframes = ptb.waitframes;

% Flip to clear
Screen('Flip', window);

%----------------------------------------------------------------------
%                     Stimuli - Load Config File
%----------------------------------------------------------------------

% Read specified config file for this experiment
config_dat = readtable(config_fname);


% -- LOAD IMAGES AS PSYCHTOOLBOX TEXTURES

% Retreive names of stimuli images
img_names = unique(config_dat.stim_img_name);

% Extract unique strings for a list of stims
img_names = unique(img_names);

% Load each unique stimulus image as a PsychToolbox texture handle 
% and store handles in a Map object to render during trials
stim_textures = ptb_loadtextures(stim_dir, img_names, img_formats, ptb);
% stim_textures = containers.Map;
% for i=1:length(img_names)
%     % Retrieve name of stimulus image
%     img_name = char(img_names(i));
%     img_fname = [];
%     
%     for ff = 1:length(img_formats)
%         fn = [stim_dir img_name char(img_formats{ff})];
%         if exist(fn, 'file') == 2
%             img_fname = fn;
%         end
%     end
%     
%     if isempty(img_fname)
%         sca;
%         error(['Could not load stimulus image \"%s\" \n' ...
%                'Please check your config file and stimulus directory.\n'...
%                'Exiting...'], img_name);
%     end
%     
%     % Load image using imread(), save colormap and alpha channel
%     % for further processing if needed
%     [img, map, alpha] = imread(img_fname);
%    
%     % Create white background (may be changed to another value if needed)
%     white_bg = (ones(size(img))*255);
%     
%     % Scale alpha channel and apply background to transparent part of image
%     alphaMask = im2double(alpha); %scale between 0 and 1
%     img = im2uint8(white_bg.*repmat((1-alphaMask),[1 1 3]) + double(img).*repmat(alphaMask,[1 1 3]));
%     
%     % Make the image into a texture, save object handle as property
%     h = Screen('MakeTexture', window, img);
%     
%     % Add image name and texture handle to Map for later display
%     stim_textures(img_name) = h;
% end


%----------------------------------------------------------------------
%                     Make an output data table
%----------------------------------------------------------------------

% Timing:
%       stim presentation times (requested vs actual)
%       average frame rate, max framerate error
     
% Parse and count the total number of trials specified
tot_num_trials = 0;
blocks = unique(config_dat.block_num);
num_blocks = length(blocks);
for b=1:num_blocks
     % Collect all block numbers
    bn = blocks(b);
    
    % Retrieve block rows for current block
    block_dat = config_dat(config_dat.block_num == bn, :);
    
    % Retreive list of trials
    trial_list = unique(block_dat.trial_num);

    % Calculate total number of trials
    tot_num_trials = tot_num_trials + length(trial_list);
end

% Initialize output table
data_columns = {'Subject',...
                'BlockNum',...
                'TrialNum',...
                'BadTrial',...
                'target_mask_name',...
                'target_center_x',...
                'target_center_y',...
                'target_rotation',...
                'target_mask_num',...
                'chosen_mask_num',...
                'choice_correct',...
                'touch_time',...
                'touch_point_x',...
                'touch_point_y'};

subjData = cell2table(cell(tot_num_trials, length(data_columns)),...
                      'VariableNames',data_columns);

% Counter for correct row indexing when writing output
row_ct = 1;
                  
% Set output file name
out_fname = [subj '_out.csv'];

%----------------------------------------------------------------------
%                       Experimental loop
%----------------------------------------------------------------------

% Animation Loop -- here we render and display all blocks & trials
% specified in the config file, and save results

% BLOCK LOOP
for b=1:num_blocks
    
    % If this is the very  first trial, present a start screen and wait
    % for a key-press
    if b == 1
        DrawFormattedText(window, welcomeMsg, 'center', 'center', ptb.black);
        Screen('Flip', window);

        % Listen for Esc key to abort experiment
        [~, keyCode, ~] = KbStrokeWait;
        if keyCode(ptb.escapeKey)
            sca;
            return;
        end
    end
    
    % Collect all block numbers
    bn = blocks(b);
    
    % Retrieve block rows for current block
    block_dat = config_dat(config_dat.block_num == bn, :);
    
    % Retreive list of trials
    trial_list = unique(block_dat.trial_num);

    % Calculate total number of trials
    num_trials = length(trial_list);

    % TRIAL LOOP
    for t=1:num_trials

        % -- Pre-trial initialization: event schedules, calculate masks --

        % Collect all trial numbers
        tn = trial_list(t);
        
        % Retrieve stimuli list for current trial
        trial_dat = block_dat(block_dat.trial_num == tn, :);
        
        % Convert cell array of stim names to char arrays
        stim_img_names = cell2mat(trial_dat.stim_img_name);
        
        % Calculate total number of frames for current trial
        trial_frames  = ceil(max(trial_dat.trial_max_time) / ifi);

        % Make event schedule arrays, which will allow for rapid checking
        % of which stimuli to render for each frame
        [num_stims, ~] = size(stim_img_names);
        stim_schedule = zeros(num_stims, trial_frames);
        mask_schedule = zeros(num_stims, trial_frames);

        % Initialize bounding rectangle for each stimulus
        stim_bounds = zeros(num_stims, 4);

        % Initialize mask arrays for each stimlus
        dotPosMatrix = cell(1,num_stims);
        dotSizes = cell(1,num_stims);
        dotColors = cell(1,num_stims);
        dotBoundingRect = cell(1,num_stims);

        % Declare target mask
        targetMask = NaN;
        
        % Populate stimulus schedules and mask vectors
        for s=1:num_stims
            % Extract stimulus properties from trial row
            st = trial_dat(s,:);
            
            % Store target mask if it matches stim description
            % NOTE: assumption (FOR NOW) is that only 1 stim is the target
            if st.stim_is_target
                targetMask = s;
            end
            
            % Fill all frames where this stim is displayed with a 1
            stim_img_on = ceil(st.stim_onset / ifi);
            stim_img_off = ceil((st.stim_onset+st.stim_duration) / ifi);
            stim_schedule(s, stim_img_on:stim_img_off) = 1;

            % Fill all frames where this stim's mask is displayed with a 1
            mask_img_on = ceil(st.mask_onset / ifi);
            mask_img_off = ceil((st.mask_onset+st.mask_duration) / ifi);
            mask_schedule(s, mask_img_on:mask_img_off) = 1;

            % Draw rectangle at the requested stimulus size
            rect = [0 0 st.stim_size_x st.stim_size_y] * ptb.ppcm;

            % Offset rectangle to stimulus center point
            rect = CenterRectOnPoint(rect, st.stim_cent_x, st.stim_cent_y);

            % Store rectangle in stim_bounds array for texture rendering
            stim_bounds(s,:) = rect;

            % Make a dotPositionMatrix, dotSize, dotColor array for mask
            if st.mask_fit == 1
                % Draw mask dots fitted to shape bounds
            else
                % Draw mask dots that show the rect's coordinates
                mask_pos = [rect(1), rect(1), rect(3), rect(3);...
                            rect(2), rect(4), rect(2), rect(4)];
            end

            % Calculate the number of dots For help see: help numel
            mask_ndots = length(mask_pos);

            % Set the color of our dots (alpha channel defaults to 1)
            mask_dotcolors = [ones(1,3).*st.mask_color 1];

            % Set the size of the dots
            mask_dotsizes = ones(1, mask_ndots) .* st.mask_size;


            %------ Scale and offset mask points according to input params

            % Center mask positions relative to stimulus location
            [~,c]=size(mask_pos);
            mask_pos = mask_pos - [repmat(st.stim_cent_x,1,c); repmat(st.stim_cent_y,1,c)];

            % -- Scale mask dots to meet requested margin

            % Get mask diagonal in inches
            mask_diag = sqrt((st.stim_size_x)^2 + (st.stim_size_y)^2);

            % Get margin diagonal
            margin_diag = sqrt(2 * (st.mask_margin)^2);

            % Calculate scaling factor and apply to mask points
            scale = (mask_diag + margin_diag) / mask_diag;
            mask_pos = mask_pos * scale;

            % Create a rotation matrix
            rot = [cosd(st.mask_rotation) -sind(st.mask_rotation);...
                   sind(st.mask_rotation) cosd(st.mask_rotation)];

            % Rotate mask dots around center to match requested rotation
            mask_pos = rot * mask_pos;

            % Reapply offset to mask positions
            [~,c]=size(mask_pos);
            mask_pos = mask_pos + [repmat(st.stim_cent_x,1,c); repmat(st.stim_cent_y,1,c)];

            % Store mask in cell array for rendering in trial presentation
            dotPosMatrix{s} = mask_pos;
            dotSizes{s} = mask_dotsizes;
            dotColors{s} = mask_dotcolors;
            dotBoundingRect{s} = rect;
        end % stim loop

        % -- Display trial stimuli, record responses ---------------------

        % Draw the fixation cross in black, set it to the center of our
        % screen and set good quality antialiasing
        Screen('DrawLines', window, allCoords,...
            lineWidthPix, ptb.black, [ptb.xCenter ptb.yCenter], 2);

        % Flip again to sync us to the vertical retrace
        vbl = Screen('Flip', window);

        % Display fixation cross and wait for trial to begin
        % (can be set to key press, eye on fixation, finger on screen, etc)
        trialWait = true;
        while trialWait == true
            % Draw the fixation cross in black, set it to the center of
            % our screen and set good quality antialiasing
            Screen('DrawLines', window, allCoords,...
                lineWidthPix, ptb.black, [ptb.xCenter ptb.yCenter], 2);

            % Listen for Esc key to abort experiment, Space to begin
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(ptb.escapeKey)
                sca;
                return
            elseif keyCode(ptb.spaceKey)
                trialWait = false;
            end

            % Flip again to sync us to the vertical retrace
            vbl = Screen('Flip', window);
        end

        % Init time stamp (will be overwritten on first Screen flip)
        rtStart = NaN;
        max_mouse_samples = round(10 / ifi);
        mouseData = zeros(3, max_mouse_samples);
        mouseCounter = 1;

        % Flip again to sync us to the vertical retrace
        vbl = Screen('Flip', window);

        % Cue to determine whether a response has been made
        hasSelected = false;

        % Save response time (no response = NaN)
        rt = NaN;
        
        % Display all stimuli which are scheduled for the current frame
        for frame=1:trial_frames

            % Listen for Esc key to abort experiment
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(ptb.escapeKey)
                sca;
                return
            end

            % Get the current position of the mouse
            [tx, ty, buttons] = GetMouse(window);

            % We clamp the values at the maximum values of the screen
            % in X and Y in case people have two monitors connected.
            tx = min(tx, screenXpixels);
            ty = min(ty, screenYpixels);

            % Store mouse position sample in mouseData array
            if mouseCounter < max_mouse_samples
                mouseData(1, mouseCounter) = tx;
                mouseData(2, mouseCounter) = ty;
                mouseData(3, mouseCounter) = vbl;
                mouseCounter = mouseCounter + 1;
            end

            % Check image and mask schedules and render anything to screen
            % that is scheduled to be shown on the current frame
            for s=1:num_stims
                if stim_schedule(s, frame) == 1
                    stim_name = char(trial_dat.stim_img_name(s));
                    Screen('DrawTextures', window,...
                            stim_textures(stim_name),[],...
                            stim_bounds(s,:), trial_dat.stim_rotation(s));
                end % image scheduled
                if mask_schedule(s, frame) == 1
                    % Draw mask dots
                    Screen('DrawDots', window, dotPosMatrix{s},...
                            dotSizes{s}, dotColors{s}, [], 2);
                    % Check for collisions.
                    if trial_dat.mask_touchable(s) == 1
                        if IsInRect(tx, ty, dotBoundingRect{s})
                            hasSelected = true;
                            selectedMask = s;
                        end % touched stim mask
                    end % stim mask touchable
                end % mask scheduled
            end % stim loop

            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

            % Note start time to calculate reaction time
            if frame == 1
                rtStart = vbl;
            end

            % If a touchable mask has been touched, end trial immediately
            if hasSelected == true
                rt = vbl - rtStart;
                correct_choice = (selectedMask == targetMask);
                break
            end

        end % frame loop
        
        % -- Output data for current trial --------------------------------

        % Extract target stimulus information from trial stim table
        tr = trial_dat(targetMask,:);
        
        % Record the trial data into out data matrix
        subjData.Subject{row_ct} = subj;
        subjData.BlockNum{row_ct} = bn;
        subjData.TrialNum{row_ct} = row_ct;
        subjData.BadTrial{row_ct} = 0; % TBD
        subjData.target_mask_name{row_ct} = cell2mat(tr.stim_img_name);
        subjData.target_center_x{row_ct} = tr.stim_cent_x;
        subjData.target_center_y{row_ct} = tr.stim_cent_y;
        subjData.target_rotation{row_ct} = tr.stim_rotation;
        subjData.target_mask_num{row_ct} = targetMask;
        subjData.chosen_mask_num{row_ct} = selectedMask;
        subjData.choice_correct{row_ct} = correct_choice;
        subjData.touch_time{row_ct} = rt;
        subjData.touch_point_x{row_ct} = tx;
        subjData.touch_point_y{row_ct} = ty;
        row_ct = row_ct + 1;
    end % trial loop

    % Flip again to sync us to the vertical retrace
    vbl = Screen('Flip', window);

    if bn < num_blocks
        % End of block screen. 
        % We clear the screen once they have made their response
        DrawFormattedText(window,blockEndMsg,...
                            'center', 'center', black);
        Screen('Flip', window);
        KbStrokeWait;
    end

end % block loop

% Flip again to sync us to the vertical retrace
vbl = Screen('Flip', window);

% End of experiment screen. We clear the screen once they have made their
% response
DrawFormattedText(window, expEndMsg, 'center', 'center', black);
Screen('Flip', window);
KbStrokeWait;
sca;

% Output data table to file
writetable(subjData, out_fname, 'Delimiter', ',', 'QuoteStrings', true);