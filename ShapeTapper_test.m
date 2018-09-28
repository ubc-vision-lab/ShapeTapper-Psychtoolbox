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
              'press and hold on the cross until the stimuli appear.\n',...
              'Touch the stimuli which appears to be upside down\n\n\n'...
              'Press Any Key To Begin\n',...
              'Esc to exit at any time'];
          
blockEndMsg = 'Block Finished\n\n\nPress Any Key To Continue';

expEndMsg = 'Experiment Finished\n\n\nPress Any Key To Exit';

timeoutMsg = 'Trial Time Expired.\n\n\nPress Any Key To Continue';

% Acceptable stimulus image formats, must be compatible with imread()
img_formats = {'.png', '.jpg'};


% Subject name
subj = 'test';

% Specify directory containing stimulus images
stim_dir = 'Image_Files\';

% Specify directory containing config files
switch nargin
    case 0
        [config_fname,config_path] = uigetfile('./Config_Files/*.csv',...
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


%----------------------------------------------------------------------
%                     Make an output data table
%----------------------------------------------------------------------
     
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
                'chosen_mask_name',... 
                'chosen_center_x',...
                'chosen_center_y',...
                'chosen_rotation',...
                'target_mask_num',...
                'chosen_mask_num',...
                'choice_correct',...
                'touch_point_x',...
                'touch_point_y',...
                'trial_time_elapsed',...
                'rt_target_to_choice',...
                'rt_target_to_saccade',...
                'rt_saccade_to_choice',...
                'rt_target_to_fingerlift',...
                'rt_fingerlift_to_choice'};

subjData = cell2table(cell(tot_num_trials, length(data_columns)),...
                      'VariableNames',data_columns);

% Counter for correct row indexing when writing output
row_ct = 1;
           
out_path = 'Data/';

% Set output file name
subj_data_fname = [out_path subj '_out.csv'];
fid = fopen( subj_data_fname, 'wt+' );
if fid>0
    fprintf(fid,strjoin(data_columns,','));
    fprintf(fid,'\n');
    fclose(fid);
end

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
        DrawFormattedText(window,welcomeMsg,'center','center',ptb.black);
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
        
        % Calculate total number of frames for current trial
        trial_frames  = ceil(max(trial_dat.trial_max_time) / ptb.ifi);

        % Convert cell array of stim names to char arrays
        stim_img_names = cell2mat(trial_dat.stim_img_name);
        
        % Make event schedule arrays, which will allow for rapid checking
        % of which stimuli to render for each frame
        [num_stims, ~] = size(stim_img_names);
        stim_schedule = zeros(num_stims, trial_frames);
        mask_schedule = zeros(num_stims, trial_frames);
        fixation_schedule = zeros(num_stims, trial_frames);
        
        % Initialize counters for fixation durations
        fixation_counter = zeros(num_stims, 3);

        % Initialize bounding rectangle for each stimulus
        stim_bounds = zeros(num_stims, 4);

        % Initialize mask arrays for each stimlus
        dotPosMatrix = cell(1,num_stims);
        dotSizes = cell(1,num_stims);
        dotColors = cell(1,num_stims);
        dotBoundingRect = cell(1,num_stims);

        % Declare target mask
        targetMask = NaN;
        selectedMask = NaN;
        correct_choice = false;
        
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
            stim_img_on = ceil(st.stim_onset / ifi) + 1;
            stim_img_off = ceil((st.stim_onset+st.stim_duration) / ifi) + 1;
            if stim_img_on ~= stim_img_off
                stim_schedule(s, stim_img_on:stim_img_off) = 1;
            end

            % Fill all frames where this stim's mask is displayed with a 1
            mask_img_on = ceil(st.mask_onset / ifi) + 1;
            mask_img_off = ceil((st.mask_onset+st.mask_duration) / ifi) + 1;
            if mask_img_on ~= mask_img_off
                mask_schedule(s, mask_img_on:mask_img_off) = 1;
            end
            
            % If stim is a fixation event, save in fixation schedule
            if st.subj_fixation_type ~= 0
                fixation_schedule(s, st.subj_fixation_onset + 1) = st.subj_fixation_type;
                fixation_frames = ceil(st.subj_fixation_duration / ifi);
                fixation_counter(s, 1) = fixation_frames;
            end
            
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
            mask_dotsizes = ones(1, mask_ndots) .* max(1, st.mask_size);


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

        % -- Initialize trial variables, timing table ---------------------

        % Init time stamps (will be overwritten on first Screen flip)
        trial_timestamps = nan(4, trial_frames);
        
        % Vectors for touch and gaze data
        trial_touch_samples = nan(3, trial_frames); % x, y, is_touching
        trial_gaze_samples = nan(2, trial_frames);

        % Record of which stims were displayed at each frame
        stim_presentations = zeros(num_stims, trial_frames);
        
        % Response time - depends on fixation type, see README 
        rt_start_target_onset = NaN;
        
        rt_target_to_choice = NaN;
        
        rt_start_saccade = NaN;
        rt_target_to_saccade = NaN;
        rt_saccade_to_choice = NaN;
        
        rt_start_fingerlift = NaN;
        rt_target_to_fingerlift = NaN;
        rt_fingerlift_to_choice = NaN;
        
        % State variables which control response time markers
        [start_saccade, start_fingerlift, start_target_onset] = deal(false);
%         
%         % State variables which control response time markers
%         [start_rt1, end_rt1, start_rt2, end_rt2] = deal(false);
        
        % Cue to determine whether a response has been made
        hasSelected = false;
        
        % Cue to determine if a fixation pause is happening
        isFixation = false;
        
        % Frame counters: one for total elapsed time (incl. fixation pause)
        % and one for stimuli schedule
        frame = 1;
        stim_sched_frame = 1;
        
        % If we have passed fixation requirement, save this state as we
        % will begin reactime time on saccade/finger lift
        post_fix = 0;
        
        % -- Display trial stimuli, record responses ----------------------

        % Flip again to sync us to the vertical retrace
        vbl = Screen('Flip', window);
        
        % Display all stimuli which are scheduled for the current frame
        while frame < trial_frames + 1

            % Listen for Esc key to abort experiment
            [keyIsDown, ~, keyCode] = KbCheck;
            if keyIsDown && keyCode(ptb.escapeKey)
                sca;
                return
            end

            % Get the current position of the mouse
            [tx, ty, buttons] = GetMouse(window);
            is_touch = any(buttons ~= 0);

            % We clamp the values at the maximum values of the screen
            % in X and Y in case people have two monitors connected.
            tx = min(tx, ptb.screenXpixels);
            ty = min(ty, ptb.screenYpixels);

            % Store touch location and whether finger is touching
            trial_touch_samples(1, frame) = tx;
            trial_touch_samples(2, frame) = ty;
            trial_touch_samples(3, frame) = is_touch;

            % Get current gaze value
            % TBD
            
            % Test if saccade has occurred, start response period 2
            if post_fix == 1
                start_saccade = true;
                post_fix = 0;
            end    
               
            % Test if finger lift has occurred, start response period 2
            if post_fix == 2 && ~is_touch
                start_fingerlift = true;
                post_fix = 0;
            end
            
            % Check if fixation appears in current frame
            fixations = find(fixation_schedule(:,frame) ~= 0, 1);
            if ~isempty(fixations)
                isFixation = true;
                post_fix = 0;
                for f=1:length(fixations)
                    fidx = fixations(f);
                    % Set fixation to active
                    fixation_counter(fidx,3) = true;
                    % Note fixation type in counter
                    fixation_counter(fidx,4) = fixation_schedule(fidx,frame);
                end
            end
            
            % Test each active fixation for gaze/touch and iterate counters       
            active_fixations = find(fixation_counter(:,3) == true, 1);
            if isempty(active_fixations)
            	isFixation = false;
                post_fix = trial_dat.stim_is_target(targetMask);
            else
                for f=1:length(active_fixations)
                    fidx = active_fixations(f);
                    
                    % Render fixation image
                    fix_name = char(trial_dat.stim_img_name(fidx));
                    Screen('DrawTextures', window,...
                            stim_textures(fix_name),[],...
                            stim_bounds(fidx,:), ...
                            trial_dat.stim_rotation(fidx));
                    stim_presentations(fidx, frame) = 1;
                        
                    % If fixation is gaze type, check for gaze, reset if
                    % gaze lost
                    if fixation_counter(fidx,4) == 1
                    end
                    
                    % If fixation is touch type check for touch, reset if
                    % touch lifted
                    if fixation_counter(fidx,4) == 2
                        if IsInRect(tx, ty, dotBoundingRect{fidx}) && is_touch
                            fixation_counter(fidx,2) = fixation_counter(fidx,2) + 1;
                        else
                            fixation_counter(fidx,2) = 0;
                        end % touched fixation mask
                    end
                    
                    % If subject fixation duration exceeds requirement,
                    % deactivate fixation
                    if fixation_counter(fidx,1) < fixation_counter(fidx,2)
                        fixation_counter(fidx,3) = false;
                    end
                end
            end
            
            % Check image and mask schedules and render anything to screen
            % that is scheduled to be shown on the current frame
            for s=1:num_stims
                if stim_schedule(s, stim_sched_frame) == 1
                    stim_name = char(trial_dat.stim_img_name(s));
                    
                    % Display stim
                    Screen('DrawTextures', window,...
                            stim_textures(stim_name),[],...
                            stim_bounds(s,:), trial_dat.stim_rotation(s));
                    
                    % Record stim presentation in global frame counter 
                    stim_presentations(s, frame) = 1;
                    
                    % If stim is target, begin reaction time counter
                    if s == targetMask
                        start_target_onset = true;
                    end
                end % image scheduled
                if mask_schedule(s, stim_sched_frame) == 1
                    % Draw mask dots
                    Screen('DrawDots', window, dotPosMatrix{s},...
                            dotSizes{s}, dotColors{s}, [], 2);
                    
                    % Record mask presentation in global frame counter
                    % Stim only = 1, mask only = 2, stim & mask = 3
                    stim_presentations(s, frame) = stim_presentations(s, frame) + 2;
                        
                    % Check for collisions.
                    if trial_dat.mask_touchable(s) == 1
                        if IsInRect(tx, ty, dotBoundingRect{s}) && is_touch
                            hasSelected = true;
                            selectedMask = s;
                        end % touched stim mask
                    end % stim mask touchable
                end % mask scheduled
            end % stim loop

            % Flip to the screen
            [vbl, sot, flip, miss, ~] = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
            
            % Store exact timestamp and max error for current frame
            trial_timestamps(:, frame) = [vbl, sot, flip, miss];
            
            % If a touchable mask has been touched, store response time
            if hasSelected == true
                rt_target_to_choice = frame - rt_start_target_onset;
                rt_saccade_to_choice = frame - rt_start_saccade;
                rt_fingerlift_to_choice = frame - rt_start_fingerlift;
            end

            % Note start time to calculate reaction time 1
            if start_target_onset
                rt_start_target_onset = frame;
                start_target_onset = false;
            end
            
            % If saccade occurs, note time & start rt counter
            if start_saccade
                rt_start_saccade = frame;
                rt_target_to_saccade = frame - rt_start_target_onset;
                start_saccade = false;
            end
            
            % If finger lift occurs, note time & start rt counter
            if start_fingerlift
                rt_start_fingerlift = frame;
                rt_target_to_fingerlift = frame - rt_start_target_onset;
                start_fingerlift = false;
            end
            
            % End trial if selection has been made
            if hasSelected == true
                correct_choice = (selectedMask == targetMask);
                break
            end
                
            % If not in a fixation pause, iterate stim schedule
            if ~isFixation
                stim_sched_frame = stim_sched_frame + 1;
            end
            
            % Iterate trial frame counter
            frame = frame + 1;
            
        end % frame loop
        
        % Save total trial time
        trial_frames_elapsed = min(frame, trial_frames);
        
        % If trial timed out, display time out message
        if trial_frames_elapsed >= trial_frames
            % We clear the screen once they have made their response
            DrawFormattedText(window, timeoutMsg, ...
                                'center', 'center', ptb.black);
            Screen('Flip', window);
            KbStrokeWait;
        end
        
        % Output frame-by-frame timestamp/presentation record
        timestamp_record = [1:trial_frames; trial_timestamps; stim_presentations; trial_touch_samples; trial_gaze_samples]';
        timestamp_labels = {'frame_number','timestamp',...
                            'stimulus_onset_time',...
                            'flip_timestamp','missed'};
        stim_ons = {};
        for i=1:num_stims
            stim_ons = [stim_ons {['stim_on_' num2str(i) '_']}];
        end
        stim_sched_labels = strcat(stim_ons, trial_dat.stim_img_name');
        touch_labels = {'touch_loc_x', 'touch_loc_y', 'is_touching'};
        gaze_labels = {'gaze_loc_x', 'gaze_loc_y'};
        timerecord_out_labels = [timestamp_labels stim_sched_labels touch_labels gaze_labels];
        timestamp_out = array2table(timestamp_record, 'VariableNames', timerecord_out_labels);
        
        ts_out_name = [out_path subj '_timestamp_schedule_' 'b' num2str(bn, '%02d') '_t' num2str(tn, '%02d') '.csv'];
        writetable(timestamp_out, ts_out_name, 'Delimiter', ',', 'QuoteStrings', true);
        
        % -- Output data for current trial --------------------------------

        % Record the trial data into out data matrix
        subjData.Subject{row_ct} = subj;
        subjData.BlockNum{row_ct} = bn;
        subjData.TrialNum{row_ct} = row_ct;
        subjData.BadTrial{row_ct} = 0; % TBD
        
        % Extract target stimulus info from trial stim table if possible
        if ~isnan(targetMask)
            tr = trial_dat(targetMask,:);
            subjData.target_mask_name{row_ct} = cell2mat(tr.stim_img_name);
            subjData.target_center_x{row_ct} = tr.stim_cent_x;
            subjData.target_center_y{row_ct} = tr.stim_cent_y;
            subjData.target_rotation{row_ct} = tr.stim_rotation;
            subjData.target_mask_num{row_ct} = targetMask;
        end
        
        % Extract chosen stimulus info from trial stim table if possible
        if ~isnan(selectedMask)
            ch = trial_dat(selectedMask,:);
            subjData.chosen_mask_name{row_ct} = cell2mat(ch.stim_img_name);
            subjData.chosen_center_x{row_ct} = ch.stim_cent_x;
            subjData.chosen_center_y{row_ct} = ch.stim_cent_y;
            subjData.chosen_rotation{row_ct} = ch.stim_rotation;
            subjData.chosen_mask_num{row_ct} = selectedMask;
        end
        
        subjData.choice_correct{row_ct} = correct_choice;
        
        subjData.touch_point_x{row_ct} = tx;
        subjData.touch_point_y{row_ct} = ty;
        subjData.trial_time_elapsed{row_ct} = trial_frames_elapsed * ifi;
        subjData.rt_target_to_choice{row_ct} = rt_target_to_choice * ifi;
        subjData.rt_target_to_saccade{row_ct} = rt_target_to_saccade * ifi;
        subjData.rt_saccade_to_choice{row_ct} = rt_saccade_to_choice * ifi;
        subjData.rt_target_to_fingerlift{row_ct} = rt_target_to_fingerlift * ifi;
        subjData.rt_fingerlift_to_choice{row_ct} = rt_fingerlift_to_choice * ifi;
        
        % -- Write trial data out to one line of subject data file --------
        % Extract row corresponding to current trial
        trial_data = subjData{row_ct,:};
        
        % Fill any empty entries with NaN
        for idx = 1:numel(trial_data)
            if isempty(trial_data{idx})
                trial_data{idx} = NaN;
            end
        end
        
        % Convert all cells to strings
        trial_data = cellfun(@num2str,trial_data,'UniformOutput',false);
        
        % Append strings
        trial_row = strjoin(trial_data,',');
        
        % Write line to subject data output file
        fid = fopen(subj_data_fname,'at');
        if fid>0
            fprintf(fid,'%s,',trial_row);
            fprintf(fid,'\n');
            fclose(fid);
        end
        
        % Iterate row counter for subject data
        row_ct = row_ct + 1;
    end % trial loop

    % Flip again to sync us to the vertical retrace
    Screen('Flip', window);

    if bn < num_blocks
        % End of block screen. 
        % We clear the screen once they have made their response
        DrawFormattedText(window, blockEndMsg, ...
                            'center', 'center', ptb.black);
        Screen('Flip', window);
        [~, keyCode, ~] = KbStrokeWait;
        if keyCode(ptb.escapeKey)
            sca;
            return;
        end
    end

end % block loop

% Flip again to sync us to the vertical retrace
Screen('Flip', window);

% End of experiment screen. We clear the screen once they have made their
% response
DrawFormattedText(window, expEndMsg, 'center', 'center', ptb.black);
Screen('Flip', window);
KbStrokeWait;
sca;