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

function ShapeTapper_RunExperiment(config_fname)

% Clear the workspace and the screen
sca;
close all;
clearvars;

% Add source code path
addpath('source');

% Seed random number generator for participant ID
rand('seed', sum(100 * clock));

% Set messages to participant
welcomeMsg = ['Welcome to the experiment!\n\n\n',...
              'Please wait for instructions \n',...
              'from the experimenter to continue.\n\n\n'...
              'Press Any Key To Begin\n',...
              'Esc to exit at any time'];
          
practiceStartMsg = ['Now beginning practice block.\n\n\n',...
                    'Press Any Key To Continue'];
practiceEndMsg = 'Practice Block Finished\n\n\nPress Any Key To Continue';
          
blockEndMsg = 'Block Finished\n\n\nPress Any Key To Continue';

expEndMsg = 'Experiment Finished\n\n\nPress Any Key To Exit';

timeoutMsg = 'Trial Time Expired.\n\n\nPress Any Key To Continue';

feedback_message_correct = 'Correct!\n\n\nPress Any Key To Continue';
feedback_message_incorrect = 'Incorrect!\n\n\nPress Any Key To Continue';

% Acceptable stimulus image formats, must be compatible with imread()
img_formats = {'.png', '.jpg'};

% Specify directory containing stimulus images
stim_dir = 'Image_Files\';

% Data output path
out_path = 'Data/';

% Safety margin added to touchable area around stimulus images 
% (touchable area is the defined as the smallest circle containing image)
TOUCH_MARGIN = 1;

%----------------------------------------------------------------------
%                       Experimenter Setup
%----------------------------------------------------------------------
% Get demographics info
part_dems = GetDemographics();

% Specify directory containing config files
switch nargin
   case 0
       [config_fname,config_path] = uigetfile('./Config_Files/*.csv',...
                                      'Select an experiment %%config file');
       config_fname = [config_path config_fname];
   case 1
       if isempty(config_fname) || isnan(config_fname)
           config_fname = 'Config_Files\test_config.csv';
       end
end

% Extract config file name for demographics file
[~, fname, ext] = fileparts(config_fname); 
part_dems.config_name = [fname ext];

% Save date for demographics file
part_dems.experiment_date = char(datetime);

% Read specified config file for this experiment
config_dat = readtable(config_fname);

% Read demographics file to check if participant ID has already bene used
demog_fname = 'Data\Demographics_ShapeTapper.csv';
if ~exist(demog_fname, 'file')
    demog_header = 'ParticipantID,Handedness,Gender,Age,Config_File,Experiment_Time';
    fid = fopen(demog_fname,'at');
    if fid>0
       fprintf(fid,'%s,',demog_header);
       fprintf(fid,'\n');
       fclose(fid);
    end
end

demog_dat = readtable(demog_fname);
matches_participant = strcmp(part_dems.id, demog_dat.ParticipantID);
matches_config =  strcmp(part_dems.config_name, demog_dat.Config_File);
if any(all([matches_participant matches_config],2))
    match_idx = find(all([matches_participant matches_config],2));
    
    % Include the desired Default answer
    opts.Interpreter = 'tex';
    opts.Default = 'Cancel';
    part_message = ['Participant \bf' part_dems.id ...
                   '\rm has existing data for config \bf' ...
                   strrep(part_dems.config_name,'_','\_') ...
                   '\rm \newline \newline' ...
                   'Would you like to continue this experiment, or start over?'];
    answer = questdlg(part_message, 'Participant ID Already Exists!', ...
                      'Continue', 'Start Over', 'Cancel', opts);
    switch answer
        case 'Continue'
            old_data_fname = [out_path part_dems.id '/' fname '/' part_dems.id '_results_out.csv'];
            old_data = readtable(old_data_fname,'Header', 1, 'Delimiter',',');
            last_block = old_data{end, 2};
            last_trial = old_data{end, 3};
            opts2.Interpreter = 'tex';
            opts2.Default = 'Cancel';
            resume_message = ['Resuming experiment with config file ' ...
                         strrep(part_dems.config_name,'_','\_') ...
                          '\rm \newline' ...
                          'Start from block ' num2str(last_block) ...
                          ', trial ' num2str(last_trial) ' ?'];
            resume_answer = questdlg(resume_message, 'Resume Experiment', ...
                             'Resume', 'Cancel', opts2);
            switch resume_answer
                case 'Resume'
                    dat_block_rows = config_dat.block_num == last_block;
                    dat_trial_rows = config_dat.trial_num == last_trial;
                    resume_row = find(all([dat_block_rows dat_trial_rows],2), 1, 'last');
                    config_dat = config_dat(resume_row+1:end,:);
                case 'Cancel'
                    return
            end
        case 'Start Over'
            old_path = [out_path part_dems.id '/' fname];
            movefile([old_path '/'],[old_path '_old/']);
        case 'Cancel'
            return
    end
end
%----------------------------------------------------------------------
%                       PyschToolbox Setup
%----------------------------------------------------------------------
% Get background color from config file
bg_color = config_dat.background_color(1);

% Run function to init PsychToolbox screen and store relevant vars for
% later rendering and timing
ptb = ptb_initscreen( bg_color );

% Get window handle for PTB rendering
window = ptb.window;

% Get monitors Inter-Frame Interval and waittime (# of frames b/t renders)
ifi = ptb.ifi;
waitframes = ptb.waitframes;

% Flip to clear
Screen('Flip', window);

%----------------------------------------------------------------------
%                     Load Stimuli Images
%----------------------------------------------------------------------

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
% Output demographics information
%Define output name for demographics file
demog_fname = [out_path 'Demographics_ShapeTapper.csv'];

% Convert demographic info to cell array
demog_data = struct2cell(part_dems);

% Convert all cells to strings
part_demog_data = cellfun(@num2str,demog_data,'UniformOutput',false);

% Append strings
demographic_row = strjoin(part_demog_data,',');

% Write line to participant data output file
fid = fopen(demog_fname,'at');
if fid>0
    fprintf(fid,'%s,',demographic_row);
    fprintf(fid,'\n');
    fclose(fid);
end

% Participant data output path
part_out_path = [out_path part_dems.id '/' fname '/'];
mkdir(part_out_path);

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
data_columns = {'participant',...
                'BlockNum',...
                'TrialNum',...
                'BadTrial',...
                'target_image_name',...
                'target_center_x',...
                'target_center_y',...
                'target_size_x',...
                'target_size_y',...
                'target_rotation',...
                'chosen_image_name',... 
                'chosen_center_x',...
                'chosen_center_y',...
                'chosen_size_x',...
                'chosen_size_y',...
                'chosen_rotation',...
                'target_image_num',...
                'chosen_image_num',...
                'choice_correct',...
                'touch_point_x',...
                'touch_point_y',...
                'touch_point_relative_x',...
                'touch_point_relative_y',...
                'trial_time_elapsed',...
                'rt_target_to_choice',...
                'rt_target_to_saccade',...
                'rt_saccade_to_choice',...
                'rt_target_to_fingerlift',...
                'rt_fingerlift_to_choice'};

partData = cell2table(cell(tot_num_trials, length(data_columns)),...
                      'VariableNames',data_columns);

% Counter for correct row indexing when writing output
row_ct = 1;

% Set output file name
part_data_fname = [part_out_path part_dems.id '_results_out.csv'];
if ~exist(part_data_fname, 'file')
    fid = fopen( part_data_fname, 'wt+' );
    if fid>0
        fprintf(fid,strjoin(data_columns,','));
        fprintf(fid,'\n');
        fclose(fid);
    end
end

%----------------------------------------------------------------------
%                       Experimental loop
%----------------------------------------------------------------------

abort_experiment = false;

% Animation Loop -- here we render and display all blocks & trials
% specified in the config file, and save results

% BLOCK LOOP
for b=1:num_blocks
    
    % If this is the very  first trial, present a start screen and wait
    % for a key-press
    if b == 1
        % Set text and background color to match first trial
        welcome_bg = ptb.white;
        welcome_text = ptb.black;
        stim_bg_color = cell2mat(config_dat.background_color(1));
        if ischar(stim_bg_color)
            if strcmp(stim_bg_color,'black')
                welcome_bg = ptb.black;
                welcome_text = ptb.white;
            end
        end
        Screen('FillRect', window, welcome_bg);
        
        % Draw welcome text
        DrawFormattedText(window,welcomeMsg,'center','center', welcome_text);
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
    
    % Display practice block message, if applicable
    if block_dat.is_practice_block(1)
        block_bg = ptb.white;
        block_text = ptb.black;
        stim_bg_color = cell2mat(block_dat.background_color(1));
        if ischar(stim_bg_color)
            if strcmp(stim_bg_color,'black')
                block_bg = ptb.black;
                block_text = ptb.white;
            end
        end
        Screen('FillRect', window, block_bg);

        % Draw practice block or generic block end message
        DrawFormattedText(window, practiceStartMsg, ...
                        'center', 'center', block_text);
                    
        Screen('Flip', window);
        
        % Listen for Esc key to abort experiment
        [~, keyCode, ~] = KbStrokeWait;
        if keyCode(ptb.escapeKey)
            sca;
            return;
        end
    end
    
    % Create vector to flag bad trials for repeat at end of block
    bad_trials = zeros(num_trials,1);
    
    % Set trial counter
    t = 1;
    block_complete = false;
    rpt_bad_trials = false;
    
    % TRIAL LOOP
    while ~block_complete

        % If all trials have been presented
        if t > num_trials
           % If no bad trials, then complete block
           if isempty(find(bad_trials))
               rpt_bad_trials = false;
               block_complete = true;
               continue
           else
               % Begin looping through bad trials
               rpt_bad_trials = true;
           end
        end
        
        % If we are repeating bad trials, then randomly present a bad trial
        if rpt_bad_trials
            bad_trials_remaining = find(bad_trials);
            % If no bad trials, then complete block
            if isempty(bad_trials_remaining)
               rpt_bad_trials = false;
               block_complete = true;
               continue
            else
                rand_bad_trial = randi([1 length(bad_trials_remaining)]);
                t = bad_trials_remaining(rand_bad_trial);
            end   
        end
        
        if abort_experiment
            break
        end
        
        % -- Pre-trial initialization: event schedules, calculate masks --

        % Collect all trial numbers
        tn = trial_list(t);
        
        % Retrieve stimuli list for current trial
        trial_dat = block_dat(block_dat.trial_num == tn, :);
        
        % Calculate total number of frames for current trial
        trial_frames  = ceil(max(trial_dat.trial_max_time) / ptb.ifi);

        % Extract array of stim names
        stim_img_names = trial_dat.stim_img_name;
        
        % Make event schedule arrays, which will allow for rapid checking
        % of which stimuli to render for each frame
        num_stims = size(stim_img_names, 1);
        stim_schedule = zeros(num_stims, trial_frames);
        mask_schedule = zeros(num_stims, trial_frames);
        fixation_schedule = zeros(num_stims, trial_frames);
        
        % Initialize counters for fixation durations
        fixation_counter = zeros(num_stims, 3);

        % Initialize record of stim displayed
        stim_displayed = zeros(num_stims, 1);
        
        % Initialize bounding rectangle for each stimulus
        stim_bounds = zeros(num_stims, 4);
        stim_centers = zeros(num_stims, 2);
        stim_radius = zeros(num_stims, 1);

        % Initialize mask arrays for each stimlus
        dotPosMatrix = cell(1,num_stims);
        dotSizes = cell(1,num_stims);
        dotColors = cell(1,num_stims);
        dotBoundingRect = cell(1,num_stims);

        % Declare target mask
        targetMask = NaN;
        selectedStim = NaN;
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
            stim_cx_px = st.stim_cent_x * ptb.screenXpixels;
            stim_cy_px = st.stim_cent_y * ptb.screenYpixels;
            rect = CenterRectOnPoint(rect, stim_cx_px, stim_cy_px);

            % Store rectangle in stim_bounds array for texture rendering
            stim_bounds(s,:) = rect;

            % Calculate diagonal and center in pixels for touch detection
            stim_size_x_px = st.stim_size_x * ptb.ppcm;
            stim_size_y_px = st.stim_size_y * ptb.ppcm;
            stim_centers(s,:) = [stim_cx_px stim_cy_px];
            stim_radius(s) = sqrt(stim_size_x_px^2 + stim_size_y_px^2)/2 + ...
                                    (TOUCH_MARGIN*ptb.ppcm);
            
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
            mask_pos = mask_pos - [repmat(stim_cx_px,1,c); repmat(stim_cy_px,1,c)];

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
            mask_pos = mask_pos + [repmat(stim_cx_px,1,c); repmat(stim_cy_px,1,c)];

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
        
        % Cue to determine whether a response has been made
        hasSelected = false;
        
        % Cue to determine if a fixation pause is happening
        isFixation = false;
        
        % Frame counters: one for total elapsed time (incl. fixation pause)
        % and one for stimuli schedule
        frame = 1;
        stim_sched_frame = 1;
        
        % This state is true prior to target onset, false after (used for
        % reaction time calculation)
        pre_target = 1;
        
        % If we have passed fixation requirement, save this state as we
        % will begin reactime time on saccade/finger lift
        post_fix = 0;
        
        % Detect starting touch/gaze and require lift/saccade to begin
        % recording data
        touch_on_start = 0;
        gaze_on_start = 0;
        
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
            
            % Note when a background color change has occurred in the
            % current frame
            bg_color_change = false;

            % Get the current position of the mouse
            [tx, ty, buttons] = GetMouse(window);
            is_touch = any(buttons ~= 0);

            % Detect touch/gaze on trial start
            if frame == 1 && is_touch
                touch_on_start = true;
            end
            
            % If trial started with touch, disable is_touch and require
            % a finger lift to register a new touch
            if touch_on_start
                if is_touch
                    is_touch = false;
                else
                    touch_on_start = false;
                end
            end

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
            	if isFixation
                    isFixation = false;
                    post_fix = 1;
                end
            else
                for f=1:length(active_fixations)
                    fidx = active_fixations(f);
                    
                    % Detect if backgroudnd color changed
                    if ~bg_color_change
                        stim_bg = ptb.white;
                        stim_bg_color = cell2mat(trial_dat.background_color(fidx));
                        if ischar(stim_bg_color)
                            if strcmp(stim_bg_color,'black')
                                stim_bg = ptb.black;
                            end
                        end
                        Screen('FillRect', window, stim_bg);
                        bg_color_change = true;
                    end
                    
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
                        if IsInRect(tx, ty, dotBoundingRect{fidx})% && is_touch
                            fixation_counter(fidx,2) = fixation_counter(fidx,2) + 1;
                        else
                            fixation_counter(fidx,2) = 0;
                        end % touched fixation mask
                    end
                    
                    % If participant fixation duration exceeds requirement,
                    % deactivate fixation
                    if fixation_counter(fidx,1) < fixation_counter(fidx,2)
                        fixation_counter(fidx,3) = false;
                    end
                end
            end
            
            
            % Check image and mask schedules and render anything to screen
            % that is scheduled to be shown on the current frame
            if ~isFixation
                for s=1:num_stims
                    if stim_schedule(s, stim_sched_frame) == 1
                        stim_name = char(trial_dat.stim_img_name(s));

                        % Detect if backgroudnd color changed
                        if ~bg_color_change
                            stim_bg = ptb.white;
                            stim_bg_color = cell2mat(trial_dat.background_color(s));
                            if ischar(stim_bg_color)
                                if strcmp(stim_bg_color,'black')
                                    stim_bg = ptb.black;
                                end
                            end
                            Screen('FillRect', window, stim_bg);
                            bg_color_change = true;
                        end

                        % Display stim
                        Screen('DrawTextures', window,...
                                stim_textures(stim_name),[],...
                                stim_bounds(s,:), trial_dat.stim_rotation(s));

                        % Record stim presentation in global frame counter 
                        stim_presentations(s, frame) = 1;
                        stim_displayed(s) = 1;
                        
                        % If stim is target, begin reaction time counter
                        if s == targetMask && pre_target
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
                    end % mask scheduled
                end % stim loop
                % Check for collisions
                for s=1:num_stims
                    if stim_displayed(s) && trial_dat.stim_is_touchable(s) == 1
                        if ptb_in_circ(stim_centers(s,:), [tx ty], stim_radius(s)) && is_touch
                            hasSelected = true;
                            selectedStim = s;
                            if bad_trials(t)
                                bad_trials(t) = false;
                            end
                        end% touched stim
                    end % image displayed & touchable
                end
                % Check for any non-stim touches, add trial to lsit of bad
                % trials to repeat at end of block
                if any(stim_displayed) && ~hasSelected && is_touch
                    hasSelected = true;
                    selectedStim = NaN;
                    bad_trials(t) = true;
                end
            end % if fixation pause
                
            % Flip to the screen
            [vbl, sot, flip, miss, ~] = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                        
            % Store exact timestamp and max error for current frame
            trial_timestamps(:, frame) = [vbl, sot, flip, miss];
            
            % If missed frame, iterate 2 frames to stay on track.
            missed_frame = (miss > 0);
            next_frame = 1 + missed_frame;
            
            % If missed frame, copy fixation events to subsequent frame so
            % as not to skip them.
            if missed_frame && (frame + 1 < trial_frames)
                fixation_schedule(:,frame+2) = fixation_schedule(:,frame+1);
            end
          
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
                pre_target = false;
            end
            
            % If a fixation has occured, record finger lift or saccade
            % reaction time
            if post_fix
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
            end % if post_fix
            
            % End trial if selection has been made
            if hasSelected == true
                correct_choice = (selectedStim == targetMask);
                break
            end
                
            % If not in a fixation pause, iterate stim schedule
            if ~isFixation
                stim_sched_frame = stim_sched_frame + next_frame;
            end
            
            % Iterate trial frame counter
            frame = frame + next_frame;
            
        end % frame loop
        
        % Save total trial time
        trial_frames_elapsed = min(frame, trial_frames);
        
        % If trial timed out, display time out message
        if trial_frames_elapsed >= trial_frames
        	% Set text and background color to match first trial
            timeout_bg = ptb.white;
            timeout_text = ptb.black;
            stim_bg_color = cell2mat(trial_dat.background_color(num_stims));
            if ischar(stim_bg_color)
                if strcmp(stim_bg_color,'black')
                    timeout_bg = ptb.black;
                    timeout_text = ptb.white;
                end
            end
            Screen('FillRect', window, timeout_bg);

            % Draw timeout text
            Screen('Flip', window);
            DrawFormattedText(window, timeoutMsg, ...
                                'center', 'center', timeout_text);
            Screen('Flip', window);
            [~, keyCode, ~] = KbStrokeWait;
            if keyCode(ptb.escapeKey)
                abort_experiment = true;
            end
        elseif block_dat.is_practice_block(1) && trial_dat.trial_feedback(end)
            % if in practice block, with feedback, display feedback
            % Set text and background color to match first trial
            feedback_bg = ptb.white;
            timeout_text = ptb.black;
            stim_bg_color = cell2mat(trial_dat.background_color(num_stims));
            if ischar(stim_bg_color)
                if strcmp(stim_bg_color,'black')
                    feedback_bg = ptb.black;
                    feedback_text = ptb.white;
                end
            end
            
            if correct_choice
                feedbackMsg = feedback_message_correct;
            else
                feedbackMsg = feedback_message_incorrect;
            end
            
            Screen('FillRect', window, feedback_bg);

            % Draw timeout text
            Screen('Flip', window);
            DrawFormattedText(window, feedbackMsg, ...
                                'center', 'center', feedback_text);
            Screen('Flip', window);
            [~, keyCode, ~] = KbStrokeWait;
            if keyCode(ptb.escapeKey)
                abort_experiment = true;
            end
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
        
        ts_out_name = [part_out_path part_dems.id '_timestamp_schedule_' 'b' num2str(bn, '%02d') '_t' num2str(tn, '%02d') '.csv'];
        writetable(timestamp_out, ts_out_name, 'Delimiter', ',', 'QuoteStrings', true);
        
        % -- Output data for current trial --------------------------------

        % Record the trial data into out data matrix
        partData.participant{row_ct} = part_dems.id;
        partData.BlockNum{row_ct} = bn;
        partData.TrialNum{row_ct} = tn;
        partData.BadTrial{row_ct} = 0; % TBD
        
        % Extract target stimulus info from trial stim table if possible
        if ~isnan(targetMask)
            tr = trial_dat(targetMask,:);
            partData.target_image_name{row_ct} = cell2mat(tr.stim_img_name);
            partData.target_center_x{row_ct} = tr.stim_cent_x*ptb.screenXpixels;
            partData.target_center_y{row_ct} = tr.stim_cent_y*ptb.screenYpixels;
            partData.target_size_x{row_ct} = tr.stim_size_x*ptb.ppcm;
            partData.target_size_y{row_ct} = tr.stim_size_y*ptb.ppcm;
            partData.target_rotation{row_ct} = tr.stim_rotation;
            partData.target_image_num{row_ct} = targetMask;
        end
        
        % Extract chosen stimulus info from trial stim table if possible
        if ~isnan(selectedStim)
            %Extract chosen stimulus info
            ch = trial_dat(selectedStim,:);
            
            % Convert center to pixels
            ch_cent_x_px = ch.stim_cent_x*ptb.screenXpixels;
            ch_cent_y_px = ch.stim_cent_y*ptb.screenYpixels;
            
            % Get chosen stim rotation
            c_rot = ch.stim_rotation;
            
            % Save necessary stim info
            partData.chosen_image_name{row_ct} = cell2mat(ch.stim_img_name);
            partData.chosen_center_x{row_ct} = ch_cent_x_px;
            partData.chosen_center_y{row_ct} = ch_cent_y_px;
            partData.chosen_size_x{row_ct} = ch.stim_size_x*ptb.ppcm;
            partData.chosen_size_y{row_ct} = ch.stim_size_y*ptb.ppcm;
            partData.chosen_rotation{row_ct} = c_rot;
            partData.chosen_image_num{row_ct} = selectedStim;
            
            % -- Calculate relative touchpoint to chosen object --
            % Remove offset by chosen center
            tp_x_offset = tx - ch_cent_x_px;
            tp_y_offset = ty - ch_cent_y_px;
            % Create rotation matrix (reverse rotation of chosen stim)
            R = [cosd(-c_rot) -sind(-c_rot); sind(-c_rot) cosd(-c_rot)];
            tp_cent = [tp_x_offset; tp_y_offset];
            % Reverse stim rotation to get touch point relative to center
            tp_cent_rot = R*tp_cent;
            % Offset by stim size (in pixels) and flip y-axis to get
            % top-left = (0,0) coordinate system used in OpenCV analysis
            tp_rel_x = tp_cent_rot(1) + (ch.stim_size_x*ptb.ppcm/2);
            tp_rel_y = -(tp_cent_rot(2) - (ch.stim_size_x*ptb.ppcm/2));
            % Save relative touchpoint in participant trial data
        	partData.touch_point_relative_x{row_ct} = tp_rel_x;
            partData.touch_point_relative_y{row_ct} = tp_rel_y;
        else
            % Non-stim touch, mark as bad trial
            partData.BadTrial{row_ct} = 1;
        end
        
        partData.choice_correct{row_ct} = correct_choice;
        partData.touch_point_x{row_ct} = tx;
        partData.touch_point_y{row_ct} = ty;
        partData.trial_time_elapsed{row_ct} = (trial_frames_elapsed+1) * ifi;
        partData.rt_target_to_choice{row_ct} = rt_target_to_choice * ifi;
        partData.rt_target_to_saccade{row_ct} = rt_target_to_saccade * ifi;
        partData.rt_saccade_to_choice{row_ct} = rt_saccade_to_choice * ifi;
        partData.rt_target_to_fingerlift{row_ct} = rt_target_to_fingerlift * ifi;
        partData.rt_fingerlift_to_choice{row_ct} = rt_fingerlift_to_choice * ifi;
        
        % -- Write trial data out to one line of participant data file --------
        % Extract row corresponding to current trial
        trial_data = partData{row_ct,:};
        
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
        
        % Write line to participant data output file
        fid = fopen(part_data_fname,'at');
        if fid>0
            fprintf(fid,'%s,',trial_row);
            fprintf(fid,'\n');
            fclose(fid);
        end
        
        % Iterate row counter for participant data
        row_ct = row_ct + 1;
        
        % Iterate trial counter
        t = t + 1;
    end % trial loop

    % Flip again to sync us to the vertical retrace
    Screen('Flip', window);
    
    if ~abort_experiment
        if bn < num_blocks
            % End of block screen. 
            % Set text and background color to match first trial
            block_bg = ptb.white;
            block_text = ptb.black;
            stim_bg_color = cell2mat(trial_dat.background_color(num_stims));
            if ischar(stim_bg_color)
                if strcmp(stim_bg_color,'black')
                    block_bg = ptb.black;
                    block_text = ptb.white;
                end
            end
            Screen('FillRect', window, block_bg);

            % Draw practice block or generic block end message
            if trial_dat.is_practice_block(1)
            	DrawFormattedText(window, practiceEndMsg, ...
                                'center', 'center', block_text);
            else
                % Draw block text
                DrawFormattedText(window, blockEndMsg, ...
                                'center', 'center', block_text);
            end
            Screen('Flip', window);
            [~, keyCode, ~] = KbStrokeWait;
            if keyCode(ptb.escapeKey)
                sca;
                return;
            end
        end % bn < num_blocks
    end % abort_experiment
    
end % block loop

% Flip again to sync us to the vertical retrace
Screen('Flip', window);

% End of experiment screen. We clear the screen once they have made their
% response
% Set text and background color to match first trial
eoe_bg = ptb.white;
eoe_text = ptb.black;
stim_bg_color = cell2mat(trial_dat.background_color(num_stims));
if ischar(stim_bg_color)
    if strcmp(stim_bg_color,'black')
        eoe_bg = ptb.black;
        eoe_text = ptb.white;
    end
end
Screen('FillRect', window, eoe_bg);
% Draw End of Experiment text
DrawFormattedText(window, expEndMsg, 'center', 'center', eoe_text);
Screen('Flip', window);
KbStrokeWait;
ShowCursor('arrow', window);
sca;