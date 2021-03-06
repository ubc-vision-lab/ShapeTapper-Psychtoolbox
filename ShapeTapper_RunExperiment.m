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

function ShapeSpider_RunExperiment(useEyelink, useOptotrak)

% Clear the workspace and the screen
sca;       
close all;
% clearvars;

instrreset
opto_initialized = 0; 

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
practiceEndMsg = ['Practice Block Finished\n\n\n',...
                  'Press Any Key To Continue'];
          
blockEndMsg = ['Block Finished\n\n\n',...
               'Press Any Key To Continue\n\n',...
               'Press F3 to Calibrate EyeLink'];

expEndMsg = 'Experiment Finished\n\n\nPress Any Key To Exit';

timeoutMsg = 'Trial Time Expired.\n\n\nPress Any Key To Continue';

feedback_message_correct = 'Correct!';
feedback_message_incorrect = 'Incorrect!';

gaze_message_lost_gaze = ['Lost gaze!\n\n\n',...
                          'Press Any Key To Continue\n\n',...
                          'Press F3 to Calibrate EyeLink'];

kbrespMsg = 'Letter?';

% Acceptable stimulus image formats, must be compatible with imread()
img_formats = {'.png', '.jpg'};

% Acceptable stimulus sound formats, must be compatible with wavread()
snd_formats = {'.wav'};

% Specify directory containing stimulus images
stim_dir = 'Image_Files/';

% Specify directory containing stimulus images
snd_dir = 'Sound_Files/';

% Data output path
out_path = 'Data/';

% Enter monitor diagonal size in inches
SCREEN_DIAG_INCHES = 31.5;

% Safety margin added to touchable area around stimulus images 
% (touchable area is the defined as the smallest circle containing image)
TOUCH_MARGIN_CM = 3;

% Amount of time (ms) gaze must settle on a target to register a choice
MIN_GAZE_CHOICE_TIME = 250;

% Amount of time (ms) we may lost gaze or fixation before aborting trial
GAZE_MAX_LOST_TIME = 500;
FIX_MAX_LOST_TIME = 50;
MAX_ON_FIX_RT = 1000;

%----------------------------------------------------------------------
%                       Experimenter Setup
%----------------------------------------------------------------------
% Get demographics info
part_dems = GetDemographics();

if isempty(fieldnames(part_dems))
   return 
end

[config_fname,config_path] = uigetfile('./Config_Files/*.csv',...
                              'Select an experiment %%config file');                     
config_fname = [config_path config_fname];

if ~ischar(config_fname)
    return
end

% Specify directory containing config files
switch nargin
   case 0
       useEyelink = 1;
       useOptotrak = 1;
%        [config_fname,config_path] = uigetfile('./Config_Files/*.csv',...
%                                       'Select an experiment %%config file');
%        config_fname = [config_path config_fname];
   case 2
%        if isempty(config_fname) || isnan(config_fname)
%            [config_fname,config_path] = uigetfile('./Config_Files/*.csv',...
%                                       'Select an experiment %%config file');
%            config_fname = [config_path config_fname];
%        end
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
% Get background and text color from config file
bg_color = config_dat.background_color(1);
text_color = config_dat.text_color(1);

% Run function to init PsychToolbox screen and store relevant vars for
% later rendering and timing
ptb = ptb_init_screen( bg_color, text_color, SCREEN_DIAG_INCHES );

% Get window handle for PTB rendering
window = ptb.window;

% Get monitors Inter-Frame Interval and waittime (# of frames b/t renders)
ifi = ptb.ifi;
waitframes = ptb.waitframes;

% Flip to clear
Screen('Flip', window);

%----------------------------------------------------------------------
%                       Sound Setup
%----------------------------------------------------------------------
useSound = 0;

if ismember('stim_sound_file', config_dat.Properties.VariableNames)
    if ~all(strcmp(config_dat.stim_sound_file,'None'))
%         InitializePsychSound(1); %init sound driver, 1 = low latency
        Snd('Open');
        
        % Retreive names of stimuli images
        has_sound_stim = ~(strcmp(config_dat.stim_sound_file,'None'));
        snd_names = unique(config_dat.stim_sound_file(has_sound_stim));

        % Extract unique strings for a list of stims
        snd_names = unique(snd_names);
        
        [snd_buffers, fs] = ptb_load_sounds(snd_dir, snd_names, snd_formats);
        
%         pa_handle = PsychPortAudio('Open', [], [], 2, fs, 1, 0); % opens sound buffer
        useSound = 1;
    end
end

%----------------------------------------------------------------------
%                       EyeLink Setup
%----------------------------------------------------------------------
if useEyelink
    % Initialize and calibrate Eyelink
    el = ptb_init_eyelink(window, out_path, fname, part_dems);

    % Test if Eyelink initialization was successful
    if ~isempty(el)
        useEyelink = 1;
        % Use right eye
        eye_used = el.RIGHT_EYE;
    else
        useEyelink = 0;
    end
    
    % Load calibration test screen
    calibration_sc = ptb_load_textures(stim_dir, {'cali_all'}, img_formats, ptb);
    
end 

% Set threshold for bad EyeLink samples in a given trial
eyelink_bad_count = 0; % count the # of times Eyelink has missed an eye

% Max tolerable number of consecutive frames gaze is lost
el.lost_gaze_max_frames = ceil(GAZE_MAX_LOST_TIME / (1000*ifi));
        
%----------------------------------------------------------------------
%                       Optotrak Setup
%----------------------------------------------------------------------
if useOptotrak
    max_trial_time = max(config_dat.trial_max_time);
    warning('off','all')
    % Initialize Optotrak
    optk = ptb_init_optotrak(opto_initialized, part_dems, out_path, max_trial_time);
    warning('on','all')
end 

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
stim_textures = ptb_load_textures(stim_dir, img_names, img_formats, ptb);

% Retreive names of feedback images (ones not marked 'None')
fb_img_idx =  ~strcmp(config_dat.trial_feedback_img,'None') ;
fb_img_names = unique(config_dat.trial_feedback_img(fb_img_idx));

% Extract unique strings for a list of imgs
fb_img_names = unique(fb_img_names);

fb_textures = ptb_load_textures(stim_dir, fb_img_names, img_formats, ptb);

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
                'keyboard_user_response',...
                'keyboard_corr_response',...
                'keyboard_correct',...
                'trial_time_elapsed',...
                'rt_target_to_choice',...
                'rt_target_to_saccade',...
                'rt_saccade_to_choice',...
                'rt_target_to_fingerlift',...
                'rt_fingerlift_to_choice'};

partData = cell2table(cell(tot_num_trials, length(data_columns)),...
                      'VariableNames',data_columns);

% Preallocate blank data row to expand data output if necessary
blank_data_row = array2table(nan(1,length(data_columns)),...
                        'VariableNames',partData.Properties.VariableNames);
                  
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

% Flag to detect experiment end (ESC)
abort_experiment = false;

% Max number of consecutive frames touch fixation is lost
lost_fix_max_frames = ceil(FIX_MAX_LOST_TIME / (1000*ifi));

% Max number of consective frames that user may remain on blank area after
% fixation is satsfied
on_fix_rt_max = ceil(MAX_ON_FIX_RT / (1000*ifi));

% Min number of consectutive frames to register gaze choice
gaze_choice_min_frames = ceil(MIN_GAZE_CHOICE_TIME / (1000*ifi));

pre_trial_calibration = false;

% Animation Loop -- here we render and display all blocks & trials
% specified in the config file, and save results

% BLOCK LOOP
for b=1:num_blocks
    
    % If this is the very first trial, present a start screen and wait
    % for a key-press
    if b == 1
        if useEyelink
            ptb_calibrate_eyelink( calibration_sc, el, ptb );
        end
        
        % Set text and background color to match first trial
        Screen('FillRect', window, ptb.bg_color);
        
        % Draw welcome text
        DrawFormattedText(window,welcomeMsg,'center','center', ptb.text_color);
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
        % Set text and background color to match first trial
        if ~bg_color_change
        	[ptb, bg_color_change] = ptb_set_bg_color(cell2mat(block_dat.background_color(1)), cell2mat(block_dat.text_color(1)), ptb);
        end
                    
        Screen('FillRect', window, ptb.bg_color);

        % Draw practice block or generic block end message
        DrawFormattedText(window, practiceStartMsg, ...
                            'center', 'center', ptb.text_color);
                    
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

        % If trial was aborted for re-calibration (press F3), then
        % recalibrate before starting next trial
        if pre_trial_calibration
            ptb_calibrate_eyelink( calibration_sc, el, ptb );
            pre_trial_calibration = false;
        end
        
        % If all trials have been presented
        if t > num_trials
           % If no bad trials, then complete block
           if isempty(find(bad_trials > 0, 1))
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
            bad_trials_remaining = find(bad_trials > 0);
            % If no bad trials, then complete block
            if isempty(bad_trials_remaining)
               rpt_bad_trials = false;
               block_complete = true;
               continue
            else
                rand_bad_trial = randi([1 length(bad_trials_remaining)]);
                t = bad_trials_remaining(rand_bad_trial);
                bad_trials(t) = 0;
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
        fixation_counter = zeros(num_stims, 5);

        % Initialize record of stim displayed
        stim_displayed = zeros(num_stims, 1);
        
        % Initialize record of sounds played
        snd_played = zeros(num_stims, 1);
        
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
        targetStims = find(trial_dat.stim_is_target ~= 0);
        selectedStim = NaN;
        correct_choice = false;
        
        % Populate stimulus schedules and mask vectors
        for s=1:num_stims
            % Extract stimulus properties from trial row
            st = trial_dat(s,:);
            
            % Check if stimulus is activated by user action (finger lift,
            % saccade) or scheduled event, then store timing info in sched
            if st.stim_onset < 0
                % If stim onset time is -1, it is activated by finger lift
                % if stim onset time is -2, it is activated by saccade
                
                % If stim is activated & deactivated by same action => null
                if st.stim_onset == st.stim_duration
                    stim_schedule(s, :) = 0;
                % If stim is finger lift activated and saccade deactivated
                elseif st.stim_onset == -1 && st.stim_duration == -2
                    stim_schedule(s, :) = -12;
                % If stim is saccade activated and finger lift deactivated
                elseif st.stim_onset == -2 && st.stim_duration == -1
                    stim_schedule(s, :) = -21;
                % If -3, stim stays on until end of trial after activated
                elseif st.stim_duration == -3
                    stim_schedule(s, :) = st.stim_onset;
                % If stim has duration after activated, stim_schedule will 
                % be updated once stim is activated during the trial
                elseif st.stim_duration > 0
                    % Duration Flag: Finger lift = -101, Saccade = -102
                    stim_schedule(s, :) = -100 + st.stim_onset;
                end
                
            elseif st.stim_onset >= 0
                % Find first frame stim is displayed
                stim_img_on = ceil(st.stim_onset / ifi) + 1;
                
                % If stim is finger lift deactivated
                if st.stim_duration == -1
                    stim_schedule(s, stim_img_on:end) = -11;
                % If stim is saccade deactivated
                elseif st.stim_duration == -2
                    stim_schedule(s, stim_img_on:end) = -22;
                % If -3, stim stays on until end of trial after activated
                elseif st.stim_duration == -3
                    stim_schedule(s, stim_img_on:end) = 1;
                % If stim has duration after activated, fill all frames
                % where this stim is displayed with a 1
                elseif st.stim_duration > 0
                    stim_img_off = ceil((st.stim_onset+st.stim_duration) / ifi) + 1;
                    if stim_img_on ~= stim_img_off
                        stim_schedule(s, stim_img_on:stim_img_off) = 1;
                    end
                end

                % Fill all frames where this stim's mask is displayed with a 1
                mask_img_on = ceil(st.mask_onset / ifi) + 1;
                mask_img_off = ceil((st.mask_onset+st.mask_duration) / ifi) + 1;
                if mask_img_on ~= mask_img_off
                    mask_schedule(s, mask_img_on:mask_img_off) = 1;
                end
            end
            
            % If stim is a fixation event, save in fixation schedule
            if st.subj_fixation_type ~= 0
                stim_fix_on = max(ceil(st.subj_fixation_onset / ifi),1);
                fixation_schedule(s, stim_fix_on) = st.subj_fixation_type;
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
                                    (TOUCH_MARGIN_CM*ptb.ppcm);
            
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
            mask_dotsizes = ones(1, mask_ndots) .* max(1, st.mask_size*ptb.ppcm);


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

        % Nominal trial frames
        nominal_frames = nan(1, trial_frames);
        
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
        isFixationPause = false;
        
        % Keyboard response
        kb_resp = 'None';
        
        % Frame counters: (frame) for stimuli schedule
        % and (glob_frame) for total elapsed (incl. fixation pause)
        frame = 1;
        glob_frame = 1;
        init_frame = 1; % actual first frame (if fixation pause on frame 1)
        
        % This state is true prior to target onset, false after (used for
        % reaction time calculation)
        pre_target = 1;
        
        % If we have passed fixation requirement, save this state as we
        % will begin reactime time on saccade/finger lift
%         post_fix = 0;
        post_fix_gaze = 0;
        post_fix_touch = 0;
        
        % Flags if fingerlift and/or saccade has occured, used to display 
        % stims which depend on these actions (stim_onset = -1 or -2)
        is_fingerlift = false;
        is_saccade = false;
        
        % Detect starting touch/gaze and require lift/saccade to begin
        % recording data
        touch_on_start = 0;
        gaze_on_start = 0;
        
        % Set gaze counter for gaze targets
        gaze_counter = [NaN, 0];
        
        % Detect if finger is on touch fixation
        onFix = 0;
        on_fix_rt_delay = 0;
        
        % Initialize Optotrak recording
        if useOptotrak
            puRealtimeData=0;puSpoolComplete=0;
            puSpoolStatus=0;pulFramesBuffered=0;
            
            %determine the trial number text
            txtTrial = [num2str(tn, '%03d') '_']; % pad left side with zeros
            
            %determine block number text
            block_txt = [num2str(bn, '%02d') '_']; % pad left side with zeros

%             optotrack_trial_mapping = [optotrack_trial_mapping; [int2str(tn) current_trial_name]];
            %navigate to the NDI Raw data file
            NDIFPath = [out_path part_dems.id '/' fname '/' optk.FNRAWDATA '/'];

            %The NDI dat file name for the current trial
            NDIFName = [part_dems.id optk.OPTO txtTrial block_txt '.dat'];

            if t==1
                %re-initialize the reorganized data
                dataReorg = zeros(optk.fTotPerTrial, optk.nMarkers*4);
            end
    % -------------- Start the optotrack recording here! -------------------- %
            %start collecting data, the number of frames to collect was
            %pre-specified by the function OPTOTRAKSETUPCOLLECTION.

            %initialize the file for spooling
            OptotrakActivateMarkers();
            WaitSecs(0.010);
            DataBufferInitializeFile(0,[NDIFPath NDIFName]);
            DataBufferStart();
        end

        % Begin recording EyeLink data
        if useEyelink
            Eyelink('StartRecording');
        end
        
        % -- Display trial stimuli, record responses ----------------------

        % Flip again to sync us to the vertical retrace
        vbl = Screen('Flip', window);
        
        % Display all stimuli which are scheduled for the current frame
        while frame < trial_frames + 1

            % Listen for Esc key to abort experiment
            [keyIsDown, ~, keyCode] = KbCheck;
            
            % If F3 pressed during trial, abort current trial (recycle)
            % and calibrate EyeLink before starting next trial
            if useEyelink && keyIsDown && keyCode(ptb.F3Key)
                pre_trial_calibration = true;
                bad_trials(t) = 3;
                break
            end
            
            if keyIsDown && keyCode(ptb.escapeKey)
                abort_experiment = true;
                break
            end
            
            % Note when a background color change has occurred in the
            % current frame
            bg_color_change = false;

            % Get the current position of the mouse
            [tx, ty, buttons] = GetMouse(window);
            is_touch = any(buttons ~= 0);

            % Detect touch/gaze on trial start
            if init_frame == 1
                if is_touch
                    touch_on_start = true;
                end
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
            
            % If trial started with gaze fixation, require a saccade
            % to register a new touch
            if gaze_on_start
                %ADDGAZE
            end
            
            % We clamp the values at the maximum values of the screen
            % in X and Y in case people have two monitors connected.
            tx = min(tx, ptb.screenXpixels);
            ty = min(ty, ptb.screenYpixels);
% 
%             % Store touch location and whether finger is touching
%             trial_touch_samples(1, glob_frame) = tx;
%             trial_touch_samples(2, glob_frame) = ty;
%             trial_touch_samples(3, glob_frame) = is_touch;
            
            % Get current gaze value
            if useEyelink
                % Check recording status, restart EyeLink if error
                el_error=Eyelink('CheckRecording');
                if(el_error~=0)
                    eyelink_bad_count = eyelink_bad_count + 1;
                    % don't let one bad data point kill ya
                    if(eyelink_bad_count > el.lost_gaze_max_frames)
%                         disp('Eyelink Error!');
%                         % attempt to restart Eyelink recording
%                         Eyelink('StopRecording');
%                         WaitSecs(0.005);
%                         Eyelink('StartRecording');
                    end
                end
                
                % check for presence of a new sample update
                if Eyelink('NewFloatSampleAvailable') > 0
                    
                    % get the sample in the form of an event structure
                    evt = Eyelink( 'NewestFloatSample');
                    
                    % get current gaze position from sample
                    gx = evt.gx(eye_used+1); % +1 indexing MATLAB array
                    gy = evt.gy(eye_used+1);
                    
                    % do we have valid data and is the pupil visible?
                    if gx~=el.MISSING_DATA && gy~=el.MISSING_DATA %|| ~(evt.pa(eye_used+1)>0)
%                         gaze_good = true; 
                        eyelink_bad_count = 0; % reset the counter for bad 
                    else
                        gx = NaN;
                        gy = NaN;
                        eyelink_bad_count = eyelink_bad_count + 1;
                    end
                end
                      
                if eyelink_bad_count > el.lost_gaze_max_frames
                    bad_trials(t) = -1;
                    break
                end
                
                % Test if saccade has occurred, start response period 2
                if post_fix_gaze == 1 && isnan(rt_start_saccade)
                    % If gaze no longer persists on pfidx (previous fixaition)
                    if ~ptb_in_circ(stim_centers(gpfidx,:), [gx gy], stim_radius(gpfidx)*1.5)
                        start_saccade = true;
                        post_fix_gaze = 0;
                    end
                end   
            else
                gx = NaN;
                gy = NaN;
            end
%             
%             % Add gaze sample to trial array
%             trial_gaze_samples(1, glob_frame) = gx;
%             trial_gaze_samples(2, glob_frame) = gy;
               
            % Test if finger lift has occurred, start response period 2
            if post_fix_touch == 1 && ~is_touch && isnan(rt_start_fingerlift)
                start_fingerlift = true;
                post_fix_touch = 0;
            end
            
            % Check if fixation appears in current frame
            active_fixations = find(fixation_schedule(:,frame) ~= 0);
            
            % If in fixation pause for current frame
            if isFixationPause
                % Get pause fixation indices
                pause_fix = (trial_dat.subj_fixation_pause(active_fixations) == 1);
                
                % Find any satisfied pause fixations and deactivate
                pause_fix_done = fixation_counter(pause_fix,3) == 0;
                active_fixations(pause_fix_done) = [];
                
                % If all pause fixations are satisfied, deactivate pause
                if all(pause_fix_done)
                    isFixationPause = false;
                end
            end
            
            % Populate fixation counter with any currently active fixations
            if ~isempty(active_fixations)
                if any(trial_dat.subj_fixation_pause(active_fixations))
                    isFixationPause = true;
                end
%                 post_fix = 0;
                for f=1:length(active_fixations)
                    fidx = active_fixations(f);
                    % Set fixation to active
                    fixation_counter(fidx,3) = true;
                    % Note fixation type in counter
                    fixation_counter(fidx,4) = fixation_schedule(fidx,frame);
                    if fixation_counter(fidx,4) == 1
                        post_fix_touch = 0;
                    elseif fixation_counter(fidx,4) == 2
                        post_fix_gaze = 0;
                    end
                end
            end
            
            % Test each active fixation for gaze/touch and iterate counters       
            active_fixations = find(fixation_counter(:,3) ~= 0);
            if isempty(active_fixations)
                if isFixationPause
                    isFixationPause = false;
                end
                post_fix_touch = 1;
                post_fix_gaze = 1;
            else
                if (~any(trial_dat.subj_fixation_type(active_fixations) == 1))
                    post_fix_touch = 1;
                end
                
                if (~any(trial_dat.subj_fixation_type(active_fixations) == 2))
                    post_fix_gaze = 1;
                end
                
                if isFixationPause
                    active_p_fixations = fixation_schedule(:,frame) ~= 0;
                    pause_fix = find(trial_dat.subj_fixation_pause(active_p_fixations) == 1);
                    for pf=1:length(pause_fix)
                        pfidx = pause_fix(pf);
                        % Render fixation image
                        pfix_name = char(trial_dat.stim_img_name(pfidx));
                        Screen('DrawTextures', window,...
                                stim_textures(pfix_name),[],...
                                stim_bounds(pfidx,:), ...
                                trial_dat.stim_rotation(pfidx));
                        stim_presentations(pfidx, glob_frame) = 1;
                    end
                end
                
                for f=1:length(active_fixations)
                    fidx = active_fixations(f);
                    
                    % Detect if background color changed
                    if ~bg_color_change
                        [ptb, bg_color_change] = ptb_set_bg_color(cell2mat(trial_dat.background_color(fidx)), cell2mat(trial_dat.text_color(fidx)), ptb);
                    end
                    
                    % If fixation is touch type check for touch, reset if
                    % touch lifted
                    if fixation_counter(fidx,4) == 1
                        if is_touch && ptb_in_circ(stim_centers(fidx,:), [tx ty], stim_radius(fidx))
                            fixation_counter(fidx,2) = fixation_counter(fidx,2) + 1;
                            fixation_counter(fidx,5) = 0;
                        else
                            fixation_counter(fidx,2) = 0;
                            if ~isFixationPause
                                fixation_counter(fidx,5) = fixation_counter(fidx,5) + 1;
                            end
                        end % touched fixation mask
                    end
                        
                    % If fixation is gaze type, check for gaze, reset if
                    % gaze lost
                    if fixation_counter(fidx,4) == 2
                        if ~isnan(gx) && ~isnan(gy) 
                            if ptb_in_circ(stim_centers(fidx,:), [gx gy], stim_radius(fidx)*1.5)
                                fixation_counter(fidx,2) = fixation_counter(fidx,2) + 1;
                                fixation_counter(fidx,5) = 0;
                            else
                                fixation_counter(fidx,2) = 0;
                                if ~isFixationPause
                                    fixation_counter(fidx,5) = fixation_counter(fidx,5) + 1;
                                end
                            end % touched fixation mask
                        end
                    end
                    
                    if fixation_counter(fidx,5) > lost_fix_max_frames
                        bad_trials(t) = 1;
                        break
                    end
                    
                    % If participant fixation duration exceeds requirement,
                    % deactivate fixation
                    if fixation_counter(fidx,1) < fixation_counter(fidx,2)
                        fixation_counter(fidx,3) = false;
                        % Set fixation ID for gaze persistences
                        if fixation_counter(fidx,4) == 2
                            gpfidx = fidx;
                        end
                    end
                end
            end
            
            if bad_trials(t) == 1
                break
            end
            
            % Check image and mask schedules and render anything to screen
            % that is scheduled to be shown on the current frame
            if ~isFixationPause
                for s=1:num_stims
                    % Get current entry in stim_schedule
                    stim_sch_val = stim_schedule(s, frame);
                    
                    % If schedule entry conditions are filled, then display
                    % stim
                    display_stim = false;
                    if stim_sch_val == 1
                        display_stim = true;
                    elseif stim_sch_val == -1 && is_fingerlift
                        display_stim = true;
                    elseif stim_sch_val == -2 && is_saccade
                        display_stim = true;
                    elseif stim_sch_val == -11 && ~is_fingerlift
                        display_stim = true;
                    elseif stim_sch_val == -22 && ~is_saccade
                        display_stim = true;
                    elseif stim_sch_val == -12 && is_fingerlift && ~is_saccade
                        display_stim = true;
                    elseif stim_sch_val == -21 && ~is_fingerlift && is_saccade
                        display_stim = true;  
                    elseif stim_sch_val == -101 && is_fingerlift
                        display_stim = true;  
                        % Duration flag invoked - complete stim schedule
                        stim_img_off = ceil((frame+trial_dat.stim_duration(s)) / ifi) + 1;
                        if frame < trial_frames
                        %if stim_sched_frame < trial_frames    
                            if stim_img_off > trial_frames 
                                stim_schedule(s, frame+1:end) = 1;
                            else
                                stim_schedule(s, frame+1:stim_img_off-1) = 1;
                                stim_schedule(s, stim_img_off:end) = 0;
                            end
                        end
                    elseif stim_sch_val == -102 && is_saccade
                        display_stim = true;  
                        % Duration flag invoked - complete stim schedule
                        stim_img_off = ceil((frame+trial_dat.stim_duration(s)) / ifi) + 1;
                        if frame < trial_frames
                        %if stim_sched_frame < trial_frames    
                            if stim_img_off > trial_frames 
                                stim_schedule(s, frame+1:end) = 1;
                            else
                                stim_schedule(s, frame+1:stim_img_off-1) = 1;
                                stim_schedule(s, stim_img_off:end) = 0;
                            end
                        end
                    end
                    
                    if display_stim
                        stim_name = char(trial_dat.stim_img_name(s));

                        % Detect if backgroudnd color changed
                        if ~bg_color_change
                            [ptb, bg_color_change] = ptb_set_bg_color(cell2mat(trial_dat.background_color(s)), cell2mat(trial_dat.text_color(s)), ptb);
                        end

                        % Display stim
                        Screen('DrawTextures', window,...
                                stim_textures(stim_name),[],...
                                stim_bounds(s,:), trial_dat.stim_rotation(s));

                        % Play associated sound if applicable
                        if useSound && ~snd_played(s) && ~strcmp(trial_dat.stim_sound_file(s), 'None')
                            snd_name = char(trial_dat.stim_sound_file(s));
                            Snd('Play', snd_buffers(snd_name), fs);
%                           PsychPortAudio('FillBuffer', pa_handle, snd_buffers(snd_name)); % loads data into buffer
%                             PsychPortAudio('Start', pa_handle, 1, 0); %starts sound immediatley
                            snd_played(s) = 1;
                        end
                        
                        % Record stim presentation in global frame counter 
                        stim_presentations(s, glob_frame) = 1;
                        stim_displayed(s) = 1;
                        
                        % If stim is a target, begin reaction time counter
                        if any(targetStims == s) && pre_target
                            start_target_onset = true;
                        end
                    end % image scheduled
                    if mask_schedule(s, frame) == 1
                    %if mask_schedule(s, stim_sched_frame) == 1    
                        % Draw mask dots
                        Screen('DrawDots', window, dotPosMatrix{s},...
                                dotSizes{s}, dotColors{s}, [], 2);

                        % Record mask presentation in global frame counter
                        % Stim only = 1, mask only = 2, stim & mask = 3
                        stim_presentations(s, glob_frame) = stim_presentations(s, glob_frame) + 2;
                    end % mask scheduled
                end % stim loop
                
                % Check for collisions
                stim_is_touched = false;
                gaze_on_object = false;
                for s=1:num_stims
                    if stim_displayed(s) && ~stim_is_touched
                        % If stim is touchable, detect if touched
                        if trial_dat.stim_is_touchable(s) == 1
                            if ptb_in_circ(stim_centers(s,:), [tx ty], stim_radius(s)) && is_touch
                                stim_is_touched = true;
                            end% touched stim
                        end
                        
                        % If stim is gaze target, detect gaze
                        if trial_dat.stim_is_touchable(s) == 2
                            if ptb_in_circ(stim_centers(s,:), [gx gy], stim_radius(s)) && ~gaze_on_start
                                gaze_on_object = true;
                                if gaze_counter(1) == s
                                    gaze_counter(2) = gaze_counter(2) + 1;
                                else
                                    gaze_counter(1) = s;
                                    gaze_counter(2) = 0;
                                end
                                if gaze_counter(2) > gaze_choice_min_frames
                                    stim_is_touched = true;
                                end
                            end% touched stim
                        end
                        
                        % If stim is touch or gaze target, detect either
                        if trial_dat.stim_is_touchable(s) == 3
                            if ptb_in_circ(stim_centers(s,:), [tx ty], stim_radius(s)) && is_touch
                                stim_is_touched = true;
                            end% touched stim
                            if ptb_in_circ(stim_centers(s,:), [gx gy], stim_radius(s)) && ~gaze_on_start
                                gaze_on_object = true;
                                if gaze_counter(1) == s
                                    gaze_counter(2) = gaze_counter(2) + 1;
                                else
                                    gaze_counter(1) = s;
                                    gaze_counter(2) = 0;
                                end
                                if gaze_counter(2) > gaze_choice_min_frames
                                    stim_is_touched = true;
                                end
                            end% touched stim
                        end
                        
                        % If stim is touch and gaze target, detect both
                        if trial_dat.stim_is_touchable(s) == 4
                            if ~gaze_on_start && ptb_in_circ(stim_centers(s,:), [gx gy], stim_radius(s))
                                stim_is_touched = true;
                                gaze_on_object = true;
                                if gaze_counter(1) == s
                                    gaze_counter(2) = gaze_counter(2) + 1;
                                else
                                    gaze_counter(1) = s;
                                    gaze_counter(2) = 0;
                                end
                                if (gaze_counter(2) > gaze_choice_min_frames && ...
                                         is_touch && ptb_in_circ(stim_centers(s,:), [tx ty], stim_radius(s)))
                                    stim_is_touched = true;
                                end
                            end% touched stim
                        end
                        
                        % If stim is touch fixation, detect if touched
                        if trial_dat.subj_fixation_type(s) == 1
                            if ptb_in_circ(stim_centers(s,:), [tx ty], stim_radius(s)) && is_touch
                                stim_is_touched = true;
                            end% touched stim
                        end
                        
                        if stim_is_touched
                            if trial_dat.subj_fixation_type(s) == 1
                                if ~(fixation_counter(s,2)<fixation_counter(s,1)) && ~trial_dat.subj_fixation_pause(s)
                                    on_fix_rt_delay = on_fix_rt_delay + 1;
                                    if on_fix_rt_delay > on_fix_rt_max
                                        onFix = false;
                                    else
                                        onFix = true;
                                    end
                                else
                                    onFix = true;
                                    on_fix_rt_delay = 0;
                                end
                            else
                                hasSelected = true;
                                onFix = false;
                                selectedStim = s;
                                if bad_trials(t)
                                    bad_trials(t) = false;
                                end
                            end
                        elseif is_fingerlift && is_touch
                            onFix = false;
                        else
                            onFix = true;
                        end
                    end % image displayed & touchable
                end
                
                % Reset gaze counter if no gaze found on any object
                if ~gaze_on_object
                    gaze_counter(1) = NaN;
                    gaze_counter(2) = 0;
                end
                
                % Check for any non-stim touches, add trial to list of bad
                % trials to repeat at end of block
                if any(stim_displayed) && ~hasSelected && is_touch && ~onFix
                    hasSelected = true;
                    selectedStim = NaN;
                    bad_trials(t) = 2;
                end
            end % if fixation pause
            
            % Flip to the screen
            [vbl, sot, flip, miss, ~] = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
                        
            % Store exact timestamp and max error for current frame
            trial_timestamps(:, glob_frame) = [vbl, sot, flip, miss];
            
            % If missed frame, iterate 2 frames to stay on track.
            missed_frame = (miss > 0);
            next_frame = 1 + missed_frame;
            
            % If missed frame, copy fixation events to subsequent frame so
            % as not to skip them (only if they are blank)
            if missed_frame && (frame + 1 < trial_frames)
                if fixation_schedule(:,frame+1) ~= 0
                    fixation_schedule(:,frame+2) = fixation_schedule(:,frame+1);
                end
            end
          
            % If a touchable mask has been touched, store response time
            if hasSelected == true
                rt_target_to_choice = frame - rt_start_target_onset;
                rt_saccade_to_choice = frame - rt_start_saccade;
                rt_fingerlift_to_choice = frame - rt_start_fingerlift;
            end

            % Note start time to calculate reaction time 1
            if start_target_onset
%                 KbStrokeWait();
                rt_start_target_onset = frame;
                start_target_onset = false;
                pre_target = false;
            end
            
            % If a fixation has occured, record finger lift or saccade
            % reaction time
            % If saccade occurs, note time & start rt counter
            if start_saccade
                is_saccade = true;
                rt_start_saccade = frame;
                rt_target_to_saccade = frame - rt_start_target_onset;
                start_saccade = false;
            end
            
            % If a fixation has occured, record finger lift or saccade
            % reaction time
            % If finger lift occurs, note time & start rt counter
            if start_fingerlift
                is_fingerlift = true;
                rt_start_fingerlift = frame;
                rt_target_to_fingerlift = frame - rt_start_target_onset;
                start_fingerlift = false;
            end
            
            [~, ~, trial_keyCode] = KbCheck();
            
            % End trial if selection has been made
            if hasSelected == true || (trial_dat.trial_feedback_type(end)>1 && any(trial_keyCode))
                nominal_frames(glob_frame) = frame;
                if ~isnan(selectedStim) 
                    correct_choice = any(targetStims == selectedStim);
                else
                    correct_choice = NaN;
                end
                break
            end
            
            % Store touch location and whether finger is touching
            trial_touch_samples(1, glob_frame) = tx;
            trial_touch_samples(2, glob_frame) = ty;
            trial_touch_samples(3, glob_frame) = is_touch;
            
            % Add gaze sample to trial array
            trial_gaze_samples(1, glob_frame) = gx;
            trial_gaze_samples(2, glob_frame) = gy;
            
            % Turn off initial frame
            if init_frame
                init_frame = false;
            end
            
            % If not in a fixation pause, iterate stim schedule
            if ~isFixationPause
                % Iterate trial frame counter
                frame = frame + next_frame;
            end
            
            % Save nominal frame for data output 
            nominal_frames(glob_frame) = frame;
            glob_frame = glob_frame + 1;
            
        end % frame loop
        
%         if useSound
%             Snd('Stop');
% %             PsychPortAudio('Stop', pa_handle);% Stop sound playback
%         end
        
        % Stop recording EyeLink data
        if useEyelink
            Eyelink('StopRecording');
        end
        
        if useOptotrak
            DataBufferStop();
            puSpoolComplete = 0;
            fprintf('spooling data\n');
            %transfer data from the optotrak to the computer
            while (puSpoolComplete == 0)
                % call C library function here. See PDF
                [puRealtimeData,puSpoolComplete,puSpoolStatus,pulFramesBuffered]=DataBufferWriteData(puRealtimeData,puSpoolComplete,puSpoolStatus,pulFramesBuffered);
                WaitSecs(.1);
                if puSpoolStatus ~= 0
                    disp('Spooling Error');
                    break;
                end
            end
            disp('spooling complete. deactivating markers');
            
            %Deactivate the Markers
            OptotrakDeActivateMarkers();
            
            disp(['attempting to read file: ' deblank([NDIFPath NDIFName])]);
            %open the .dat file stored by Optotrak
            [fid, ~] = fopen(deblank([NDIFPath NDIFName]),'r+','l'); % little endian byte ordering

            %Read the header portion of the file
            [~,items,subitems,numframes] = read_nd_c_file_header(fid);

            %Read the data portion of the file
            rawData = read_nd_c_file_data(fid, items, subitems, numframes);
            fclose(fid);
            
            disp('converting data');
            %Convert the data to a format OTCTextWrite accepts
            for i = 1:optk.nMarkers
                %append the marker data to dataReorg and delete the 
                dataReorg(1:size(rawData,3),((i+1)*3)-2:(i+1)*3) = transpose(squeeze(rawData(i,:,:)));
            end
    
            %Convert missing data to NaNs
            dataReorg(dataReorg < optk.MISSINGDATACUTOFFVALUE) = NaN;
        
            %add the data collection information to the data
            dataReorg(:,1) = tn;
            dataReorg(:,2) = optk.smpRt;
        
            %navigate to the 'OTCReorganized' file folder
            NDIDatPath = [out_path part_dems.id '/' fname '/' optk.FNOTCREORG '/'];
            NDIDatFile = [part_dems.id optk.OPTO txtTrial block_txt '.txt'];
            %ALWAYS STORE A RAW DATA SET
            disp(['Saving file: ' NDIDatFile ' ...'])
            fid = fopen([NDIDatPath NDIDatFile], 'w');
            %fprintf(fid, hdr);
            fclose(fid);
            dlmwrite([NDIDatPath NDIDatFile],...
                        dataReorg,'-append','delimiter',...
                        '\t','newline','pc','precision',12);
            disp('Success!')
        end
        
        % Save total trial time
        trial_frames_elapsed = min(frame, trial_frames);
        
        % If trial timed out, display time out message
        if trial_frames_elapsed >= trial_frames
            if trial_dat.trial_timeout_msg(num_stims)
                % Set text and background color to match first trial
                if ~bg_color_change
                    [ptb, bg_color_change] = ptb_set_bg_color(cell2mat(trial_dat.background_color(num_stims)), cell2mat(trial_dat.text_color(num_stims)), ptb);
                end

                % Draw timeout text
                Screen('Flip', window);
                DrawFormattedText(window, timeoutMsg, ...
                                    'center', 'center', ptb.text_color);
                Screen('Flip', window);
                [~, keyCode, ~] = KbStrokeWait;
                if keyCode(ptb.escapeKey)
                    abort_experiment = true;
                end
            end
        elseif ~abort_experiment
            % Display bad trial message
            if bad_trials(t) == 1 || bad_trials(t) == 2
                % Set text and background color to match first trial
                if ~bg_color_change
                    [ptb, bg_color_change] = ptb_set_bg_color(cell2mat(trial_dat.background_color(num_stims)), cell2mat(trial_dat.text_color(num_stims)), ptb);
                end

                % Clear screen
                Screen('Flip', window);
                
                if strcmp(trial_dat.trial_feedback_img(end), 'None')
                    feedbackMsg = feedback_message_incorrect;
                    Screen('Flip', window);
                    DrawFormattedText(window, feedbackMsg, 'center', 'center', ptb.text_color);
                else
                    fb_fname = char(trial_dat.trial_feedback_img(num_stims));
                    Screen('DrawTextures', window, fb_textures(fb_fname), [], [], 0);
                end
                
                Screen('Flip', window);
                pause(1);
            elseif bad_trials(t) == 3
                Screen('Flip', window);
            elseif bad_trials(t) == -1
                % Set text and background color to match first trial
                if ~bg_color_change
                    [ptb, bg_color_change] = ptb_set_bg_color(cell2mat(trial_dat.background_color(num_stims)), cell2mat(trial_dat.text_color(num_stims)), ptb);
                end

                % Clear screen
                Screen('Flip', window);
                
                gazeLostMsg = gaze_message_lost_gaze;
                Screen('Flip', window);
                DrawFormattedText(window, gazeLostMsg, 'center', 'center', ptb.text_color);
                
                Screen('Flip', window);
                [~, keyCode, ~] = KbStrokeWait;
                if keyCode(ptb.escapeKey)
                    abort_experiment = true;
                    break
                end
                if useEyelink && keyCode(ptb.F3Key)
                    ptb_calibrate_eyelink( calibration_sc, el, ptb );
                end
            % If trial keyboard response is on, prompt for keypress
            else
                if mod(trial_dat.trial_kb_resp(end),2)==1 % trial_dat.trial_kb_resp(end)
                    % Set text and background color to match first trial
                    if ~bg_color_change
                        [ptb, bg_color_change] = ptb_set_bg_color(cell2mat(trial_dat.background_color(num_stims)), cell2mat(trial_dat.text_color(num_stims)), ptb);
                    end

                    % Draw timeout text
                    Screen('Flip', window);
                    DrawFormattedText(window, kbrespMsg, ...
                                        'center', 'center', ptb.text_color);
                    Screen('Flip', window);
                    [~, keyCode, ~] = KbStrokeWait;
                    if keyCode(ptb.escapeKey)
                        abort_experiment = true;
                        break
                    end
                    kb_resp = KbName(keyCode);
                elseif trial_dat.trial_kb_resp(end) == 2
                    kb_resp = KbName(trial_keyCode);
                end
            end
                
            % Display feedback
            if trial_dat.trial_feedback_type(end) && bad_trials(t) == 0
                % Set text and background color to match first trial
                if ~bg_color_change
                    [ptb, bg_color_change] = ptb_set_bg_color(cell2mat(trial_dat.background_color(num_stims)), cell2mat(trial_dat.text_color(num_stims)), ptb);
                end

                % Clear screen
                Screen('Flip', window);
                
                % Render feedback image or text
                if ~(strcmp(kb_resp, trial_dat.correct_kb_resp(num_stims)))
                    if strcmp(trial_dat.trial_feedback_img(end), 'None')
                        feedbackMsg = feedback_message_incorrect;
                        Screen('Flip', window);
                        DrawFormattedText(window, feedbackMsg, 'center', 'center', ptb.text_color);
                    else
                        fb_fname = char(trial_dat.trial_feedback_img(num_stims));
                        Screen('DrawTextures', window, fb_textures(fb_fname), [], [], 0);
                    end
                else
                    if strcmp(trial_dat.trial_feedback_img(end), 'None')
                        feedbackMsg = feedback_message_correct;
                        Screen('Flip', window);
                        DrawFormattedText(window, feedbackMsg, 'center', 'center', ptb.text_color);
                    end
                end
                
                Screen('Flip', window);
                pause(0.6);
            end
        end
        
        % Output frame-by-frame timestamp/presentation record
        data_num_frames_min = min([length(nominal_frames), ...
                                   length(trial_timestamps), ...
                                   length(stim_presentations), ...
                                   length(trial_touch_samples), ... 
                                   length(trial_gaze_samples)]);
        rfr = 1:data_num_frames_min;
        
        timestamp_record = [nominal_frames(rfr); trial_timestamps(:,rfr); stim_presentations(:,rfr); trial_touch_samples(:,rfr); trial_gaze_samples(:,rfr)]';
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

        % Allocate new otuput table row if necessary
        if row_ct > size(partData,1)
            partData = [partData; blank_data_row];
        end
        
        % Record the trial data into out data matrix
        partData.participant{row_ct} = part_dems.id;
        partData.BlockNum{row_ct} = bn;
        partData.TrialNum{row_ct} = tn;
        partData.BadTrial{row_ct} = bad_trials(t); % TBD
        
        % Extract target stimulus info from trial stim table if possible
        if ~isempty(targetStims)
            tr = trial_dat(targetStims,:);
            partData.target_image_name{row_ct} = strjoin((tr.stim_img_name),' ');
            partData.target_center_x{row_ct} = strjoin(cellstr(num2str(tr.stim_cent_x*ptb.screenXpixels)),' ');
            partData.target_center_y{row_ct} = strjoin(cellstr(num2str(tr.stim_cent_y*ptb.screenYpixels)),' ');
            partData.target_size_x{row_ct} = strjoin(cellstr(num2str(tr.stim_size_x*ptb.ppcm)),' ');
            partData.target_size_y{row_ct} = strjoin(cellstr(num2str(tr.stim_size_y*ptb.ppcm)),' ');
            partData.target_rotation{row_ct} = strjoin(cellstr(num2str(tr.stim_rotation)),' ');
            partData.target_image_num{row_ct} = strjoin(cellstr(num2str(targetStims)),' ');
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
            
            if ~strcmp(ch.correct_kb_resp, 'None')
                kb_correct = any(strcmp(ch.correct_kb_resp, kb_resp));
            else
                kb_correct = NaN;
            end
            partData.keyboard_corr_response{row_ct} = cell2mat(ch.correct_kb_resp);
            partData.keyboard_correct{row_ct} = kb_correct;
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
        
        if iscell(kb_resp)
            partData.keyboard_user_response{row_ct} = [kb_resp{:}];
        else
            partData.keyboard_user_response{row_ct} = kb_resp;
        end
        
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
            if ~bg_color_change
            	[ptb, bg_color_change] = ptb_set_bg_color(cell2mat(trial_dat.background_color(num_stims)), cell2mat(trial_dat.text_color(num_stims)), ptb);
            end

            % Draw practice block or generic block end message
            if trial_dat.is_practice_block(1)
            	DrawFormattedText(window, practiceEndMsg, ...
                                'center', 'center', ptb.text_color);
            else
                DrawFormattedText(window, blockEndMsg, ...
                                'center', 'center', ptb.text_color);
            end
            Screen('Flip', window);
            [~, keyCode, ~] = KbStrokeWait;
            if keyCode(ptb.escapeKey)
                abort_experiment = true;
                break
            end
            
            if useEyelink && keyCode(ptb.F3Key)
                ptb_calibrate_eyelink( calibration_sc, el, ptb );
            end
            
        end % bn < num_blocks
    else
        break
    end % abort_experiment
    
end % block loop

% Flip again to sync us to the vertical retrace
Screen('Flip', window);

% End of experiment screen. We clear the screen once they have made their
% response
% Set text and background color to match first trial
if ~bg_color_change
    ptb = ptb_set_bg_color(cell2mat(trial_dat.background_color(num_stims)), cell2mat(trial_dat.text_color(num_stims)), ptb);
end

% Draw End of Experiment text
DrawFormattedText(window, expEndMsg, 'center', 'center', ptb.text_color);
Screen('Flip', window);
KbStrokeWait;
ShowCursor('arrow', window);

if useSound
    Snd('Close');
%     PsychPortAudio('Close', pa_handle);% Close the audio device
end

% End EyeLink recording and save data file
if useEyelink
    edf_file_out = [part_out_path fname '_' el.edfFile];
    ptb_end_eyelink(el.edfFile, edf_file_out);
end
    
sca;