function end_code = ptb_calibrate_eyelink( calibration_sc, el, ptb )
%PTB_CALIBRATE_EYELINK Perform EyeLink calibration and verify on screen

    eyelink_bad_count = 0;

    eye_used = el.RIGHT_EYE;
    
    EyelinkDoTrackerSetup(el); % do EyeLink calibration, validation, etc.
    
    Eyelink('StartRecording');
    
    % Clear screen
    Screen('Flip', ptb.window);
    
    isLoop = 1;
    
    keyIsDown = KbCheck;
    if keyIsDown
        waitforKeyRelease = true;
    else 
        waitforKeyRelease = false;
    end   
     
    while isLoop
        
        % Check recording status, restart EyeLink if error
        el_error=Eyelink('CheckRecording');
        
        if(el_error~=0)
            eyelink_bad_count = eyelink_bad_count + 1;
            % don't let one bad data point kill ya
            if(eyelink_bad_count > el.lost_gaze_max_frames)
                disp('Eyelink Error!');
                % attempt to restart Eyelink recording
                Eyelink('StopRecording');
                WaitSecs(0.005);
                Eyelink('StartRecording');
            end
        end

        % check for presence of a new sample update
        if Eyelink('NewFloatSampleAvailable') > 0

            % get the sample in the form of an event structure
            evt = Eyelink('NewestFloatSample');

            % get current gaze position from sample
            gx = evt.gx(eye_used+1); % +1 indexing MATLAB array
            gy = evt.gy(eye_used+1);

            % do we have valid data and is the pupil visible?
            if gx~=el.MISSING_DATA && gy~=el.MISSING_DATA 
                eyelink_bad_count = 0; % reset the counter for bad 
            else
                gx = NaN;
                gy = NaN;
                eyelink_bad_count = eyelink_bad_count + 1;
            end
        else
            gx = NaN;
            gy = NaN;
        end

        if eyelink_bad_count > el.lost_gaze_max_frames
            disp('Eyelink Error!');
            % attempt to restart Eyelink recording
            Eyelink('StopRecording');
            WaitSecs(0.005);
            Eyelink('StartRecording');
        end
    
        Screen('DrawTextures', ptb.window,...
                calibration_sc('cali_all'),[],...
                [0 0 ptb.screenXpixels ptb.screenYpixels], 0);
            
    	gazeRect=[gx-7 gy-7 gx+8 gy+8];         
        gx_str = num2str(gx);
        gy_str = num2str(gy);     
        mos_pos = strcat(gx_str, ',', gy_str);         
        Screen('FrameOval', ptb.window, ptb.black, gazeRect, 6, 6);       
        Screen('DrawText', ptb.window, mos_pos, gx+50, gy+50);
            
        Screen('Flip', ptb.window);
        
        [keyIsDown, ~, keyCode, ~] = KbCheck;
        
        if ~keyIsDown && waitforKeyRelease
            waitforKeyRelease = false;
        end
        
        if keyIsDown && ~waitforKeyRelease
            if keyCode(ptb.escapeKey)
                end_code = 1;
            else
                end_code = 0;
            end
            isLoop = false;
        end
        
    end % while loop
    
    Eyelink('StopRecording');
    
    Screen('FillRect', ptb.window, ptb.bg_color);
    
    % Clear screen
    Screen('Flip', ptb.window);
end

