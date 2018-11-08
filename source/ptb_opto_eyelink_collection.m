function  ptb_opto_eyelink_collection(  )
%ptb_opto_eyelink_collection Data collection after Optotrack & EyeLink initialized

%CONSTANTS  
OPTO = '_opto_';
TXT = '.txt';
DAT = '.dat'; 
FNRAWDATA = 'RawNDIDatFiles';
FNOTCREORG = 'OTCReorganized';
fixation_threshold = 75;
MISSINGDATACUTOFFVALUE = -3.6972e28; %an 'anything less than' cutoff value for missing data, 

run_experiment = true;
current_state = 1; % state machine behaviour
% 1 = waiting for unity (finger on homebutton)
% 2 = eyetracking, waiting for optotrack signal
% 3 = eyetracking + optotrack, waiting for end/exit  
% for 2 and 3, when eyetracking sees that eye out of range, it will
% stop recording for all devices and send the 'Restart' signal to Unity
time_between_eyelink_samples = 0.001;
countbad = 0; % count the number of times Eyelink has missed an eye;
error_threshold = 100; % if Eyelink misreads this many times (in a row), then quit
current_trial_name = '';
disp('Optotrack ready, switch to ShapeSpider now.');
optotrack_trial_mapping = [];
trial=1; %trial counter
while(run_experiment)
    % check is server connection still there
    [keyIsDown,secs,keyCode] = KbCheck;
    if keyCode(escKEY)
        current_state = 1;
    elseif keyCode(F3Key)
        fixation_pos = VisionLabEyelinkSetup();
    end
    switch current_state 
        case 1 % waiting for Unity to signal start
            %spool data from buffer to file
            countbad = 0;
            puRealtimeData=0;puSpoolComplete=0;puSpoolStatus=0;pulFramesBuffered=0;
            if(UnityPort.BytesAvailable)
                disp(['previous state: ' num2str(current_state)]);
                unityMessage = fscanf(UnityPort);
                unityMessage = strtrim(unityMessage);
                if ~isempty(strfind(unityMessage, 'Eyelink'))
                    block_txt = strsplit(unityMessage);
                    block_txt = block_txt{end};
                    current_trial_name = strsplit(unityMessage);
                    current_trial_name = current_trial_name(2);
                    % fetch the trial number also?
                    Eyelink('StartRecording'); % start recording
                    eyetrack_check = 0;
                    eyetrack_threshold = 1000;
                    while(~eyeOnFixation(fixation_pos, fixation_threshold, el, eye_used))
                        eyetrack_check = eyetrack_check + 1;
                        if(eyetrack_check > eyetrack_threshold)
                            fprintf(UnityPort,'Restart');
                            break;
                        end
                        error = Eyelink('CheckRecording');
                        if(error~=0)
                            disp('Eyelink Error!');
                            % attempt to restart Eyelink recording
                            Eyelink('StopRecording');
                            WaitSecs(0.005);
                            Eyelink('StartRecording');
                        end
                        WaitSecs(0.005);
                    end
                    disp('Sending Fixation');
                    fprintf(UnityPort,'Fixation');
                    current_state = 2;
                elseif strcmp(unityMessage,'Exit') % quit, maybe done experiment
                    run_experiment = false;
                end
                disp(['Unity sent: ' unityMessage]);
            end
        case 2
            % eyetrack, matlab waiting to start optotrack recording -- in a sense, MATLAB is in more control here
            error=Eyelink('CheckRecording');
            if(error~=0)
                disp('something goes wrong');
            end
            if(eyeOnFixation(fixation_pos, fixation_threshold,el,eye_used)) % eye is okay
                countbad = 0; % reset the counter for bad -- maybe we're good now?
                % check for message to start recording optotrack
                if(get(UnityPort,'BytesAvailable')>0)
                    disp(['previous state: ' num2str(current_state)]);
                    unityMessage = fscanf(UnityPort);
                    unityMessage = strtrim(unityMessage);
                    disp(['Unity sent: ' unityMessage]);
                    if strcmp(unityMessage,'Optotrack')
                        
                        %determine the trial number text
                        txtTrial = [num2str(trial, '%03d') '_']; % pad left side with zeros
                        
                        optotrack_trial_mapping = [optotrack_trial_mapping; [int2str(trial) current_trial_name]];
                        %navigate to the NDI Raw data file
                        cd([expDir '\' part_dems.id '\' FNRAWDATA]);

                        %The NDI dat file name for the current trial
                        NDIFName = [part_dems.id OPTO txtTrial block_txt DAT];

                        if trial==1
                            %re-initialize the reorganized data
                            dataReorg = zeros(fTotPerTrial, nMarkers*4);
                        end
% -------------- Start the optotrack recording here! -------------------- %
                        %start collecting data, the number of frames to collect was
                        %pre-specified by the function OPTOTRAKSETUPCOLLECTION.
                        
                        %initialize the file for spooling
                        OptotrakActivateMarkers();
                        WaitSecs(0.010);
                        DataBufferInitializeFile(0,NDIFName);
                        DataBufferStart();
                        current_state = 3;
                    elseif strcmp(unityMessage, 'Restart')
                        Eyelink('StopRecording');
                        current_state = 1;
                    elseif strcmp(unityMessage,'End') %Unity says to kill?!
                        Eyelink('StopRecording');
                        current_state = 1;
                    elseif strcmp(unityMessage,'Exit') % quit, for whatever reason
                        Eyelink('StopRecording');
                        current_state = 1; % this doesn't matter cause it won't loop again
                        run_experiment = false;
                    end
                    disp(['new state: ' num2str(current_state)]);
                end
            else % eye moved, let's start again
                countbad = countbad + 1;  
                if(countbad > error_threshold) % don't let one bad datapoint kill ya
                    Eyelink('StopRecording');
                    fprintf(UnityPort,'Restart');
                    disp('Eye issue');
                    current_state = 1;
                end
            end
            WaitSecs(time_between_eyelink_samples); % wait for next sample
        case 3 % eyetracking and  optotrack both on
            error=Eyelink('CheckRecording');
            if(error~=0)
                disp('something goes wrong');
            end
            if(UnityPort.BytesAvailable > 0)
                disp(['previous state: ' num2str(current_state)]);
                unityMessage = fgetl(UnityPort);
                unityMessage = strtrim(unityMessage);
                disp(['Unity sent: ' unityMessage]);
                if strcmp(unityMessage,'End')
% -------------------- stop optotrack here ---------------------------- %
                    Eyelink('StopRecording');
                    current_state = 4;
                elseif strcmp(unityMessage,'Exit') % quit, for whatever reason
% -------------------- stop optotrack here ---------------------------- %
                    Eyelink('StopRecording');
                    current_state = 4; % this doesn't matter cause it won't loop again
                    run_experiment = false;
                elseif strcmp(unityMessage,'Restart')
                    Eyelink('StopRecording');
                    current_state = 4;
                end
                disp(['new state: ' num2str(current_state)]);
            end
            if(~eyeOnFixation(fixation_pos,fixation_threshold,el,eye_used))
% -------------------- stop optotrack here ---------------------------- %
                countbad = countbad + 1;  
                if(countbad > error_threshold) % don't let one bad datapoint kill ya
                    Eyelink('StopRecording');
                    fprintf(UnityPort,'Restart');
                    disp('Eye issue');
                    current_state = 4;
                  end
            else
                countbad = 0;  
            end
            WaitSecs(time_between_eyelink_samples);
        case 4 % Completion of trial
% -------------------- spooling happens here -------------------------- %
            Eyelink('StopRecording');
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
            
            %open the .dat file stored by Optotrak
            [fid, message] = fopen(deblank(NDIFName),'r+','l'); % little endian byte ordering

            %Read the header portion of the file
            [filetype,items,subitems,numframes,frequency,user_comments,sys_comments,descrip_file,cutoff,coll_time,coll_date,frame_start,extended_header,character_subitems,integer_subitems,double_subitems,item_size,padding] = read_nd_c_file_header(fid);

            %Read the data portion of the file
            rawData = read_nd_c_file_data(fid, items, subitems, numframes);
            fclose(fid);
            disp('converting data');
            %Convert the data to a format OTCTextWrite accepts
            for i = 1:nMarkers
    
                %append the marker data to dataReorg and delete the 
                dataReorg(:,((i+1)*3)-2:(i+1)*3) = transpose(squeeze(rawData(i,:,:)));
            end
    
            %Convert missing data to NaNs
            dataReorg(dataReorg < MISSINGDATACUTOFFVALUE) = NaN;
        
            %add the data collection information to the data
            dataReorg(:,1) = trial;
            dataReorg(:,2) = smpRt;
        
            %navigate to the 'OTCReorganized' file folder
            cd([expDir '\' part_dems.id '\' FNOTCREORG]);

            %ALWAYS STORE A RAW DATA SET
            disp(['Saving file: ' part_dems.id OPTO txtTrial block_txt TXT '...'])
            fid = fopen([part_dems.id OPTO txtTrial block_txt TXT], 'w');
            %fprintf(fid, hdr);
            fclose(fid);
            dlmwrite([part_dems.id OPTO txtTrial block_txt TXT],dataReorg,'-append','delimiter',...
                '\t','newline','pc','precision',12);
            disp('Success!')

            current_state = 1;
            trial=trial+1; %advance to the next trial
    end
end

disp(optotrack_trial_mapping);

WaitSecs(0.1);
Eyelink('CloseFile');
% download data file
try
    fprintf('Receiving data file ''%s''\n', edfFile );
    status=Eyelink('ReceiveFile');
    if status > 0
        fprintf('ReceiveFile status %d\n', status);
    end
    if 2==exist (edfFile, 'file')
        fprintf('Data file ''%s'' can be found in ''%s''\n', edfFile, pwd );
    end
catch
    fprintf('Problem receiving data file ''%s''\n', edfFile );
end
Eyelink('ShutDown');
Screen('CloseAll');

fclose(t);
delete(t);


fclose(UnityPort);
cd(olddir);
end

