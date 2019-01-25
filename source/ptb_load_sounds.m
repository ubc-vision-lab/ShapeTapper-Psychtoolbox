function [ snd_buffers, fs ] = ptb_load_sounds(snd_dir, snd_names, snd_formats)
%ptb_loadtextures % Load each stimulus image as a PsychToolbox texture

% store handles in a Map object
snd_buffers = containers.Map;

if isempty(snd_names)
    return;
end

for i=1:length(snd_names)
    % Retrieve name of stimulus image
    snd_name = char(snd_names(i));
    snd_fname = [];
    
    for sf = 1:length(snd_formats)
        sn = [snd_dir snd_name char(snd_formats{sf})];
        if exist(sn, 'file') == 2
            snd_fname = sn;
        end
    end
    
    if isempty(snd_fname)
        sca;
        error(['Could not load stimulus sound \"%s\" \n' ...
               'Please check your config file and stimulus directory.\n'...
               'Exiting...'], snd_name);
    end
    
    % Load sound using wavread(), save sample rate for further processing
    [snd_data, fs] = audioread(snd_fname);
    
    % Add sound name and buffer handle to Map for later playback
    snd_buffers(snd_name) = snd_data';
end

end

