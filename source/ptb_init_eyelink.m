function el = ptb_init_eyelink(window, out_path, fname, part_dems)
%init_eyelink Initialize EyeLink and perform calibration

% Eyelink Setup routine
if EyelinkInit()~= 1 %
    el = struct();
    return;
end

% Initialize EyeLink defaults with Psychtoolbox window
el=EyelinkInitDefaults(window);

Eyelink('Command', 'link_sample_data = LEFT,RIGHT,GAZE,AREA');

% Get Eyelink Tracker Version - NOT USED
% [v vs]=Eyelink('GetTrackerVersion');

% Create EDF file (EyeLink recording file)
el.edfFile = [part_dems.id '.edf'];
Eyelink('Openfile', el.edfFile);

EyelinkDoTrackerSetup(el); % do EyeLink calibration, validation, etc.

end

