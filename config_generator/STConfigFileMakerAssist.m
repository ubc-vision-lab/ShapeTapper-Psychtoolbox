function [A,B,C,D] = STConfigFileMakerAssist

A = {}; %the field names (cell array)
B = []; %the header line (character array)
C = []; %the file specifier for saving a config file (character array)
D = []; %the file specifier for loading in a config file (character array)

%Field Names (column headers to be saved to the file.
A{1,1} = 'block_num'; %u
A{2,1} = 'is_practice_block'; %u
A{3,1} = 'block_has_feedback'; %u
A{4,1} = 'block_feedback_thresh'; %f
A{5,1} = 'trial_num'; %u
A{6,1} = 'trial_max_time'; %f
A{7,1} = 'trial_feedback'; %u
A{8,1} = 'background_color'; %s
A{9,1} = 'stim_img_name'; %s
A{10,1} = 'stim_onset'; %f
A{11,1} = 'stim_duration'; %f
A{12,1} = 'stim_cent_x'; %f
A{13,1} = 'stim_cent_y'; %f
A{14,1} = 'stim_size_x'; %f
A{15,1} = 'stim_size_y'; %f
A{16,1} = 'stim_rotation'; %f
A{17,1} = 'stim_is_touchable'; %u
A{18,1} = 'stim_is_target'; %u
A{19,1} = 'subj_fixation_type'; %u
A{20,1} = 'subj_fixation_onset'; %f
A{21,1} = 'subj_fixation_duration'; %f
A{22,1} = 'mask_onset'; %f
A{23,1} = 'mask_duration'; %f
A{24,1} = 'mask_size'; %f
A{25,1} = 'mask_color'; %f
A{26,1} = 'mask_rotation'; %f
A{27,1} = 'mask_fit'; %u
A{28,1} = 'mask_margin'; %f

%Create the header line 
for i = 1:length(A)
    if i < length(A)
        B = [B A{i,1} ','];
    else
        B = [B A{i,1} '\r\n'];
    end
end
        
%format specifier for saving the data
C = '%u,%u,%u,%f,%u,%f,%u,%s,%s,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%f,%f,%f,%f,%f,%f,%f,%u,%f\r\n';

%format specifier for reading the data in.
D = '%u%u%u%f%u%f%u%s%s%f%f%f%f%f%f%f%u%u%u%f%f%f%f%f%f%f%u%f';