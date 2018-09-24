function [A,B,C] = STConfigFileMakerAssist

A = {};
B = [];
C = [];

%Headers (column headers to be saved to the file.
A{1,1} = 'block_num'; %u
A{2,1} = 'is_practice_block'; %u
A{3,1} = 'block_has_feedback'; %u
A{4,1} = 'block_feedback_thresh'; %f
A{5,1} = 'trial_num'; %u
A{6,1} = 'trial_max_time'; %f
A{7,1} = 'trial_feedback'; %u
A{8,1} = 'stim_img_name'; %s
A{9,1} = 'stim_onset'; %f
A{10,1} = 'stim_duration'; %f
A{11,1} = 'stim_cent_x'; %f
A{12,1} = 'stim_cent_y'; %f
A{13,1} = 'stim_size_x'; %f
A{14,1} = 'stim_size_y'; %f
A{15,1} = 'stim_rotation'; %f
A{16,1} = 'stim_is_target'; %u
A{17,1} = 'mask_onset'; %f
A{18,1} = 'mask_duration'; %f
A{19,1} = 'mask_size'; %f
A{20,1} = 'mask_color'; %f
A{21,1} = 'mask_rotation'; %f
A{22,1} = 'mask_fit'; %f
A{23,1} = 'mask_margin'; %f
A{24,1} = 'mask_touchable'; %u

%Create the header line 
for i = 1:length(A)
    if i < length(A)
        B = [B A{i,1} ','];
    else
        B = [B A{i,1} '\r\n'];
    end
end
        
%format specifier for the data
C = '%u,%u,%u,%f,%u,%f,%u,%s,%f,%f,%f,%f,%f,%f,%f,%u,%f,%f,%f,%f,%f,%f,%f,%u\r\n';