function [A,B,C,D,E,F,G,H] = SpiderStudyConfigFileMakerAssist

A = {}; %the field names for the config file itself (cell array)
B = []; %the header line (character array)
C = []; %the file specifier for saving a config file (character array)
D = []; %the file specifier for loading in a config file (character array)

E = {}; %the field names for condition file (cell array)
F = []; %the header-line for the condition order file (character array)
G = []; %the file specifier for saving a condition order file (character array)
H = []; %the file specifier for reading-in a condition order file (character array)

%Field Names (column headers to be saved to the file.
A{1,1} = 'block_num'; %u
A{2,1} = 'is_practice_block'; %u
A{3,1} = 'block_has_feedback'; %u
A{4,1} = 'block_feedback_thresh'; %f
A{5,1} = 'trial_num'; %u
A{6,1} = 'trial_max_time'; %f
A{7,1} = 'trial_timeout_msg'; %u
A{8,1} = 'trial_kb_resp'; %u
A{9,1} = 'correct_kb_resp'; %s
A{10,1} = 'trial_feedback_type'; %u
A{11,1} = 'trial_feedback_img'; %s
A{12,1} = 'background_color'; %s
A{13,1} = 'text_color'; %s
A{14,1} = 'stim_img_name'; %s
A{15,1} = 'stim_onset'; %f
A{16,1} = 'stim_duration'; %f
A{17,1} = 'stim_cent_x'; %f
A{18,1} = 'stim_cent_y'; %f
A{19,1} = 'stim_size_x'; %f
A{20,1} = 'stim_size_y'; %f
A{21,1} = 'stim_rotation'; %f
A{22,1} = 'stim_is_touchable'; %u
A{23,1} = 'stim_is_target'; %u
A{24,1} = 'subj_fixation_type'; %u
A{25,1} = 'subj_fixation_pause'; %u
A{26,1} = 'subj_fixation_onset'; %f
A{27,1} = 'subj_fixation_duration'; %f
A{28,1} = 'mask_onset'; %f
A{29,1} = 'mask_duration'; %f
A{30,1} = 'mask_size'; %f
A{31,1} = 'mask_color'; %f
A{32,1} = 'mask_rotation'; %f
A{33,1} = 'mask_fit'; %u
A{34,1} = 'mask_margin'; %f

%Create the header line 
for i = 1:length(A)
    if i < length(A)
        B = [B A{i,1} ','];
    else
        B = [B A{i,1} '\r\n'];
    end
end
        
%format specifier for saving the config file
C = ['%u,%u,%u,%f,%u,%f,%u,%u,%s,%u,%s,%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%f,%f,%f,'...
    '%f,%f,%f,%f,%u,%f\r\n'];

%format specifier for reading-in the config file
D = '%u%u%u%f%u%f%u%u%s%u%s%s%s%s%f%f%f%f%f%f%f%u%u%u%u%f%f%f%f%f%f%f%u%f';

%Field Names (column headers to be saved to the file.
E{1,1} = 'Trial #'; %u
E{2,1} = 'Block #'; %u
E{3,1} = 'Trial # w/n Block'; %u
E{4,1} = 'Jump?'; %u
E{5,1} = 'Img1 Val.'; %s
E{6,1} = 'Img2 Val.'; %s
E{7,1} = 'Valence Switch?'; %u
E{8,1} = 'Initial Position'; %u
E{9,1} = 'Final Position'; %u
E{10,1} = 'Initial Img'; %s
E{11,1} = 'Final Img'; %s
E{12,1} = 'Target Letter'; %s

%Field Names (column headers to be saved to the file) for the condition order file
F = ['Trial #,Block #,Trial # w/n Block,Jump?,Img1 Val.,Img2 Val.,Valence Switch?,'...
    'Initial Position,Final Position,Initial Img,Final Img,Target Letter\r\n'];

%the file specifier for saving a condition order file (character array)
G = '%u,%u,%u,%u,%s,%s,%u,%u,%u,%s,%s,%s\r\n';

%the file specifier for reading-in a condition order file (character array)
H = '%u%u%u%u%s%s%u%u%u%s%s%s';

return;
