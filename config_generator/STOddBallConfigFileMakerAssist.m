function [A,B,C,D,E,F,G,H] = STOddBallConfigFileMakerAssist

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
        
%format specifier for saving the data
C = ['%u,%u,%u,%f,%u,%f,%u,%u,%s,%u,%s,%s,%s,%s,%f,%f,%f,%f,%f,%f,%f,%u,%u,%u,%u,%f,%f,%f,'...
    '%f,%f,%f,%f,%u,%f\r\n'];

%format specifier for reading-in the data.
D = '%u%u%u%f%u%f%u%u%s%u%s%s%s%s%f%f%f%f%f%f%f%u%u%u%u%f%f%f%f%f%f%f%u%f';


%Field names for the condition file
E{1,1} = 'Trial #'; %u
E{2,1} = 'Block #'; %u
E{3,1} = 'Trial # w/n Block'; %u
E{4,1} = 'Target Img'; %s
E{5,1} = 'X-Pos. Target Img'; %f
E{6,1} = 'Y-Pos. Target Img'; %f
E{7,1} = 'Rotation Target Img'; %f
E{8,1} = 'Img1 Distractor'; %s
E{9,1} = 'X-Pos. Img1 Distractor'; %f
E{10,1} = 'Y-Pos. Img1 Distractor'; %f
E{11,1} = 'Rotation Img1 Distractor'; %f
E{12,1} = 'Img2 Distractor'; %s
E{13,1} = 'X-Pos. Img2 Distractor'; %f
E{14,1} = 'Y-Pos. Img2 Distractor'; %f
E{15,1} = 'Rotation Img2 Distractor'; %f
E{16,1} = 'Img Onset (s)'; %f
E{17,1} = 'Img Duration (s)'; %f

%Field Names (column headers to be saved to the file) for the condition order file
F = ['Trial #,Block #,Trial # w/n Block,Target Img,X-Pos. Target Img,Y-Pos. Target Img,'...
        'Rotation Target Img,Img1 Distractor,X-Pos. Img1 Distractor,Y-Pos. Img1 Distractor,'...
        'Rotation Img1 Distractor,Img2 Distractor,X-Pos. Img2 Distractor,'...
        'Y-Pos. Img2 Distractor,Rotation Img2 Distractor,Img Onset (s),Img Duration (s)\r\n'];

%the file specifier for saving a condition order file (character array)
G = '%u,%u,%u,%s,%f,%f,%f,%s,%f,%f,%f,%s,%f,%f,%f,%f,%f\r\n';

%the file specifier for reading-in a condition order file (character array)
H = '%u%u%u%s%f%f%f%s%f%f%f%s%f%f%f%f%f';

return;