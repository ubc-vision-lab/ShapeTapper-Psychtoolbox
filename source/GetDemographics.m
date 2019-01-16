function part_dems = GetDemographics()
%GetDemographics Get participant demographics (handedness, gender, age)

% Init output structure
part_dems = struct();

% Generate subject ID
symbols = ['a':'z' 'A':'Z' '0':'9'];
MAX_ST_LENGTH = 6;
nums = randi(numel(symbols),[1 MAX_ST_LENGTH]);
subj_id = symbols (nums);
part_id = inputdlg('Participant ID (press Ok to continue):','Participant ID',[1 40],{subj_id});

if isempty(part_id)
    return
end

% Get participant handedness
handedness = {'Left','Right'};
[part_h,tf] = listdlg('PromptString','Select Participant''s Handedness:',...
                           'SelectionMode','single',...
                           'ListString',handedness);
if isempty(part_h) || ~tf
    return
end

% Get participant gender
genders = {'Male','Female'};
[part_g,tf] = listdlg('PromptString','Select Participant''s Gender:',...
   'SelectionMode','single',...
   'ListString',{'Male','Female'});

if isempty(part_g) || ~tf
    return
end

% Get participant age
part_age_resp = inputdlg('Enter Participant''s Age:','Age',[1 30],{'99'});

if isempty(part_age_resp)
    return
end

part_dems.id = part_id{1};
part_dems.handedness = handedness{part_h};
part_dems.gender = genders{part_g};
part_dems.age = part_age_resp{1};

end

