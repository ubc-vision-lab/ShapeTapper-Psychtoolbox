function [ptb, bg_color_change] = ptb_set_bg_color(stim_bg_color, stim_text_color, ptb)
%PTB_SET_BG_COLOR Changes background and text display color, if necessary
%
%   Input: ptb = psychtoolbox struct returned by ptb_initscreen.m
%          stim_bg_color = background color specified in config file row
%          stim_text_color = text color specified in config file row
%
%   Output: ptb = updated struct with current background & text colors
%           bg_color_change = bool, 0 if no change occured, 1 if change


% Initialize output var
bg_color_change = false;

% Convert stimulus row colors from RGB strings to vectors
new_bg_color = cellfun(@str2double, strsplit(stim_bg_color));
new_text_color = cellfun(@str2double, strsplit(stim_text_color));

% Scale RGB values (assumed [0-1]) to PTB screen black and white range
bw_scale = ptb.white - ptb.black;
new_bg_color = (new_bg_color .* bw_scale) + ptb.black;
new_text_color =  (new_text_color .* bw_scale) + ptb.black;

% If a change in background color is detected, then redraw background
if ptb.bg_color ~= new_bg_color
    ptb.bg_color = new_bg_color;
    Screen('FillRect', ptb.window, ptb.bg_color);
    bg_color_change = true;
end

% If a change in text color is detected, then save text color for output
if ptb.text_color ~= new_text_color
    ptb.text_color = new_text_color;
    bg_color_change = true;
end

end

