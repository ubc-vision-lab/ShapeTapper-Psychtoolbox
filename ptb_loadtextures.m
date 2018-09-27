function [ stim_textures ] = ptb_loadtextures(stim_dir, img_names, img_formats, ptb)
%ptb_loadtextures % Load each stimulus image as a PsychToolbox texture

% store handles in a Map object
stim_textures = containers.Map;

for i=1:length(img_names)
    % Retrieve name of stimulus image
    img_name = char(img_names(i));
    img_fname = [];
    
    for ff = 1:length(img_formats)
        fn = [stim_dir img_name char(img_formats{ff})];
        if exist(fn, 'file') == 2
            img_fname = fn;
        end
    end
    
    if isempty(img_fname)
        sca;
        error(['Could not load stimulus image \"%s\" \n' ...
               'Please check your config file and stimulus directory.\n'...
               'Exiting...'], img_name);
    end
    
    % Load image using imread(), save colormap and alpha channel
    % for further processing if needed
    [img, map, alpha] = imread(img_fname);
   
    % Create white background (may be changed to another value if needed)
    white_bg = (ones(size(img))*255);
    
    if ~isempty(alpha)
        % Scale alpha channel and apply background to transparent part of image
        alphaMask = im2double(alpha); %scale between 0 and 1
        img = im2uint8(white_bg.*repmat((1-alphaMask),[1 1 3]) + ...
                       double(img).*repmat(alphaMask,[1 1 3]));
    end
                   
    % Make the image into a texture, save object handle as property
    h = Screen('MakeTexture', ptb.window, img);
    
    % Add image name and texture handle to Map for later display
    stim_textures(img_name) = h;
end

end

