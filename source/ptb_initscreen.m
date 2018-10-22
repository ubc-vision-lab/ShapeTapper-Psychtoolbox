function [ ptb ] = ptb_initscreen( bg_color )
%ptb_initscreen Inits Psychtoolbox and returns struct with window vars

%----------------------------------------------------------------------
%                       PyschToolbox Setup
%----------------------------------------------------------------------

ptb = struct();

%bypass psychtoolbox sync tests
Screen('Preference', 'SkipSyncTests', 2 );

% Reduce error messages to increase frame rate
Screen('Preference', 'Verbosity', 1);

% Setup PTB with some default values
PsychDefaultSetup(2);

% Seed the random number generator. Here we use the an older way to be
% compatible with older systems. Newer syntax would be rng('shuffle'). 
% Look at the help function of rand "help rand" for more information
rand('seed', sum(100 * clock));


%----------------------------------------------------------------------
%                       Display Information
%----------------------------------------------------------------------

% Get Pixels per Inch of current monitor, convert to centimeters
ppi = get(groot,'ScreenPixelsPerInch');
ptb.ppcm = ppi / 2.54;

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer. For help see: Screen Screens?
screens = Screen('Screens');

% Draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen. When only one screen is attached to the monitor we will draw to
% this. For help see: help max
ptb.screenNumber = max(screens);

% Hide cursor for touchscreen display (will use ShowCursor on program exit)
HideCursor(ptb.screenNumber);

% Define black and white (white will be 1 and black 0). This is because
% luminace values are (in general) defined between 0 and 1. For help see:
% help WhiteIndex and help BlackIndex
ptb.white = WhiteIndex(ptb.screenNumber);
ptb.black = BlackIndex(ptb.screenNumber);

% Test if bg_color has been defined as 'black' or 'white', or if an rgb
% array has been specified. If not, then default to white.
if ischar(bg_color{:})
	if strcmp(bg_color,'white')
        ptb.bg = ptb.white;
    elseif strcmp(bg_color,'black')
        ptb.bg = ptb.black;
    else
        ptb.bg = ptb.white;
    end
elseif isvector(bg_color) && length(bg_color) == 3
    ptb.bg = bg_color;
else
    ptb.bg = ptb.white;
end

% Open an on screen window and color the background For help see: Screen
% OpenWindow?
[ptb.window, ptb.windowRect] = PsychImaging('OpenWindow', ...
                                              ptb.screenNumber, ptb.bg);

% Retreive the maximum priority number and set max priority
topPriorityLevel = MaxPriority(ptb.window);
Priority(topPriorityLevel);

% Get the size of the on screen window in pixels For help see: Screen
% WindowSize?
[ptb.screenXpixels, ptb.screenYpixels] = Screen('WindowSize', ptb.window);

% Get the centre coordinate of the window in pixels For help see: help
% RectCenter
[ptb.xCenter, ptb.yCenter] = RectCenter(ptb.windowRect);

% Enable alpha blending for anti-aliasing For help see: Screen
% BlendFunction? Also see: Chapter 6 of the OpenGL programming guide
Screen('BlendFunction', ptb.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Set the text size
Screen('TextSize', ptb.window, 40);

%----------------------------------------------------------------------
%                       Mouse & Keyboard information
%----------------------------------------------------------------------

% Use OSX internal naming scheme, to increase portability of this script
KbName('UnifyKeyNames');

% Define the keyboard keys that are listened for.
ptb.escapeKey = KbName('ESCAPE');
ptb.spaceKey = KbName('SPACE');
ptb.pKey = KbName('p');
ptb.F3Key = KbName('F3');
ptb.F12Key = KbName('F12');

% % Here we set the initial position of the mouse to a random position on
% the % screen SetMouse(round(rand * screenXpixels), round(rand *
% screenYpixels), window);


%----------------------------------------------------------------------
%                       Timing Information
%----------------------------------------------------------------------

% Measure the vertical refresh rate of the monitor
ptb.ifi = Screen('GetFlipInterval', ptb.window);

% Numer of monitor frames to wait between screen renders (default: 1)
ptb.waitframes = 1;


end

