# ShapeTapper-Psychtoolbox
Shape Tapper Experiment Framework - Psychtoolbox Version

# Requirements:
* For Optotrack compatability:
  1. MATLAB 2015b
  2. Psychtoolbox 3.0.11 beta
  3. Subversion

* For use without Optotrack:
  1. MATLAB 2016a or later
  2. Psychtoolbox 3.0.14 beta or later
  3. Subversion

For Psychtoolbox Installation Instructions see here:
http://psychtoolbox.org/download.html

# Instructions:
1. Ensure that stimulus images are in the `Image_Files/` folder. Images must be in `.png` or `.jpg` formats.
2. Use the `STConfigFileMaker.m` or variant to generate an experiment config file if not already present. Place config files in the `Config_Files/` subfolder.
3. Open `ShapeTapper_RunExperiment.m` in MATLAB and run the script, entering participant demographic info and selecting the appropriate config file where necessary.
4. Data will be saved to the `Data/` subfolder. For incomplete experiments, progress is automatically saved on a trial-by-trial basis. There is an option to resume experiment by re-entering the participant ID and selecting the same config file during script startup.

# Note on Touchscreens
Psychtoolbox has poor support for touchscreen drivers (according to their developers this issue is solved in Linux installations). There is often an issue having Psychtoolbox detect touches unless the touchscreen has been configured to act as a mouse.

Using the ELO Touch 2293L, we were able to bypass this issue by re-installing the EloMultiTouch drivers with the `ForceMouse = 1` option set in `\Common\EloOptions.ini` installation file. This forces Windows to recognize touches as mouse clicks, rather than a proprietary signal that Psychtoolbox cannot parse.

See the `Elo Touch Solutions Multi-Touch Driver Package User Manual` for more details.
