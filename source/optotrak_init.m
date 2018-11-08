 function optotrak_init()
    %UNTITLED1 Summary of this function goes here
    %   Detailed explanation goes here
    
    %Ensure that the appropriate API files (oapi.dll and ndopto.h) are
    %accessible on the Matlab path
    addpath('C:\ndigital\programs\', 'C:\NDIoapi\ndlib\include');
    
    %load the library into the matlab environment
    
    %test if the library loaded
    libloaded = libisloaded('oapi');
    
    %if it isn't loaded, load it
    if libloaded~=1
        % Original Call: If you HAVE NEVER loaded this library into Matlab
        % before, you first need to PROTOTYPE the header file.  The
        % following call will do that for you, saving the prototype file as ndoptoPrototype.m:
        %%loadlibrary('C:\ndigital\programs\oapi.dll', 'C:\NDIoapi\ndlib\include\ndopto.h', 'addheader','C:\NDIoapi\ndlib\include\ndtypes.h', 'mfilename','ndoptoPrototype');
        
        %IMPORTANT: If you are using Matlab version 2008 or newer, you need
        %to CHANGE the prototype file by finding all instances of 'cdecl' 
        %and changing it to 'stdcall', explained here:
        %http://www.mathworks.com/support/solutions/en/data/1-671ZZL/index.html?product=ML&solution=1-671ZZL
 
        
        % Once you have the prototype file, you load the library using this
        % function call:
        loadlibrary('C:\ndigital\programs\oapi.dll', @ndoptoPrototype)
    end
    
    %test if the library loaded
    libloaded = libisloaded('oapi');
    if libloaded~=1
        error('oapi.dll did not load');
    end
    
    %Call if cabling configuration of system has changed (anything
    %unplugged or reconnected since last use).  Creates camera file.
    %see API manual for more details
    TransputerDetermineSystemCfg('sysconf.log');
    
    %Dowloads code about system configuration.  Must call Initialize system
    %after this function - only need to call once after power-up. See API
    %manual for more details
    TransputerLoadSystem('system');
    pause(1);
    
    %Establish communication with Optotrak.  Must be called after
    %LoadSystem.  See API for more details
    TransputerInitializeSystem(1);
    
    %Read and send camera file info to Optotrak.  Camera file name is
    %specified in the function call.  Search path is the current directory,
    %and the the /realtime folder.
    OptotrakLoadCameraParameters('standard');

    %get the status, returns information about all components of your
    %Optotrak setup.
    pnNumSensors=0;
    pnNumOdaus=0;
    pnNumRigidBodies=0;
    pnMarkers=0;
    pfFrameFrequency=0;
    pfMarkerFrequency=0;
    pnThreshold=0;
    pnMinimumGain=0;
    pnStreamData=0;
    pfDutyCycle=0;
    pfVoltage=0;
    pfCollectionTime=0;
    pfPreTriggerTime=0;
    pnFlags=0;
    [pnNumSensors,pnNumOdaus,pnNumRigidBodies,pnMarkers,pfFrameFrequency,pfMarkerFrequency,pnThreshold,pnMinimumGain,pnStreamData,pfDutyCycle,pfVoltage,pfCollectionTime,pfPreTriggerTime,pnFlags]=OptotrakGetStatus(pnNumSensors,pnNumOdaus,pnNumRigidBodies,pnMarkers,pfFrameFrequency,pfMarkerFrequency,pnThreshold,pnMinimumGain,pnStreamData,pfDutyCycle,pfVoltage,pfCollectionTime,pfPreTriggerTime,pnFlags);

 end
