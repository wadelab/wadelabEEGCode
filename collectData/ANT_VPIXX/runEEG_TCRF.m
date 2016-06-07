% This version of the code sweeps over tf and contrast space with no mask.
% It is originally used for fly ERG however this version is used for human
% EEG.
% ARW June 11 2014
% ARW Feb 2016: This version sweeps over contrast and TF (not SF)
% MMH May 5th: Changes for EEG
%


close all;
clear all;
startTime=tic;
DUMMYRUN=1; % Do we want to run without displaying anything (useful for testing)
DPIXXPRESENT = 0; % Are we physically connected to a DPIXX?

if (~DUMMYRUN)
    
    
    Screen('Preference', 'VisualDebuglevel', 0)% disables welcome and warning screens
    HideCursor % Hides the mouse cursor
    
    
    % Get the calibration and compute the gamma table
    igt=fly_computeInverseGammaFromCalibFile('CalibrationVPIxxFake.mat')
    dpy.gamma.inverse=igt;
end
if DPIXXPRESENT        % using a ViewPixx or ProPixx
    try
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        Datapixx('Open');
        Datapixx('DisableVideoScanningBacklight');      % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');       % Only required if a VIEWPixx.
        
        Datapixx('RegWr'); % Make sure the commands above are executed
    catch
        lasterr
        error('Error using Datapixx - are you sure you are connected to one?');
    end
    
end




datadir='.';
EEG_startTime=now;

%we have to edit this based on display
dpy.res = [1920 1080]; % screen resoloution
dpy.size = [.53 .3] % Meters
dpy.distance = [.22]; % Metres back from screen
dpy.frameRate=144;
% dpy will eventually contain all the info about the display e.g. size,
% distance, refresh rate, spectra, gamma.
% For now if just has the gamma function (inverse) in it.

tfList=[1,8,36;1,8,36]'; % This is in Hz.
sfList=[.056]'; % cpd
contList=[2 16 64;0 0 0]'/100;  % Contrasts to probe
nTF=size(tfList,1);% number of TFs
nSF=size(sfList,1);
nCont=size(contList,1);%number of Contrasts


ordered=1:(nTF*nCont);

% This fully shuffles the order
shuffleSeq=Shuffle(ordered); % Shuffle all the possible presentation conditions


stim.spatial.internalRotation = 0; % Does the grating rotate within the envelope?
stim.rotateMode = []; % rotation of mask grating (1= horizontal, 2= vertical, etc?)

stim.spatial.angle = [0 0]  ; % angle of gratings on screen
stim.temporal.duration=11; % how long to flicker for - remember we remove first second

% In principle we can have both mask and probe here. For now though we set
% the mask to zero. probeCont isn't used in this code (yet)
probeCont=[99]/100;
maskCont =[0];

nConds=length(ordered); %number of total conditions

nRepeats=1; %how many times to repeat each condition

% Initialize the screen just once. Store all the parameters in the dpy
% structure


%Now we run the Plaid
if (~DUMMYRUN)
Screen('Preference', 'SkipSyncTests', 0);

% Select Screen, 0 for main
WhichScreen = 0; %
Screen('Preference', 'VisualDebuglevel', 1)% disables welcome and warning screens
HideCursor % Hides the mouse cursor

% Open a fullscreen onscreen window on that display, choose a background
% color of 128 = gray, i.e. 50% max intensity:
dpy.win = Screen('OpenWindow', WhichScreen, 128);

% Set the gamma tables.
% Set the CLUTS to the calibrated values
oldClut = LoadIdentityClut(dpy.win, 1);
disp('Loaded identify clut ******');


% Make sure the GLSL shading language is supported:
AssertGLSL;

% Retrieve video redraw interval for later control of our animation timing:
dpy.ifi = Screen('GetFlipInterval', dpy.win);
 % Compute increment of phase shift per redraw:
dpy.phaseincrement = [stim.temporal.frequency] * 360 * dpy.ifi;
        
        
else
    disp('Dummy run - not initializing the screen');
    
    dpy.win=-1;
end

tic
for thisRun=1:nRepeats  % 3 repeats
    % Write something to the DataPixx to say that we are about to start the
    % sequence
    % Write out a trigger to the DataPixx if required
    if DPIXXPRESENT
        for n = 1:3  % 3 trigger pulses of '1' to indicate start of experiment
            Datapixx('SetDoutValues', transformindex(thisRun));
            Datapixx('RegWrRd');
            Datapixx('SetDoutValues', 0);
            Datapixx('RegWrRd');
            
        end
    else
        fprintf('\nCurrent condition trigger code %d\n',thisRun);
    end
    
    
    
    for thisCond=1:nConds
        stim.spatial.frequency=[sfList(1), sfList(1)]; % Set the same spatial freq for both components
        % Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
        stim.spatial.phase=[0 0 ]; %[rand(1)*360 rand(1)*360];
        stim.spatial.pOffset=rand(2,1)*360;
        
        
        thisaction= shuffleSeq(thisCond);
        tIndex=ceil(thisaction/ nTF); % Index into the list of temporal frequencies
        cIndex=1+rem(thisaction, nCont); %  Should be 1+rem(thisAction-1,nTF)  - index into the list of contrasts
        
        thisTemp(thisRun,thisCond)=tIndex;
        thisCont(thisRun,thisCond)=cIndex;
        
        stim.temporal.frequency=tfList(tIndex,:);
        stim.cont=contList(cIndex,:);
        disp(thisRun)
        disp(thisCond)
        
        fprintf('\nRunning %.2d %.2d',stim.cont(1),stim.temporal.frequency(1));
        
       
        % Write out a trigger to the DataPixx if required
        if DPIXXPRESENT
            for n = 1:3                         % 3 trigger pulses to indicate start of block
                Datapixx('SetDoutValues', transformindex(thisaction+100));
                Datapixx('RegWrRd');
                
                Datapixx('SetDoutValues', 0);
                Datapixx('RegWrRd');
                
            end
        else
            fprintf('\nCurrent condition trigger code %d\n',thisaction+100);
        end
        
        % Call the main display loop - the screen is already open so this
        % should be fast
        if (~DUMMYRUN) % Do we display something or not?
            d=EEG_runPlaid(dpy,stim);
        end % End check on dummy run
        
        if DPIXXPRESENT
            for n = 1:3                         % 3 trigger pulses to indicate end of block
                Datapixx('SetDoutValues', transformindex(255));
                Datapixx('RegWrRd');
                
                Datapixx('SetDoutValues', 0);
                Datapixx('RegWrRd');
                
            end
        else
            fprintf('\nCurrent condition trigger code %d\n',255);
        end
        
        
    end % Next contrast and temporal frequency pair
end % Next repetition


% We're done. Close the window. This will also release all other ressources:
% Bye bye!
totalSessionTime=toc;

if (~DUMMYRUN) % Save as soon as possible in case it crashes
    filename=fullfile(datadir,['EEG_CRF',datestr(EEG_startTime,30),'.mat'])
    save(filename)
end

if DPIXXPRESENT       % Close the DataPixx
    Datapixx('Close');
    Datapixx('RegWr'); % Make sure the commands above are executed
end

Screen('CloseAll');
sca
pause(2);






