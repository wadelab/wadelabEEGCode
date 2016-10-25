% This version of the code sweeps over tf and contrast space with no
% mask.qqq
% It is originally used for fly ERG however this version is used for human
% EEG.        
% ARW June 11 2014
% ARW Feb 2016: This version sweeps over contrast and TF (not SF)
% MMH May 5th: Changes for EEG incl triggers, fixation, isi.
% MMH August 1st: edited so we can see which tf/cont cond trigger is being
%sent out as we w ere getting mixed up contrasts - was starting at 2nd in
%order rather than 1st.
%MMH September 9th: Added intro and annulus  
%MMH OCT 12: We can now set inner and outer radius of a mask
%(mask.radiusDeg) to split the stimuli with an annulus. We also randomise
%the grating orientation each trial (0 or 90), and we have full randomisation of pres order.

close all
clear all;
startTime=tic;
DUMMYRUN=0; % Do we want to run without displaying anything (useful for testing)
DPIXXPRESENT=1; % Are we physically connected to a DPIXX?
WhichScreen = 0;

%% get session details and create associated save file for results
% R = input('Enter participant details\n');
% C = input('Enter Condition\n');
% myfile = sprintf('%s_%s_EEG_TCRF_%s',R,C,datestr(now,'ddmmyy'));
% [~,mydir] = uiputfile(myfile,'Choose file directory');
% out_file = [mydir,myfile,'.mat'];

if (~DUMMYRUN)
        
    Screen('Preference', 'VisualDebuglevel', 0)% disables welcome and warning screens
    HideCursor % Hides the mouse cursor   
end

if DPIXXPRESENT        % using a ViewPixx or ProPixx
    try
        AssertOpenGL;
        % Open PTB onscreen window: We request a 32 bit per colour component
        % floating point framebuffer if it supports alpha-blending. Otherwise
        % the system shall fall back to a 16 bit per colour component framebuffer:
        PsychImaging('PrepareConfiguration');
        PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible')
        
        % required for gamma correction through the PsychImaging pipeline:
        PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');
       
        % Open an on screen (grey) window and configure the imaging pipeline
        [dpy.win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);

        PsychImaging('PrepareConfiguration');
        
        % Initialise the Vpixx device:
        PsychImaging('AddTask', 'General', 'UseDataPixx');
        dataPixxOk=Datapixx('Open');
        % The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
        Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
        Datapixx('EnableVideoStereoBlueline');
        
        %Gamma Correction
        R_gamma_L = 3.9964;
        G_gamma_L = 3.6888;
        B_gamma_L = 6.0477;
        %K_gamma_L = 3.6886;
    
        % The right eye values:
        R_gamma_R = 3.9964;
        G_gamma_R = 3.6888;
        B_gamma_R = 6.0477;
        %K_gamma_R = 3.6886;
    
        % Do the correction. We can enter either 1 value (across all guns) or 3 (1 for each gun)
        % Note that the gammas are inverted here.
        PsychColorCorrection('SetEncodingGamma', dpy.win, 1./[mean([R_gamma_L, R_gamma_R]), mean([G_gamma_L, G_gamma_R]), mean([B_gamma_L, B_gamma_R])]);
    
        %Datapixx('EnableVideoLcd3D60Hz');
        Datapixx('DisableVideoLcd3D60Hz'); %> seems to work better, according to Daniel.
        
        Datapixx('RegWr'); % Make sure the commands above are executed
        
    catch
        lasterr
        error('Error using Datapixx - are you sure you are connected to one?');
    end
end

%we have to edit this based on display
dpy.res = [1920 1080]; % screen resoloution
dpy.size = [.522 .292]; % Meters
dpy.distance = [.57]; % Metres back from screen
dpy.frameRate=120;
dpy.center = dpy.res/2;
% dpy will eventually contain all the info about the display e.g. size,
% distance, refresh rate, spectra, gamma.
% For now if just has the gamma function (inverse) in it.
qqq
%% Here we set out experimental size paramaters
%FOV = [0 0] and 2
%ANU = [0 3] 5 (3-5 degree annulus)
%FULL = [2 3] and 5 (1 degree grey annulus)

mask.radiusDeg=[2 3]; % Inner and outer radii of the mask (this can be zero if you like)
radiusDeg=5; % Total stimulus radius

%here we set our other paramaters
tfList=[1,2,4;1,2,4]'; % This is in Hz. Note integer numbers of frames on a 120Hz mon
sfList=[1]'; %cpd
contList=[3,6,12,24,48;0 0 0 0 0]'/100;  % Contrasts to probe
nTF=size(tfList,1);% number of TFs
nSF=size(sfList,1);
nCont=size(contList,1);%number of Contrasts
ordered=1:(nTF*nCont);

stim.spatial.internalRotation = 0; % Does the grating rotate within the envelope?
stim.rotateMode = []; % rotation of mask grating (1= horizontal, 2= vertical, etc?)
stim.spatial.angle = [0,90]; % possible angle of gratings on screen
stim.temporal.duration=11; % how long to flicker for - remember we remove first second
stim.temporal.modtype = 'onoff'; %this is our stim flicker type - can be on/off or square

stim.mask=mask;
stim.radiusDeg=radiusDeg;

% In principle we can have both mask and probe here. For now though we set
% the mask to zero. probeCont isn't used in this code (yet)
probeCont=[99]/100;
maskCont =[];

nConds=length(ordered); %number of total conditions
shuffleSeq=Shuffle(ordered);% This fully shuffles the order, we have the same order in all 4 runs.

nRepeats=4; %how many times to repeat each condition

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

else
    disp('Dummy run - not initializing the screen');
    
    dpy.win=-1;
end
tic

% %Present Intro Screen
% imDir = '/Users/tyrion/Documents/MATLAB/Marc/EEG_ANT_VPIXX/V2/WadeLabScreen.jpg';
% im = imread(imDir);
% IntroTexture = Screen('MakeTexture', dpy.win, im);
% Screen('DrawTexture', dpy.win, IntroTexture, [], [], 0)
% Screen('Flip', dpy.win);
% KbWait;

% %Present instructions
% pause(1);
% Screen('TextSize', dpy.win, 30);
% DrawFormattedText(dpy.win, 'Instructions:\n\n In this experiment you are required to fixate on the black central fixation cross whilst a pattern is flashed for 11 seconds. \n\n Please try not to blink during this time. You can blink when the cross changes to white for 3 seconds between trials.\n \n There will be rest points during the experiment. \n\n\n Do NOT press spacebar yet. \n\n The experimenter will tell you when to press spacebar to begin. ','center','center');
% Screen('Flip', dpy.win);
% KbWait;
% %change this after piloting to 5s
% pause(1)

for thisRun=1:nRepeats  % 5 mins per rep

    % Write something to the DataPixx to say that we are about to start the
    % sequence
    % Write out a trigger to the DataPixx if required
    if DPIXXPRESENT
        Datapixx('SetDoutValues', transformindex(thisRun));
        Datapixx('RegWrRd'); %This flushes the register
        
    else
        fprintf('\nCurrent condition trigger code %d\n',thisRun);
    end
    
    pause(.01);
      
    for thisCond=1:nConds
        
        stim.spatial.angle = shuffle(stim.spatial.angle); %shuffle our orientations

        stim.spatial.frequency=[sfList(1), sfList(1)]; % Set the same spatial freq for both components
        % Phase is the phase shift in degrees (0-360 etc.)applied to the sine grating:
        stim.spatial.phase=[0 0 ]; %[rand(1)*360 rand(1)*360];
        stim.spatial.pOffset=rand(2,1)*360;
        
        thisaction= shuffleSeq(thisCond);
        tIndex=ceil(thisaction/ nCont); % Index into the list of temporal frequencies
        cIndex=1+rem((thisaction-1), nCont); %  Should be 1+rem(thisAction-1,nTF)  - index into the list of contrasts
        
        thisTemp(thisRun,thisCond)=tIndex;
        thisCont(thisRun,thisCond)=cIndex;
        
        stim.temporal.frequency=tfList(tIndex,:);
        
        % Compute increment of phase shift per redraw:
        dpy.phaseincrement = [stim.temporal.frequency] * 360 * dpy.ifi;
        
        stim.cont=contList(cIndex,:);
        
        
        disp(thisRun)
        disp(thisCond)
        
        fprintf('\nRunning %.2d %.2d',stim.cont(1),stim.temporal.frequency(1));
        
        % Write out a trigger to the DataPixx if required
        if DPIXXPRESENT
            Datapixx('SetDoutValues', transformindex(thisaction+100));
            Datapixx('RegWrRd');
            pause(.1);
        else
            fprintf('\nCurrent condition trigger code %d\n',thisaction+100);
        end
        
        % Call the main display loop - the screen is already open so this
        % should be fast
        if (~DUMMYRUN) % Do we display something or not?
            d=EEG_runPlaid_V5(dpy,stim,dataPixxOk);
            [keyIsDown, secs, keyCode, deltaSecs] = KbCheck;
            
            if keyIsDown % This is a vector of key codes
                keyDown=find(keyCode); % Which things are down? We only take the first one..
                if (KbName(keyDown(1))=='q')
                    sca
                    error('abort experiment')
                end% End check on 'q' key
            end % End check on key down
                
        end % End check on dummy run
        
        if DPIXXPRESENT
            % Send 255 to indicate end of presentation
            Datapixx('SetDoutValues', transformindex(255));
            Datapixx('RegWrRd');
            
            Screen('TextSize', dpy.win, 25);
            DrawFormattedText(dpy.win,'Blink\n\n\n + \n\n\nBlink','center','center', [255, 255, 255]);
            Screen('Flip', dpy.win);
            WaitSecs(3) %3s ISI
                     
        else
            fprintf('\nCurrent condition trigger code %d\n',255);
        end
       
    end % Next contrast and temporal frequency pair
     
    %end of first set of reps, press keyboard to cont, if 4, present end
    %msg
 
        if thisRun==nRepeats;    
            Screen('TextSize', dpy.win, 45);
            DrawFormattedText(dpy.win, 'The experiment is finished. \n\n Please press spacebar to exit.', 'center', 'center');
            Screen('Flip', dpy.win);
            KbWait;
        else
            Screen('TextSize', dpy.win, 45);
            DrawFormattedText(dpy.win, 'Take a break. Please press spacebar when you are ready to continue','center','center');
            Screen('Flip', dpy.win);
            KbWait;
            pause(1);
        end 
end % Next repetition

if DPIXXPRESENT
    % Send 255 to indicate end of presentation
    Datapixx('SetDoutValues', transformindex(254)); % End of all reps
    Datapixx('RegWrRd');
    
else
    fprintf('\nCurrent condition trigger code %d\n',254);
end

% We're done. Close the window. This will also release all other ressources:
% Bye
totalSessionTime=toc;

if (~DUMMYRUN) % Save as soon as possible in case it crashes
    save(out_file)
end

if DPIXXPRESENT       % Close the DataPixx
    if(Datapixx('IsReady'))
        Datapixx('DisableVideoScanningBacklight');      % optionally, turn it off first, in case the refresh rate has changed since startup
        Datapixx('Close');
    end 
end

Screen('CloseAll');
sca
  
   