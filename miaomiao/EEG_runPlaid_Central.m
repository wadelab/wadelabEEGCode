function [dataOut]=EEG_runPlaid_Central(dpy,stim,dpixxStatus)
% function dataOut=flytv_PlaidDemo2(cyclespersecond, sfreq,contrast)
% Generates a 2-component plaid
% Grating 1 is orthogonal to grating 2
% All the inputs are now 2 element vectors
% that specify the parameters for grating 1 and grating 2
% So for example
% cyclespersecond=[0.1 0.2]; % The second grating has twice the sf of the
% first one
% Contrast is the contrast of each component. [.5 .5] gives equal
% mixtures of gr 1 and gr 2 with a max screen level of 100%
% You can go higher : so 0.7 .3 is okay
% The conversion from contrast to amplitudes in the code is computed by a
% separate function flytv_computeAlphaAmps
%

% History:
% 3/1/9  mk   Written.
% 25/04/14 rw531 edited for flytv
% MMH May 2015 edited to remove DAQ and include EEG trigger
% Make sure this is running on OpenGL Psychtoolbox:
%AssertOpenGL;

% Build a procedural sine grating texture for a grating with a support of
% res(1) x res(2) pixels and a RGB color offset of 0.5 -- a 50% gray.

% Begin data acquisition in the background
disp('Running');
tic


Screen('BlendFunction', dpy.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

disp(stim.cont)
toc

% Compute the alpha and amplitudes that we will use
[amps,alpha]=flytv_computeAlphaAmps(stim.cont);

gratingtex1 = CreateProceduralSineGrating(dpy.win, dpy.res(1), dpy.res(2),[.5,.5,.5, 1], 73.191); % Bottom grating
gratingtex2 = CreateProceduralSineGrating(dpy.win, dpy.res(1), dpy.res(2),[.5 .5 .5 alpha], 73.191); % Top grating blend 50%

% Wait for release of all keys on keyboard, then sync us to retrace:

vbl = Screen('Flip', dpy.win);

% We run at most 'movieDurationSecs' seconds if user doesn't abort via keypress.
vblendtime = vbl + stim.temporal.duration;
i=0;
% Update some grating animation parameters:
phase=stim.spatial.phase;
degToRad=pi/180;

pixelsPerMeter=dpy.res(1)/dpy.size(1);
metersPerDegree=dpy.distance*tan(degToRad);

pixPerDegree=pixelsPerMeter*metersPerDegree;

stim.spatial.frequencyCPerPixel=stim.spatial.frequency/pixPerDegree;
if (dpixxStatus)
%     Datapixx('StopAllSchedules');
%     Datapixx('RegWrRd');    % Synchronize DATAPixx registers to local register cache
%     
%     % We'll make sure that all the TTL digital outputs are low before we start
%     Datapixx('SetDoutValues', 0);
%     Datapixx('RegWrRd');
%     doutWave=[ones(60,1)*2;zeros(60,1)];
%     bufferAddress = 8e6;
%     Datapixx('WriteDoutBuffer', doutWave, bufferAddress);
%     
%     % Define the schedule which will play the wave.
%     samplesPerTrigger = length(doutWave)
%     Datapixx('SetDoutSchedule', 0, [120 1], 120, bufferAddress, samplesPerTrigger);
%     
%     % Insert your visual stimulus setup code here, finishing up with a Screen('Flip', window);
%     
    % Tell the trigger schedule to start on the next vertical sync pulse
%     Datapixx('StartDoutSchedule');
%     Datapixx('RegWrVideoSync');
    lTime=vbl;
    while (vbl < vblendtime)
        % Increment phase by the appropriate amount for this time period:
        phase = phase + dpy.phaseincrement;
        pMod = 180*(round(phase/180 ));
        
        
        % Draw the grating, centered on the screen, with given rotation 'angle',
        % sine grating 'phase' shift and amplitude, rotating via set
        % 'rotateMode'. Note that we pad the last argument with a 4th
        % component, which is 0. This is required, as this argument must be a
        % vector with a number of components that is an integral multiple of 4,
        % i.e. in our case it must have 4 components:
     
        
        Screen('DrawTexture', dpy.win, [gratingtex1], [], [], [stim.spatial.angle(1)], [], [0], [], [], [stim.rotateMode], [pMod(1),stim.spatial.frequencyCPerPixel(1),amps(1),0]');
        Screen('DrawTexture', dpy.win, [gratingtex2], [], [], [stim.spatial.angle(2)], [], [0], [], [], [stim.rotateMode], [pMod(2),stim.spatial.frequencyCPerPixel(2),amps(2),0]');
        Screen('TextSize', dpy.win, 35);
        DrawFormattedText(dpy.win,'+','center','center', [0, 0, 0]);
        
        % Show it at next retrace:
        vbl = Screen('Flip', dpy.win, vbl + 0.5 * dpy.ifi);
%        Sent a pulse on the DPIxx (if present) once per second
       
        if (vbl-lTime)>=1
            for t=1:5
                Datapixx('SetDoutValues', transformindex(2));
                Datapixx('RegWrRd');
            end
             for t=1:5
                Datapixx('SetDoutValues', transformindex(0));
                Datapixx('RegWrRd');
            end
            lTime=vbl;
        end
        
            
        
    end
    Datapixx('StopDoutSchedule');
    Datapixx('RegWrRd');
else
    while (vbl < vblendtime)
        % Increment phase by the appropriate amount for this time period:
        phase = phase + dpy.phaseincrement;
        pMod = 180*(round(phase/180 ));
        
        
        % Draw the grating, centered on the screen, with given rotation 'angle',
        % sine grating 'phase' shift and amplitude, rotating via set
        % 'rotateMode'. Note that we pad the last argument with a 4th
        % component, which is 0. This is required, as this argument must be a
        % vector with a number of components that is an integral multiple of 4,
        % i.e. in our case it must have 4 components:
        
        Screen('DrawTexture', dpy.win, [gratingtex1], [], [], [stim.spatial.angle(1)], [], [0], [], [], [stim.rotateMode], [pMod(1),stim.spatial.frequencyCPerPixel(1),amps(1),0]');
        Screen('DrawTexture', dpy.win, [gratingtex2], [], [], [stim.spatial.angle(2)], [], [0], [], [], [stim.rotateMode], [pMod(2),stim.spatial.frequencyCPerPixel(2),amps(2),0]');
        Screen('TextSize', dpy.win, 35);
        DrawFormattedText(dpy.win,'+','center','center',[0, 0, 0]);
        
        % Show it at next retrace:
        vbl = Screen('Flip', dpy.win, vbl + 0.5 * dpy.ifi);
        % Sent a pulse on the DPIxx (if present) once per second
    end
end


% Leave the screen in place for the next run
dataOut=0; % For the flies we actually collect data per run. Here this structure will probably be something to do with conditioins and triggers.
