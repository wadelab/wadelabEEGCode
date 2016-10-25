function [dataOut]=EEG_runPlaid_V5(dpy,stim,dpixxStatus)
% function dataOut=EEG_runplaid_donut(dpy,stim,dpixxStatus)
% Generates a 2-component plaid
% Grating 1 is orthogonal to grating 2
% All the inputs are now 2 element vectors
% that specify the parameters for ggating 1 and grating 2
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
% MMH May 2016 edited to remove DAQ and include EEG triggers
% MMH 11/11/2016 edited to have stimuli in an annulus. Spmething is wrong
% because the fixation cross is not centered, so for now we manually adjust
% the screen res so it is centered... but will have to fix it because thats
% bad science.

% Make sure this is running on OpenGL Psychtoolbox:
%AssertOpenGL;

%% Build a procedural sine grating texture for a grating with a support of
% res(1) x res(2) pixels and a RGB color offset of 0.5 -- a 50% gray.

% In this version we asume that stim has
% stim.radiusDeg  (scalar size of the checkerboard stim in degrees )
% stim.mask.radii (2x1 vector with inner and outer radius of the mask)

disp('Running');
tic

Screen('BlendFunction', dpy.win, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

disp(stim.cont)
toc

% Compute the alpha and amplitudes that we will use
[amps,alpha]=flytv_computeAlphaAmps(stim.cont);


% verticalblank = flip
vbl = Screen('Flip', dpy.win);


% We run at most 'movieDurationSecs' seconds if user doesn't abort via keypress.
vblendtime = vbl + stim.temporal.duration;
i=0;

% Update some grating animation parameters:
% Also used for Annulus set up
phase=stim.spatial.phase;
degToRad=pi/180;

pixelsPerMeter=dpy.res(1)/dpy.size(1);
metersPerDegree=dpy.distance*tan(degToRad);

pixPerDegree=pixelsPerMeter*metersPerDegree;

stim.spatial.frequencyCPerPixel=stim.spatial.frequency/pixPerDegree;
%generate our sine wave gratings (bottom is mask- we don't use this right now but we might later)
gratingtex1 = CreateProceduralSineGrating(dpy.win, dpy.res(1), dpy.res(2),[.5,.5,.5, 1], stim.radiusDeg*pixPerDegree); % Bottom grating. Note by specifying a radius we force it to be circular
gratingtex2 = CreateProceduralSineGrating(dpy.win, dpy.res(1), dpy.res(2),[.5 .5 .5 alpha], stim.radiusDeg*pixPerDegree); % Top grating blend 50%

%% Annulus set up

% Compute the inner and outer radii in pixels instead of degrees
% This comes in as part of 'stim'
% If no mask is present we just ignore this bit.
if (~isfield(stim,'mask'))
    stim.mask.radiusDeg=[1 2];
end

xScreen = dpy.res(1);% width resolution
yScreen = dpy.res(2);
halfSizeX=round(xScreen/2);
halfSizeY=round(yScreen/2);
%create meshgrid with coordinates of screen width
[xg,yg]=meshgrid(-(halfSizeX-1):halfSizeX, -(halfSizeY-1):halfSizeY);

% Trick here is to compute two masks: one from the center to the inner
% diameter, another from the inner to the outer. Then we multiply these
% together to get just the intersection.

radius_inner_pixels=stim.mask.radiusDeg(1)*pixPerDegree;
radius_outer_pixels=stim.mask.radiusDeg(2)*pixPerDegree;

transLayer=2;%for mask. I think this has to change if we're doing color but for now it's okay.

%create mask
finalMask=ones(yScreen, xScreen,2) * 128; % Make a gray thing size of the actual screen with two layers: image, alpha

dist=sqrt(xg.^2+yg.^2);%  dist from center according to pythagoras
blur=100; % lower = more blur
dm=0; % center
annulusMaskInner(:,:)=round(exp(-((dist-dm)/radius_inner_pixels).^blur)); %mask transparency 0 = 0%, 255 = 100%.
annulusMaskOuter(:,:)=round(exp(-((dist-dm)/radius_outer_pixels).^blur)); %mask transparency 0 = 0%, 255 = 100%.
finalMask(:,:,transLayer)=((1-annulusMaskInner(:,:)).*annulusMaskOuter(:,:))*255;

% Here we ask whether the modulation is square wave or on/off.
if ~isfield (stim.temporal,'modtype')
    stim.temporal.modtype = 'onoff';
end

%for square wave both sides of the duty cycle have the same magnitude (1 1)
%for on/off one side has a magnitude of 0

if strcmp(stim.temporal.modtype,'square')
    modMag = [1 1]; %square
else
    modMag = [1 0]; %on off
end

% Build a the texture:
maskTexture=Screen('MakeTexture', dpy.win, finalMask);
% Texture placement for centre.. this may be off


if (dpixxStatus)
       lTime=vbl;
    while (vbl < vblendtime)%while we are less than the stim duration
        % Increment phase by the appropriate amount for this time period:
        phase = phase + dpy.phaseincrement;
        pMod = 180*(round(phase/180)); %modulate phase
        
        thisMag = modMag(round((mod(phase,180)/180))+1);
        
        % Draw the grating, centered on the screen, with given rotation 'angle',
        % sine grating 'phase' shift and amplitude, rotating via set
        % 'rotateMode'. Note that we pad the last argument with a 4th
        % component, which is 0. This is required, as this argument must be a
        % vector with a number of components that is an integral multiple of 4,
        % i.e. in our case it must have 4 components:
        
        Screen('DrawTexture', dpy.win, [gratingtex1], [], [], [stim.spatial.angle(1)], [], [0], [], [], [stim.rotateMode], [pMod(1),stim.spatial.frequencyCPerPixel(1),amps(1)*thisMag(1),0]');
        Screen('DrawTexture', dpy.win, [gratingtex2], [], [], [stim.spatial.angle(2)], [], [0], [], [], [stim.rotateMode], [pMod(2),stim.spatial.frequencyCPerPixel(2),amps(2)*thisMag(2),0]');
        Screen('DrawTexture', dpy.win, [maskTexture],[]); % mask
        
        Screen('TextSize', dpy.win, 25);
        DrawFormattedText(dpy.win,'+','center','center', [0, 0, 0]);
        
        % Show it at next retrace:
        vbl = Screen('Flip', dpy.win, vbl + 0.5 * dpy.ifi);
        
        %       Sent a pulse on the DPIxx (if present) once per second
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
        
        %Screen('DrawTexture', winptr, txtptr [,sourceRect] [,destinationRect] [,rotationAngle] [, filterMode] [, globalAlpha] [, modulateColor] [, textureShader] [, specialFlags] [, auxParameters]);
        Screen('DrawTexture', dpy.win, [gratingtex1], [], [], [stim.spatial.angle(1)], [], [0], [], [], [stim.rotateMode], [pMod(1),stim.spatial.frequencyCPerPixel(1),amps(1)*thisMag(1),0]');
        Screen('DrawTexture', dpy.win, [gratingtex2], [], [], [stim.spatial.angle(2)], [], [0], [], [], [stim.rotateMode], [pMod(2),stim.spatial.frequencyCPerPixel(2),amps(2)*thisMag(2),0]');
        Screen('TextSize', dpy.win, 25);
        DrawFormattedText(dpy.win,'+','center','center',[0, 0, 0]);
        
        % Show it at next retrace:
        vbl = Screen('Flip', dpy.win, vbl + 0.5 * dpy.ifi);
    end
end

% Leave the screen in place for the next run
dataOut=0; % For the flies we actually collect data per run. Here this structure will probably be something to do with conditioins and triggers.
