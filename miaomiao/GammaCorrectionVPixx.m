
% Code necessary for gamma correction of the B114 Viewpixx, using measurements obtained 
% without the 3D goggles on 22 Feb 2016
% R Maloney, 22/7/16


%----------------------------
%     Set up the screen
%----------------------------

% initialization of the display
AssertOpenGL;
% Open PTB onscreen window: We request a 32 bit per colour component
% floating point framebuffer if it supports alpha-blending. Otherwise
% the system shall fall back to a 16 bit per colour component framebuffer:
PsychImaging('PrepareConfiguration');
PsychImaging('AddTask', 'General', 'FloatingPoint32BitIfPossible');

% required for gamma correction through the PsychImaging pipeline:
PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'SimpleGamma');

% Set the color range to be normalised between 0 and 1 (rather than 0-255).
% *** Note that this means your pixel intensity values that define your stimuli should all be between 0-1
PsychImaging('AddTask', 'General', 'NormalizedHighresColorRange', 1);

% Open an on screen (grey) window and configure the imaging pipeline
[win, windowRect] = PsychImaging('OpenWindow', WhichScreen, 0.5);


PsychImaging('PrepareConfiguration');

% Initialise the Vpixx device:

% *** NOTE: If just using the Viewpixx as a normal screen (ie nothing fancy or 3D)
% it may not be necessary to do the following (though it probably doesn't hurt...)
% Many of the calls have been left in, just commented out

PsychImaging('AddTask', 'General', 'UseDataPixx');

Datapixx('Open');
% The following commands are included in demos that apparently work for both the Viewpixx AND the PROpixx, though they seem specific to the Viewpixx...
Datapixx('DisableVideoScanningBacklight');    % optionally, turn it off first, in case the refresh rate has changed since startup
Datapixx('EnableVideoScanningBacklight');     % Only required if a VIEWPixx.
Datapixx('EnableVideoStereoBlueline');
%Datapixx('SetVideoStereoVesaWaveform', 2);    % If driving NVIDIA glasses

% Prepare for the M16 mode for 16 bit luminance output precision: important when varying contrast
% PsychImaging('AddTask', 'General', 'EnableDataPixxM16OutputWithOverlay'); >couldn't get working:RM

% Now also prepare the overlay. This is to allow drawing of the blue lines (I think)?
%overlay = PsychImaging('GetOverlayWindow', win);

if Datapixx('IsViewpixx3D') % If it's the Viewpixx3D
    
    % Do the gamma correction. When using the Viewpixx/PROpixx through the PsychImaging pipeline, we
    % SHOULD NOT use Screen(?LoadNormalizedGamma?) (see http://www.jennyreadresearch.com/research/lab-set-up/datapixx/)
    % The PROpixx device should have a linear LUT built in, but we will add this here for completeness.
    
    % The gamma exponents below were obtained from measurements made with the Jaz Spectrometer on 22 Feb 2016
    % at 57 cm viewing distance but with no 3D goggles present.
    % Here we'll enter 3 values, 1 for each of the R, G & B guns.
    % But we could also just enter a single value, for the K channel (all guns combined). They are also provided below.
    
    % We will average the left and right eye measurements, which *should* be the same
    
    % The left eye values:
    R_gamma_L = 3.2512;
    G_gamma_L = 3.2052;
    B_gamma_L = 4.5724;
    %K_gamma_L = 3.2052;
    
    % The right eye values:
    R_gamma_R = 3.5668;
    G_gamma_R = 3.4843;
    B_gamma_R = 5.6056;
    %K_gamma_R = 3.5217;
    
    % Do the correction. We can enter either 1 value (across all guns) or 3 (1 for each gun)
    % Note that the gammas are inverted here.
    PsychColorCorrection('SetEncodingGamma', win, 1./[mean([R_gamma_L, R_gamma_R]), mean([G_gamma_L, G_gamma_R]), mean([B_gamma_L, B_gamma_R])]);
    
    %Datapixx('EnableVideoLcd3D60Hz');
    Datapixx('DisableVideoLcd3D60Hz'); %> seems to work better, according to Daniel.
    Datapixx('RegWr');
    
end

