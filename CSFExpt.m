classdef CSFExpt
    %% Public Properties
    properties (Access = public)
        % Configuration flags
        DEBUG = 0; DEBUG_WTS = 0; DEBUG_STIMCHANGE = 0; DEBUG_PLOTFFT = 0; DEBUG_RECORDMOVIES = 0;
        preload = 1; waitframes = 1;
        DARK = 1; jitter_flag = 0;

        % Stimulus variable(s)
        % centre of the stimulus (in degrees of VA - usually set this to 0)
        % width/height of the stimulus prior to gaussian envelope
        % gaussian envelope over the grating (in degree of visual angle)
        stim = struct('condition', '', 'dur', 800, 'center', [0,0], 'size', [8,8], 'sigma', 1, 'numTrialsPerFreq', 0);
        input = {};

        % Display variables: Room 510 - Display++
        % display.width: measure the physical width of your screen in cm;
        %display.fRate: Frame rate of the display - 120Hz for ours
        % dispay.dist: in cm; measure your viewing distance
        display = struct('resolution', [1920 1080], ...
            'width', 70, 'dist', 206, ...
            'screen', 0, ...
            'fRate', 120, ...
            'pxPerDeg', 0, ...
            'anglePerPix', 0, ...
            'gray', 0);
        %% opto variables
        p = {};
        %% Keyboard variables
        homedir = cd;
        escapeKey;upKey;downKey;leftKey;rightKey;
        responseKey;keybIndices;keybNames;

        %% trial variables
        numTrials = [];
        FF = [];
        fixFreq; % which frequency do we need to fix spatial/temporal - based on the selected expt condition
        t;
        sID;
        time;
        type; % to check if 'baseline'/ 'opto'/ 'eye-movements'
        eyeData;

        %% sound variables
        beepBegin; beepHigh; beepLow;

        %% screen variables
        theScreen;
        flipHandle;
    end


    methods (Access = public)
        function this = CSFExpt(cond, sID, config)
            this.stim.condition = cond;
            this.sID = sID;
            this.DEBUG = config.DEBUG;
            this.DEBUG_STIMCHANGE = config.DEBUG_STIMCHANGE;
            this.DEBUG_RECORDMOVIES = config.DEBUG_RECORDMOVIES;
            this.DEBUG_PLOTFFT = config.DEBUG_PLOTFFT;
            if this.DEBUG_STIMCHANGE
                this.DEBUG_WTS = 1;
            else
                this.DEBUG_WTS = config.DEBUG_WTS;
            end
            this.jitter_flag = config.jitter_flag;
            if isfield(config, 'fRate') && (config.fRate ~= 0)
                this.display.fRate = config.fRate;
            end
            retCode = SetDependentLibraries(this);
            if retCode ~= 1
                sprintf('Error: The dependent libraries were not set correctly')
                return;
            end

            this.type = split(this.stim.condition, "_");
            if (strcmp(this.type{1,1},'opto') || strcmp(this.type{1,1}, 'eye')) && strcmp(config.photoswitch,'')
                sprintf('Please specify a valid photoswitch for the optogenetic protein. Bye!')
                return
            end

            this.time.startExp = datestr(now);
            this = SetExptVariables(this, config.photoswitch);
            %             if ~this.DEBUG_WTS
            this = SetScreenSettings(this);
            %             end
            setKbSettings(this);
            this = SetFbSoundSettings(this, 0.02, 0.04, 0.01);
            this = DefineInputParams(this, config.inputType);
            if strcmp(this.type{1,1},'eye')
                this = ImportEyeMovements(this);
            end
        end

        %% Routine to run the experiment
        function this = RunExpt(this)
            try

                ListenChar(2);
                ApplyKbFilter;

                % create a data structure to store all levels
                data = dictionary();
                redoTrials = []; % list of trials repeated

                sub_id = randi(29);
                for i= 1: length(this.FF)
                    qcsf = SetQCSF(this);

                    % initialize redoTrial flag to avoid skipping frames
                    redoTrial = 0;
                    breakExperiment = 0;
                    for trial = 1:this.stim.numTrialsPerFreq
                        tic
                        %% Go through trial
                        % Play trial start sound
                        PlayFbSounds(this, '0');
                        if ~this.DEBUG_WTS
                            Screen('FillRect',this.theScreen,this.display.gray)
                            DrawFormattedText(this.theScreen, '.', 'center', 'center', 0);
                            Screen('Flip',this.theScreen)
                        end
                        %%
                        if strcmp(this.type{1,1},'opto') || strcmp(this.type{1,1}, 'eye')
                            % WaitSecs(0.005);
                        else
                            WaitSecs(0.5);
                        end
                        qcsf.data.trial = trial;
                        [qcsf,VF,contrast] = runQCSF(qcsf, 'pretrial');
                        %% Determine whether contrast stim is -45 or 45 deg (or 0/90, or whatever you end up picking)
                        CR = ceil(2*rand);
                        %                         CR = randi(2);
                        if CR == 1
                            orient = -45; % in degress
                        elseif CR == 2
                            orient = 45;
                        end

                        %                         if this.DEBUG_STIMCHANGE
                        %                             orient = 0;
                        %                         end

                        if (strcmp(this.fixFreq,'temporal'))
                            %fprintf('SF = %f, TF = %f, contrast = %f, orient=%d \n',VF, this.FF(i), contrast, orient);

                            if this.DEBUG_STIMCHANGE
                                clear mex;
                                VF = 2;
                                this.FF(i) = 0;
                                contrast = 1;
                            end
                            if strcmp(this.type{1,1},'baseline')
                                grating = GenerateGrating(this, VF, this.FF(i), contrast, orient);
                                stimulus = ScaleStimulus(this, grating, this.display.gray, 1, contrast, this.FF(i));

                            elseif strcmp(this.type{1,1},'opto')
                                grating = GenerateGrating(this, VF, this.FF(i), contrast, orient);
                                opto = ApplyOptoFilter(this, grating, this.FF(i), contrast);
                                stimulus =  ScaleStimulus(this, opto, this.p.delta, this.p.offset, contrast, this.FF(i));

                            elseif strcmp(this.type{1,1},'eye')
                                break_cond = 0;
                                while (break_cond == 0)
                                    trial_eyeData = this.eyeData{randi(size(this.eyeData, 1)), sub_id};
                                    if ~isempty(trial_eyeData)
                                        break_cond = 1;
                                    end
                                end

                                if this.DEBUG_STIMCHANGE
                                    %                                     trial_eyeData = trial_eyeData(:,1:2);
                                    %                                     trial_eyeData(4:5,2) = 60;
                                    %                                     trial_eyeData(3,1) = 45;
                                    grating_ctl = GenerateGrating(this, VF, this.FF(i), contrast, orient);
                                    stim_ctl =  ScaleStimulus(this, grating_ctl, this.p.delta, this.p.offset, contrast, VF);
%                                     sprintf('Range of baseline: = %f', (max(max(stim_ctl)) - min(min(stim_ctl))))

                                    opto_ctl = ApplyOptoFilter(this, grating_ctl, VF, contrast);
                                    stim_opto =  ScaleStimulus(this, opto_ctl, this.p.delta, this.p.offset, contrast, VF);
%                                     sprintf('Range of opto: = %f', (max(max(stim_opto)) - min(min(stim_opto))))


                                    grating = GenerateEMGrating(this, VF, this.FF(i), contrast, orient, trial_eyeData); % Generates grating with eye-movement if 'eye'
                                    opto = ApplyOptoFilter(this, grating, VF, contrast);
                                    eye = CorrectEyeMovement(this, opto, trial_eyeData, VF, orient);
                                    stim_eye =  ScaleStimulus(this, eye, this.p.delta, this.p.offset, contrast, VF);
%                                     sprintf('Range of opto+eye: = %f', (max(max(stim_eye)) - min(min(stim_eye))))
                                end

                                grating = GenerateEMGrating(this, VF, this.FF(i), contrast, orient, trial_eyeData); % Generates grating with eye-movement if 'eye'
                                opto = ApplyOptoFilter(this, grating, this.FF(i), contrast);
                                eye = CorrectEyeMovement(this, grating, trial_eyeData, VF, orient);
                                stimulus =  ScaleStimulus(this, eye, this.p.delta, this.p.offset, contrast, this.FF(i));

                            end


                        elseif(strcmp(this.fixFreq,'spatial'))

                            if this.DEBUG_STIMCHANGE
                                clear mex;
                                %                                 this.FF(i) = 2;
                                %                                 VF = 5;
                                contrast = 1;
                            end
                            %fprintf('SF = %f, TF = %f, contrast = %f, orient=%d \n', this.FF(i), VF, contrast, orient);
                            if strcmp(this.type{1,1},'baseline')
                                grating = GenerateGrating(this, this.FF(i), VF, contrast,orient);
                                stimulus =  ScaleStimulus(this, grating, this.display.gray, 1, contrast, VF);

                            elseif strcmp(this.type{1,1},'opto')
                                grating = GenerateGrating(this, this.FF(i), VF, contrast,orient);
                                opto = ApplyOptoFilter(this, grating, VF, contrast); % returns the texture if opto
                                stimulus =  ScaleStimulus(this, opto, this.p.delta, this.p.offset, contrast, VF);

                            elseif strcmp(this.type{1,1},'eye')
                                break_cond = 0;
                                while (break_cond == 0)
                                    trial_eyeData = this.eyeData{randi(size(this.eyeData, 1)), sub_id};
                                    if ~isempty(trial_eyeData)
                                        break_cond = 1;
                                    end
                                end

                                if this.DEBUG_STIMCHANGE
                                    grating_opto = GenerateGrating(this, this.FF(i), VF, contrast, orient);
                                    opto_only = ApplyOptoFilter(this, grating_opto, VF, contrast);
                                    stimulus_opto =  ScaleStimulus(this, opto_only, this.p.delta, this.p.offset, contrast, VF);
                                    grating_eye = GenerateEMGrating(this, this.FF(i), VF, contrast, orient, trial_eyeData); % Generates grating with eye-movement if 'eye'
                                    opto = ApplyOptoFilter(this, grating_eye, VF, contrast);
                                    eye = CorrectEyeMovement(this, opto, trial_eyeData, this.FF(i), orient);
                                    stimulus_eye =  ScaleStimulus(this, eye, this.p.delta, this.p.offset, contrast, VF);
                                else
                                    grating = GenerateEMGrating(this, this.FF(i), VF, contrast, orient, trial_eyeData); % Generates grating with eye-movement if 'eye'
                                    opto = ApplyOptoFilter(this, grating, VF, contrast);
                                    eye = CorrectEyeMovement(this, opto, trial_eyeData, this.FF(i), orient);
                                    stimulus =  ScaleStimulus(this, eye, this.p.delta, this.p.offset, contrast, VF);
                                end
                            end

                        end
                        if this.DEBUG_STIMCHANGE
                            Tex = MakeTex(this, stim_ctl, stim_opto, stim_eye);
                        else
                            Tex = MakeTex(this, stimulus);
                        end
                        if ~this.DEBUG_WTS
                            if this.preload
                                Screen('PreloadTextures', this.theScreen, Tex);
                                % disp('Textures are preloaded');
                            end
                        end


                        % Present qCSF
                        StartT = GetSecs; %Measure start time of session
                        k = 0;
                        frame_time = zeros(length(this.t),1);
                        while (GetSecs() - StartT) < this.stim.dur/1000
                            k = ceil((GetSecs() - StartT) .* this.display.fRate); % current frame based on current time
                            if k <= (this.stim.numFrames)
                                if ~this.DEBUG_WTS
                                    Screen('DrawTextures', this.theScreen, Tex(k));
                                    Screen('Flip', this.theScreen); % show to participant
                                end
                                frame_time(k) = 1;
                            end
                        end
                        %% check if there are multiple consecutive broken frames
                        % (3 if max TF is 30)
                        % flag if it is so, and redo the trial, don't save
                        % the broken data
                        if any(movsum((diff(find(frame_time==0))==1),3) >=3)
                            % redoTrial = 1;
                            redoTrials = [redoTrials, trial];
                        else
                            redoTrial =0;
                        end

                        % Finish with background
                        if ~this.DEBUG_WTS
                            Screen('FillRect', this.theScreen, this.display.gray);
                            Screen('Flip', this.theScreen);
                        end
                        toc

                        tic
                        [keyPressed] = wait4key({'LeftArrow','RightArrow', 'ESCAPE', 'q'}, -3);
                        sprintf(keyPressed);
                        [returnCode, response] = PlayFbSounds(this,keyPressed, CR);
                        if returnCode == 0 && breakExperiment == 1
                            fprintf("Experiment quit by the user; Pressed %s key on trial = %d \n", keyPressed, trial);
                            if (strcmp(this.fixFreq,'temporal'))
                                fprintf('SF = %f, TF = %f',VF, this.FF(i));
                            else
                                fprintf('SF = %f, TF = %f',this.FF(i), VF);
                            end
                            breakExperiment = 1;
                            Screen('CloseAll');
                            break;
                        elseif returnCode == -1 && breakExperiment == 1
                            fprintf('Experiment broken by the user; Pressed an invalid key on trial = %d \n', trial);
                            breakExperiment = 1;
                            continue;
                        end

                        toc
                        WaitSecs(0.3);

                        % record and update if it is not a repeated trial
                        if ~redoTrial
                            qcsf.data.history(trial,:) = [trial VF contrast response]; % updating the experimental history
                            qcsf = runQCSF(qcsf, 'posttrial', VF, contrast, response);
                        else
                            fprintf('Trial %i is going to be repeated.', trial)
                        end
                        if ~this.DEBUG_WTS
                            Screen('Close');
                        end
                        if breakExperiment
                            break;
                        end
                    end

                    output = struct('qcsf', qcsf, 'stim', this.stim,'display', this.display, 'sID', ...
                        this.sID, 'time', this.time, 'redo_trials', redoTrials);
                    if strcmp(this.fixFreq, 'temporal')
                        data(["TF_" + this.FF(i)]) = output;
                    elseif strcmp(this.fixFreq, 'spatial')
                        data(["SF_" + this.FF(i)]) = output;
                    end

                    % if quitExperiment
                    %     break;
                    % end

                end

                this.time.endExp = datestr(now);
                diary('off');
                Screen('CloseAll');
                ListenChar(0);
                ShowCursor;
                if ~this.DEBUG_WTS
                    BitsPlusPlus('Close');
                end
                outputdir = [this.homedir filesep 'output' filesep this.sID];
                if ~(isfolder(outputdir))
                    mkdir(outputdir)
                end
                save([append(outputdir, filesep, this.sID, '_qCSF_', this.stim.condition, '_')  datestr(now, 'yy-mm-dd_HH-MM') '.mat'], 'data');

            catch ME
                this.time.endExp = datestr(now);
                diary('off');
                Screen('CloseAll');
                ListenChar(0);
                ShowCursor;
                if ~this.DEBUG_WTS
                    BitsPlusPlus('Close');
                end
                outputdir = [this.homedir filesep 'output' filesep this.sID];
                if ~(isfolder(outputdir))
                    mkdir(outputdir)
                end
                save([append(outputdir, filesep, this.sID, '_qCSF_', this.stim.condition, '_')  datestr(now, 'yy-mm-dd_HH-MM') '_error' '.mat'], 'data');
                rethrow(ME);

            end
        end

        function [this] = ImportEyeMovements(this)
            this.eyeData = ImportFixationData(this.display.fRate, append(cd, "/Data/eyemovements/"), 0 ,this.jitter_flag);
        end
    end


    methods (Access = private)
        %% Setting the required dependent library paths
        function retCode = SetDependentLibraries(this)
            try
                addpath(genpath('./qCSF_v1/'));
                addpath(genpath('../UWToolbox'));
                addpath(genpath('./utilities/'));
                retCode = 1;
            catch
                retCode = 0;
            end
        end

        %% Function to set the response keys from Keyboard input
        function this = setKbSettings(this)
            KbName('UnifyKeyNames');
            this.escapeKey = KbName('ESCAPE');
            this.upKey = KbName('UpArrow');
            this.downKey = KbName('DownArrow');
            this.leftKey = KbName('LeftArrow');
            this.rightKey = KbName('RightArrow');
            this.responseKey = KbName('Space');
        end

        %% Setting the Experiment variables
        function this = SetExptVariables(this, photoswitch)

            [this.keybIndices, this.keybNames] = GetKeyboardIndices;
            this.display.pxPerDeg = angle2pix(this.display,1);
            this.display.anglePerPix = pix2angle(this.display,1);

            % t/s CSF settings
            if length(this.type) >= 2
                this.stim.numFrames = round(this.display.fRate*this.stim.dur/1000);
                slope = 16; tails = this.stim.numFrames/slope;
                intermediate = this.stim.numFrames-2*tails;
                if mod(intermediate,2) == 0
                    this.stim.ramp = [linspace(0,1,tails) ones(1, this.stim.numFrames-2*tails) linspace(1, 0, tails)];
                else
                    this.stim.ramp = [linspace(0,1,tails) ones(1, this.stim.numFrames-2*tails + 1) linspace(1, 0, tails)];
                end

                if strcmp(this.type{1,1},'opto' ) || strcmp(this.type{1,1}, 'eye')
                    this = SetOptoProperties(this, photoswitch);
                end

                % fixing the TF fequencies based on the experiment condition
                if this.type{2,1} == 'v1'
                    this.fixFreq = 'temporal';
                    if this.DEBUG
                        this.stim.FF =  [3 10];
                        this.stim.numTrialsPerFreq = 3;
                        this.stim.seed = rng(230).Seed;
                        this.FF = this.stim.FF;
                    else
                        this.stim.seed = rng('shuffle').Seed;
                        this.stim.FF =  30*(1/2).^(4:-1:0); %[3.75 7.5 15];
                        this.stim.numTrialsPerFreq = 50;
                        this.FF = this.stim.FF(randperm(length(this.stim.FF)));
                    end
                    this.numTrials = this.stim.numTrialsPerFreq*length(this.stim.FF);

                elseif this.type{2,1} == 'v2'
                    % fixing the SF fequencies based on the experiment condition
                    this.fixFreq = 'spatial';
                    if this.DEBUG
                        this.stim.FF =  [1, 5];
                        this.stim.numTrialsPerFreq = 6;
                        this.stim.seed = rng(230).Seed;
                        this.FF = this.stim.FF;
                    else
                        this.stim.FF = exp(linspace(log(0.5), log(16),5));
                        % this.stim.FF = [0.5000,    1.4565,    4.2426,   12.3586, 30];
                        this.stim.seed = rng('shuffle').Seed;
                        this.stim.numTrialsPerFreq = 50;
                        this.FF = this.stim.FF(randperm(length(this.stim.FF)));
                    end

                    this.numTrials = this.stim.numTrialsPerFreq*length(this.stim.FF);
                else
                    sprintf("Please enter a valid condition. Thanks!")
                end
            end
        end

        %% Display screen settings for the experiment
        function this = SetScreenSettings(this)
            screenRes = Screen('Rect', this.display.screen);

            %% Open the screen and prepare stuff
            try
                % DISPLAY++ DEFAULTS
                if ~this.DEBUG_WTS
                    PsychDefaultSetup(2);
                    PsychImaging('PrepareConfiguration');
                    PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');
                    PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'ClampOnly');
                    PsychImaging('AddTask', 'General', 'EnableBits++Mono++Output');
                    BitsPlusPlus('OpenBits#');
                end
                if this.DARK
                    if strcmp(this.type{1,1},'opto') || strcmp(this.type{1,1}, 'eye')
                        this.display.gray = this.p.delta;
                    else
                        this.display.gray = 0.202;
                    end
                else
                    this.display.gray = WhiteIndex(this.display.screen)/2;
                end

                if ~this.DEBUG_WTS
                    [this.theScreen, theScreenArea] = PsychImaging('OpenWindow',this.display.screen,this.display.gray, [0 0 screenRes(3) screenRes(4)]);

                    % Generally recommended for graphics
                    Screen('BlendFunction', this.theScreen, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                    ifi = Screen('GetFlipInterval', this.theScreen);
                else
                    ifi = 1/this.display.fRate;
                end
                MaxTime = this.stim.dur/1000; % Set maximum time for all stimuli to be presented in seconds
                %                 numFrames = round(MaxTime/ifi); % Number of frames per trial
                this.t = 1/this.display.fRate :1/this.display.fRate :MaxTime;  % asssume 120Hz
                this.stim.t = this.t;
                this.stim.dt = 1/this.display.fRate;
                % throw up an errr if the frame rate is far from 120.
                if abs(this.display.fRate - 1/ifi)>1
                    error('Frame rate is %5.2f Hz',1/ifi);
                end

                if ~this.DEBUG_WTS
                    % Text settings
                    text.size = 40;
                    text.color = 1;

                    % Setting the screen text size

                    Screen('TextSize', this.theScreen, text.size);

                    % Flip (show) the screen
                    this.flipHandle = Screen('Flip', this.theScreen);

                    % Hide cursor
                    HideCursor;
                    ListenChar(2);
                end


            catch ME
                disp('Error in opening the screen.')
                Screen('CloseAll');
                rethrow(ME);
            end
        end

        %% Setting the audio frequencies for feedback sound
        function this = SetFbSoundSettings(this, Hi, Lo, Cue)
            % Feedback sound settings
            this.beepHigh = MakeBeep(2000,Hi) * 0.1;
            this.beepLow = MakeBeep(500,Lo) * 0.1;
            this.beepBegin = MakeBeep(1000,Cue)*0.1;
        end

        %% Function to play the feedback sound
        function [returnCode, response] = PlayFbSounds(this, keyPressed, CR)
            if nargin < 2
                CR = 0;
            end
            sprintf(keyPressed);
            if strcmp(keyPressed, '0')
                Snd('Play',this.beepBegin,44100);

                returnCode = 1;
            elseif strcmp(keyPressed, 'LeftArrow') || strcmp(keyPressed, '4')
                if CR == 1 % -45 degrees and correct
                    response = 1;
                    Snd('Play',this.beepHigh,44100);
                    returnCode = 1;
                elseif CR == 2
                    response = 0; % -45 and incorrect
                    Snd('Play',this.beepLow,44100);
                    returnCode = 1;
                end
            elseif strcmp(keyPressed, 'RightArrow') || strcmp(keyPressed, '6')
                if CR == 2 % 45 and correct
                    response = 1;
                    Snd('Play',this.beepHigh,44100);
                    returnCode = 1;
                elseif CR == 1
                    response = 0; % 45 and incorrect

                    Snd('Play',this.beepLow,44100);
                    returnCode = 1;
                end
            elseif strcmp(keyPressed, 'q') || strcmp(keyPressed, 'ESCAPE')
                fprintf("Experiment quit by the user; Pressed %s key \n", keyPressed);
                returnCode = 0;
                response = NaN;
            else
                fprintf('Experiment broken by the user; Pressed an invalid key\n');
                returnCode = -1;
                response = NaN;
            end
        end

        %% QCSF settings based on v1/v2
        function qcsf = SetQCSF(this)
            %% quick CSF settings
            % Set priors and setupQCSF
            % (1) peak gain
            % (2) peak frequency
            % (3) bandwidth (full width at half maximum - FWHM - in octaves
            % (4) truncation level on the low-frequency side.

            if strcmp(this.fixFreq, 'temporal')
                priors = [75 1 2.5 .1];
                qcsf = setupQCSF;
            elseif strcmp(this.fixFreq, 'spatial')
                priors = [75 5 2.5 .1];
                qcsf = setupQtCSF;
            else
                disp("Error in setupQCSF: Please check your argument, \n " + ...
                    "Options: 'spatial' or 'temporal'. Bye!");
                qcsf = {};
                return;
            end
            qcsf = setupPriors(qcsf,priors);
        end

        %% Defining input parameters based on the stimulus type such as grating, movie etc
        function [this] = DefineInputParams(this, type, filename)
            if nargin <= 2
                filename = '';
            end

            if strcmp(type, 'grating')
                % Generate a meshgrid
                [x,y] = meshgrid(linspace(-this.stim.size(1)/2,...
                    this.stim.size(1)/2,...
                    angle2pix(this.display,this.stim.size(1))),...
                    linspace(-this.stim.size(2)/2,...
                    this.stim.size(2)/2,...
                    angle2pix(this.display,this.stim.size(2))));
                param.x = x;
                param.y = y;
                Gauss = exp(-((x-this.stim.center(1)).^2+(y-this.stim.center(2)).^2)/(2*this.stim.sigma^2));
                this.input = struct('type', type, ...
                    'params', param, ...
                    'aperture', Gauss);
            elseif strcmp(type, 'movie')
                if ~strcmp(filename, '')
                    this.input = struct('type', type, 'params', filename);
                end

            end
        end

        %% Generate an input grating for a given SF and TF
        function [stimulus] = GenerateGrating(this, SF, TF, contrast, orient)

            ramp = cos(orient*pi/180)*this.input.params.x + sin(orient*pi/180)*this.input.params.y;

            grating = contrast*cos(2*pi*ramp*SF).*this.input.aperture;
            tt = cos(2*pi*this.t*TF)'.*this.stim.ramp'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).

            stimulus = tt*grating(:)';

            if this.DEBUG_STIMCHANGE
                %                 figure(5)
                %                 imagesc(ramp);
                %                 figure(6)
                %                 imagesc(grating);
                %                 figure(7)
                %                 plot(this.t, tt);
            end
        end

        function [stimulus] = GenerateEMGrating(this, SF, TF, contrast, orient, eyeData)
            stimulus = [];
            ramp = cos(orient*pi/180)*this.input.params.x + sin(orient*pi/180)*this.input.params.y;
            first_frame = 0; eye_iter = 1;
            last_frame = eyeData(3, 1);

            theta = 2*pi*SF;
            scaleFac = contrast*(0.998);

            while (last_frame < length(this.t) || (first_frame == 0))
                eye_iter = eye_iter+1;
                if eye_iter < size(eyeData,2)
                    last_frame = min(length(this.t), first_frame+eyeData(3, eye_iter));
                else
                    last_frame = length(this.t);
                end

                %                 if eyeData(3, eye_iter) > 2
                %                 scaleFac = contrast;
                if last_frame > 2
                    if orient == 45
                        orient_idx = 4;
                    else
                        orient_idx = 5;
                    end
                    grating1 = scaleFac*cos(theta*(ramp - eyeData(orient_idx,eye_iter)/3)).*this.input.aperture;
                    grating2 = scaleFac*cos(theta*(ramp - 2*eyeData(orient_idx,eye_iter)/3)).*this.input.aperture;
                    grating3 = scaleFac*cos(theta*(ramp - eyeData(orient_idx,eye_iter))).*this.input.aperture;

                    tt1 = (cos(2*pi*TF*(this.t(first_frame+1:last_frame-2))).*this.stim.ramp(first_frame+1:last_frame-2))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).
                    tt2 = (cos(2*pi*TF*(this.t(last_frame-1))).*this.stim.ramp(last_frame-1))';
                    tt3 = (cos(2*pi*TF*(this.t(last_frame))).*this.stim.ramp(last_frame))';


                    stimulus(first_frame+1 : last_frame-2, :) = tt1*grating1(:)';
                    stimulus(last_frame-1, :) = tt2*grating2(:)';
                    stimulus(last_frame, :) = tt3*grating3(:)';

                else
                    if orient == 45
                        orient_idx = 4;
                    else
                        orient_idx = 5;
                    end
                    grating = scaleFac*cos(theta*(ramp - eyeData(orient_idx,eye_iter))).*this.input.aperture;
                    tt = (cos(2*pi*TF*(this.t(first_frame+1:last_frame))).*this.stim.ramp(first_frame+1:last_frame))';
                    stimulus(first_frame+1 : last_frame, :) = tt*grating(:)';
                end
                first_frame = last_frame;
            end
        end

        function [stimulus] = GenerateEMGrating2(this, SF, TF, contrast, orient, eyeData)
            stimulus = [];
            ramp = cos(orient*pi/180)*this.input.params.x + sin(orient*pi/180)*this.input.params.y;
            first_frame = 0; eye_iter = 1;
            last_frame = eyeData(3, 1);

            theta = 2*pi*SF;
            base1 = [1;1]/(theta); % for movement along +45 deg line - for orient = -45
            base2 = [1; -1]/(theta); % for movement along -45 deg line - for orient = +45
            scaleFac = contrast*(0.998);
            temp = scaleFac*cos((2*pi*SF)*(ramp)).*this.input.aperture;

            while (last_frame < length(this.t) || (first_frame == 0))
                eye_iter = eye_iter+1;
                if eye_iter < size(eyeData,2)
                    last_frame = min(length(this.t), first_frame+eyeData(3, eye_iter));
                else
                    last_frame = length(this.t);
                end
                
                %                 if eyeData(3, eye_iter) > 2
                %                 scaleFac = contrast;
                
                if orient == 45
                    %                         grating1 = scaleFac*cos((2*pi*SF)*(ramp - eyeData(4,eye_iter)/3)).*this.input.aperture;

                    if last_frame > 3
                        movex = FindProjection(eyeData(4,eye_iter-1)*base2, 0);
                        movey = FindProjection(eyeData(4,eye_iter-1)*base2, inf);
                        grating1 = imtranslate(temp, [movex, movey]);
                        tt1 = (cos(2*pi*TF*(this.t(first_frame+1 : last_frame-3))).*this.stim.ramp(first_frame+1 : last_frame-3))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                        stimulus(first_frame+1 : last_frame-3, :) = tt1*grating1(:)';

                        for i = 1:3
                            movex = FindProjection((i/3)*eyeData(4,eye_iter)*base2, 0);
                            movey = FindProjection((i/3)*eyeData(4,eye_iter)*base2, inf);
                            grating = imtranslate(temp, [movex, movey]);
                            tt = (cos(2*pi*TF*(this.t(last_frame-3+i))).*this.stim.ramp(last_frame-3+i))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                            stimulus(last_frame-3+i, :) = tt*grating(:)';

                        end
                    else


                        movex = FindProjection(eyeData(4,eye_iter-1)*base2, 0);
                        movey = FindProjection(eyeData(4,eye_iter-1)*base2, inf);
                        grating1 = imtranslate(temp, [movex, movey]);

                        movex = FindProjection(eyeData(4,eye_iter)*base2, 0);
                        movey = FindProjection(eyeData(4,eye_iter)*base2, inf);
                        grating2 = imtranslate(temp, [movex, movey]);

                        tt1 = (cos(2*pi*TF*(this.t(first_frame+1 : last_frame-1))).*this.stim.ramp(first_frame+1 : last_frame-1))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                        tt2 = (cos(2*pi*TF*(this.t(last_frame))).*this.stim.ramp(last_frame))';
                        stimulus(first_frame+1 : last_frame-1, :) = tt1*grating1(:)';
                        stimulus(last_frame, :) = tt2*grating2(:)';

                    end

                else
                    %                         grating1 = scaleFac*cos((2*pi*SF)*(ramp - eyeData(5,eye_iter)/3)).*this.input.aperture;
                    %                         grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(5,eye_iter)/3)).*this.input.aperture;

                    if last_frame > 3
                        movex = FindProjection(eyeData(5,eye_iter-1)*base1, 0);
                        movey = FindProjection(eyeData(5,eye_iter-1)*base1, inf);
                        grating1 = imtranslate(temp, [movex, movey]);
                        tt1 = (cos(2*pi*TF*(this.t(first_frame+1 : last_frame-3))).*this.stim.ramp(first_frame+1 : last_frame-3))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                        stimulus(first_frame+1 : last_frame-3, :) = tt1*grating1(:)';


                        for i = 1:3
                            movex = FindProjection((i/3)*eyeData(5,eye_iter)*base1, 0);
                            movey = FindProjection((i/3)*eyeData(5,eye_iter)*base1, inf);

                            grating = imtranslate(temp, [movex, movey]);
                            tt = (cos(2*pi*TF*(this.t(last_frame-3+i))).*this.stim.ramp(last_frame-3+i))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                            stimulus(last_frame-3+i, :) = tt*grating(:)';
                        end
                    else

                        movex = FindProjection(eyeData(5,eye_iter-1)*base1, 0);
                        movey = FindProjection(eyeData(5,eye_iter-1)*base1, inf);
                        grating1 = imtranslate(temp, [movex, movey]);

                        movex = FindProjection(eyeData(5,eye_iter)*base1, 0);
                        movey = FindProjection(eyeData(5,eye_iter)*base1, inf);
                        grating2 = imtranslate(temp, [movex, movey]);

                        tt1 = (cos(2*pi*TF*(this.t(first_frame+1 : last_frame-1))).*this.stim.ramp(first_frame+1 : last_frame-1))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                        tt2 = (cos(2*pi*TF*(this.t(last_frame))).*this.stim.ramp(last_frame))';
                        stimulus(first_frame+1 : last_frame-1, :) = tt1*grating1(:)';
                        stimulus(last_frame, :) = tt2*grating2(:)';
                    end

                end

                %                     tt1 = (sin(2*pi*TF*(this.t(first_frame+1:last_frame-2))).*this.stim.ramp(first_frame+1:last_frame-2))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).
                %                     tt2 = (sin(2*pi*TF*(this.t(last_frame-1))).*this.stim.ramp(last_frame-1))';
                %                     tt3 = (sin(2*pi*TF*(this.t(last_frame))).*this.stim.ramp(last_frame))';


                %                     stimulus(first_frame+1 : last_frame-2, :) = tt1*grating1(:)';
                %                     stimulus(last_frame-1, :) = tt2*grating2(:)';
                %                     stimulus(last_frame, :) = tt3*grating3(:)';
                %                 stimulus(first_frame+1 : last_frame, :) = tt3*grating3(:)';
                first_frame = last_frame;
                %                 else
                %                     if orient == 45
                %                         grating = scaleFac*cos((2*pi*SF)*(ramp - eyeData(4,eye_iter))).*this.input.aperture;
                %                     else
                %                         grating = scaleFac*cos((2*pi*SF)*(ramp - eyeData(5,eye_iter))).*this.input.aperture;
                %                     end

                %                     tt = (sin(2*pi*TF*(this.t(first_frame+1:last_frame))).*this.stim.ramp(first_frame+1:last_frame))';
                % %                     tt = reshape(tt, length(tt), 1);
                %                     stimulus(first_frame+1 : last_frame, :) = tt*grating(:)';
                %                 end
            end
        end

        %% Setting the optogenetic parameters for a specified opto-protein model
        function [this] = SetOptoProperties(this, model)
            %% optogenetic model settings
            if strcmp(model, '4_BGAG_12_460_SNAPmGluR')
                this.p = struct('fitModel', model, ...
                    'baseline', 0, ...
                    'tau_on', 0.05, ...
                    'tau_off', 0.3, ...
                    'ampFac', 11, ... % setting this value to make range of baseline == range of opto
                    'bFac', 0.7, ...
                    'tau_b', 0.5, ...
                    'delta', 0.202, ...
                    'offset', 0.9899, ...
                    'opto_black', 1);

            end
        end

        %% Function to apply optogenetic filter on the given stimulus input
        function [opto_stim] = ApplyOptoFilter(this, stimulus, TF, contrast)
            opto_stim = stimulus;%this.p.delta*(zeros(size(stimulus))+this.p.offset);
            % Find the varying pixels.
            gv = ((1/(contrast*TF))*var(stimulus)>0.00001);
            % Pull out the subset of the stimulus movie containing only varying pixels.
            G = stimulus(:,gv);
            % Zero out model variables, size is only for the varying pixels.
            opto_out = zeros([length(this.t),sum(gv)]);
            b = zeros([length(this.t),sum(gv)]);
            dy = zeros(1,sum(gv));

            for j=1:length(this.t)
                if ~this.p.opto_black
                    opto_out(j+1,:) = opto_out(j,:) + this.stim.dt*(this.p.bFac*this.p.ampFac*G(j,:) - (opto_out(j,:)-this.p.baseline)/this.p.tau_b);
                else
                    % stimulus is on
                    idup = G(j,:)>0;
                    dy(idup) = this.stim.dt*(-this.p.ampFac*G(j,idup)+...
                        (this.p.baseline-opto_out(j,idup))/this.p.tau_on);

                    % stimulus is off
                    dy(~idup) = this.stim.dt*(-this.p.ampFac*G(j,~idup)+...
                        (b(j,~idup)-opto_out(j,~idup))/this.p.tau_off);

                    b(j+1,:) = b(j,:) + this.stim.dt*(this.p.bFac*this.p.ampFac*G(j,:) - (b(j,:)-this.p.baseline)/this.p.tau_b);
                    opto_out(j+1,:) = opto_out(j,:)+dy;

                    % For real time computation, co2*pi*mpute each frame by uncommenting
                    % this line:
                    % opto_stim(j,gv) = this.p.delta*(opto_out(j,:)+this.p.offset);
                    opto_stim(j,gv) = opto_out(j,:);
                end

            end
        end


        function [scaled_stimulus] = ScaleStimulus(this, stimulus, scaleFac, offset, contrast, TF)
            %gv = ((1/(contrast*TF))*var(stimulus)>0.001); %0.000001);
            % %Pull out the subset of the stimulus movie
            %G = stimulus(:,gv);
            G = stimulus;
            G = scaleFac*((G+offset));
            scaled_stimulus = scaleFac*(zeros(size(stimulus))+offset);

            for frame=1:length(this.t)
                scaled_stimulus(frame,:) = G(frame,:);
            end
        end

        function [Tex] = MakeTex(this, stimulus, stimulus2, stimulus3)
            if ~exist('stimulus2', 'var')
                stimulus2 = [];
            end
            if ~exist('stimulus3', 'var')
                stimulus3 = [];
            end
            Tex = [];
            if this.DEBUG_STIMCHANGE
                slice = ceil(size(this.input.aperture,1)/2);
                xtemp = linspace(0, 1, size(this.input.aperture,2));
            end

            for frame=1:length(this.t)
                img = reshape(stimulus(frame,:),size(this.input.aperture));
                maximg = max(stimulus(frame,:))*ones(size(xtemp));
                if ~isempty(stimulus2)
                    img2 = reshape(stimulus2(frame,:),size(this.input.aperture));
                    maximg2 = max(stimulus2(frame,:))*ones(size(xtemp));
                end
                if ~isempty(stimulus3)
                    img3 = reshape(stimulus3(frame,:),size(this.input.aperture));
                    maximg3 = max(stimulus3(frame,:))*ones(size(xtemp));
                end
                if this.DEBUG_STIMCHANGE
                    if this.DEBUG_PLOTFFT
                        x = img(slice,:);
                        N = length(x); Fs = 1000; %twice normalized freq
                        xdft = fft(x);
                        xdft = xdft(1:N/2+1);
                        psdx = (1/(Fs*N)) * abs(xdft).^2;
                        psdx(2:end-1) = 2*psdx(2:end-1);
                        freq(frame,:) = 0:Fs/length(x):Fs/2;
                        psdMat(frame,:) = pow2db(psdx);
                        bp(frame) = bandpower(img(slice,:));
                    end


                    figure(1)
                    if ~isempty(stimulus2) && ~isempty(stimulus3)

                        plot(xtemp, img(slice, :), 'Color', [0.8500 0.3250 0.0980]);
                        hold on;
                        plot(xtemp, maximg, 'Color', [0.8500 0.3250 0.0980]);
                        hold on;
                        plot(xtemp, img2(slice, :), 'Color',[0.4940 0.1840 0.5560]);
                        hold on;
                        plot(xtemp, maximg2,  'Color',[0.4940 0.1840 0.5560]);
                        hold on;
                        plot(xtemp, img3(slice, :), 'Color',[0.4660 0.6740 0.1880]);
                        hold on;
                        plot(xtemp, maximg3, 'Color',[0.4660 0.6740 0.1880]);
                        hold off;
                        legend('ctl','max-ctl', 'opto-only', 'max-opyo-only', 'opto+eye-mvt-jitter', 'max-opto+eye-mvt-jitter');

                    elseif ~isempty(stimulus2)
                        plot(xtemp, img(slice, :), 'b');
                        hold on;
                        plot(xtemp, maximg, 'b');
                        hold on;
                        plot(xtemp, img2(slice, :), 'r');
                        hold on;
                        plot(xtemp, maximg2, 'r');
                        hold off;
                        legend('opto-only','max-opto-only', 'opto+eye-mvt', 'max-opto+eye-mvt');
                    else
                        plot(xtemp, img(slice, :));
                    end
                    ylim([this.display.gray-1 ,this.display.gray+1]);
                    grid on;
                    pause(0.1);
                    %                     title('Stimulus SF = 2; TF = 5; orig_contrast = 1', 'Interpreter','none')
                    %                     xlabel('Normalized spatial frequency');
                    %                     ylabel('Contrast after opto-filtering')
                    drawnow;
                    drawframes{frame} = getframe(gcf);
                end

                if ~this.DEBUG_WTS
                    Tex(frame) = Screen('MakeTexture', this.theScreen, img(:,:),[],[],2);
                end
            end

            if this.DEBUG_RECORDMOVIES
                obj = VideoWriter("movies/optoPlusEye");
                obj.Quality = 100;
                obj.FrameRate = 10;
                open(obj);
                for i=1:length(drawframes)
                    writeVideo(obj,drawframes{i})
                end
                obj.close();
            end

            if this.DEBUG_PLOTFFT

                figure(3)
                plot3(this.t, freq,psdMat)
                grid on
                title("Periodogram Using FFT")
                zlabel("Time in s")
                ylabel("Frequency (Hz)")
                zlabel("Power/Frequency (dB/Hz)")

                unbiased_bp = bp - mean(bp);
                totAvg_bp = (mean(bp))*ones(length(this.t),1);
                figure(2)
                %                 plot(this.t, unbiased_bp);
                hold on;
                plot(this.t, totAvg_bp);
                grid on;
                ylim([0, 0.05])
                title('Average band power')
                xlabel('Across frames');
                ylabel('Across frequencies');
            end
        end


        function [eye_stim] = CorrectEyeMovement(this, opto_stim, eyeData, SF, orient)
            
            DEBUG_JITTER = 0;
            eye_iter = 1;
            last_frame = eyeData(3, 1);
%             eye_stim = zeros(length(this.t), size(this.input.aperture,1), size(this.input.aperture,2));
            TF = 0;
            eye_stim = reshape(opto_stim,[],size(this.input.aperture,1), size(this.input.aperture,2));
            ramp = cos(orient*pi/180)*this.input.params.x + sin(orient*pi/180)*this.input.params.y;
            [m] = max(opto_stim, [], 2);
            [m, i] = max(m);
            grating = squeeze(eye_stim(i,:,:));%cos(2*pi*SF*(ramp)-pi).*this.input.aperture;
%             tt = (cos(2*pi*TF*(this.t)).*this.stim.ramp)';

            for frame = 1:length(this.t)
                if ~DEBUG_JITTER
%                     eye_stim(frame, :,:) = reshape(opto_stim(frame, :), size(this.input.aperture));
                    if frame ==  1
%                         grating = squeeze(eye_stim(frame, :,:));
                        first_fft = do_fft(grating);
                    elseif (frame >= last_frame - 2) || frame <= (last_frame)
                        grating = squeeze(eye_stim(frame, :,:));
%                         grate_fft = do_fft(grating);
                        temp = changephase_fft(grating, first_fft);
                        eye_stim(frame,:,:) = temp;
                        if frame == last_frame
%                             first_frame = last_frame;
                            last_frame = eyeData(3,eye_iter);
                            eye_iter = eye_iter + 1;
                        end
                    end

                end
            end

        end

                %% Function to correct for eye movements - saccadic suppression
        function [eye_stim] = CorrectEyeMovement2(this, opto_stim2, eyeData, SF, orient)
            % aperture = this.p.delta*(this.input.aperture+this.p.offset);
            first_frame = 0; eye_iter = 1;
            last_frame = eyeData(3, 1);
%             eye_stim = zeros(size(opto_stim2));
            eye_stim = reshape(opto_stim2,[],size(this.input.aperture,1), size(this.input.aperture,2));
            
            theta = 2*pi*SF;
            base1 = [1;1]/(-theta); % for movement along +45 deg line - for orient = -45
            base2 = [1; -1]/(-theta); % for movement along -45 deg line - for orient = +45

            for frame=1:length(this.t)
                if frame == last_frame
                    eye_iter = eye_iter+1;
                    if eye_iter < size(eyeData,2)
                        last_frame = min(length(this.t), first_frame+eyeData(3, eye_iter));
                    else
                        last_frame = length(this.t);
                    end
                end


                if orient == 45
                    if last_frame > 3
                        temp = eye_stim(first_frame+1:last_frame-3,:,:);
                        movex = FindProjection(eyeData(4,eye_iter-1)*base2, 0);
                        movey = FindProjection(eyeData(4,eye_iter-1)*base2, inf);
                        eye_stim(first_frame+1:last_frame-3,:,:) = imtranslate(eye_stim(frame,:,:), [movex, movey]);

                        for i = 1:3
                            movex = FindProjection((i/3)*eyeData(4,eye_iter)*base2, 0);
                            movey = FindProjection((i/3)*eyeData(4,eye_iter)*base2, inf);
                            temp = opto_stim(last_frame-3+i,:,:);
                            grating = imtranslate(temp, [movex, movey]);
                        end
                    else

                        movex = FindProjection(eyeData(4,eye_iter-1)*base2, 0);
                        movey = FindProjection(eyeData(4,eye_iter-1)*base2, inf);
                        temp = opto_stim(last_frame-1,:,:);
                        grating1 = imtranslate(temp, [movex, movey]);

                        movex = FindProjection(eyeData(4,eye_iter)*base2, 0);
                        movey = FindProjection(eyeData(4,eye_iter)*base2, inf);
                        temp = opto_stim(last_frame,:,:);
                        grating2 = imtranslate(temp, [movex, movey]);

                    end

                else
                    if last_frame > 3
                        movex = FindProjection(eyeData(5,eye_iter-1)*base1, 0);
                        movey = FindProjection(eyeData(5,eye_iter-1)*base1, inf);
                        temp = opto_stim(last_frame-3,:,:);
                        grating1 = imtranslate(temp, [movex, movey]);
                        tt1 = (sin(2*pi*TF*(this.t(first_frame+1 : last_frame-3))).*this.stim.ramp(first_frame+1 : last_frame-3))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                        stimulus(first_frame+1 : last_frame-3, :) = tt1*grating1(:)';


                        for i = 1:3
                            movex = FindProjection((i/3)*eyeData(5,eye_iter)*base1, 0);
                            movey = FindProjection((i/3)*eyeData(5,eye_iter)*base1, inf);
                            temp = opto_stim(last_frame-3+i,:,:);
                            grating = imtranslate(temp, [movex, movey]);
                            tt = (sin(2*pi*TF*(this.t(last_frame-3+i))).*this.stim.ramp(last_frame-3+i))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                            stimulus(last_frame-3+i, :) = tt*grating(:)';
                        end
                    else

                        movex = FindProjection(eyeData(5,eye_iter-1)*base1, 0);
                        movey = FindProjection(eyeData(5,eye_iter-1)*base1, inf);
                        temp = opto_stim(last_frame-1,:,:);
                        grating1 = imtranslate(temp, [movex, movey]);

                        movex = FindProjection(eyeData(5,eye_iter)*base1, 0);
                        movey = FindProjection(eyeData(5,eye_iter)*base1, inf);
                        temp = opto_stim(last_frame,:,:);
                        grating2 = imtranslate(temp, [movex, movey]);

                        tt1 = (sin(2*pi*TF*(this.t(first_frame+1 : last_frame-1))).*this.stim.ramp(first_frame+1 : last_frame-1))'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).grating2 = scaleFac*cos((2*pi*SF)*(ramp - 2*eyeData(4,eye_iter)/3)).*this.input.aperture;
                        tt2 = (sin(2*pi*TF*(this.t(last_frame))).*this.stim.ramp(last_frame))';
                        stimulus(first_frame+1 : last_frame-1, :) = tt1*grating1(:)';
                        stimulus(last_frame, :) = tt2*grating2(:)';
                    end

                end
                first_frame = last_frame;
            end


        end

        %% Function to correct for eye movements - saccadic suppression
        function [eye_stim] = CorrectEyeMovement3(this, opto_stim, eyeData, SF, orient)
            % aperture = this.p.delta*(this.input.aperture+this.p.offset);
            first_frame = 0; eye_iter = 1;
            last_frame = eyeData(3, 1);
            eye_stim = zeros(size(opto_stim));
            if this.DEBUG_STIMCHANGE
                orient = 0;
            end
            orient_rad = pi*orient/180;
            theta = 2*pi*SF;
            base1 = [1;1]/(-theta); % for movement along +45 deg line - for orient = -45
            base2 = [1; -1]/(-theta); % for movement along -45 deg line - for orient = +45
            if orient == 45
                %                 movex = (rem(eyeData(4,1),theta))*cos(orient_rad);
                %                 movey = (rem(eyeData(4,1),theta))*sin(orient_rad);
                movex = FindProjection(eyeData(4,1)*base2, 0);
                movey = FindProjection(eyeData(4,1)*base2, inf);
            else
                %                 movex = (rem(eyeData(5,1),theta))*cos(orient_rad);
                %                 movey = (rem(eyeData(5,1),theta))*sin(orient_rad);
                movex = FindProjection(eyeData(5,1)*base1, 0);
                movey = FindProjection(eyeData(5,1)*base1, inf);
            end
            %             % move_aperture = exp(-((this.input.params.x+movex).^2+(this.input.params.y+movey).^2)/(2*this.stim.sigma^2));

            for frame=1:length(this.t)
                %                 if frame == (last_frame - 2)
                %                     if orient == 45
                %                         movex = (mod((eyeData(4,eye_iter)/3),theta))*cos(orient_rad);
                %                         movey = (mod((eyeData(4,eye_iter)/3),theta))*sin(orient_rad);
                %                     else
                %                         movex = (mod((eyeDatathis(5,eye_iter)/3),theta))*cos(orient_rad);
                %                         movey = (mod((eyeData(5,eye_iter)/3),theta))*sin(orient_rad);
                %                     end
                %                 elseif frame == (last_frame - 1)
                %                     if orient == 45
                %                         movex = (mod((2*eyeData(4,eye_iter)/3),theta))*cos(orient_rad);
                %                         movey = (mod((2*eyeData(4,eye_iter)/3),theta))*sin(orient_rad);
                %                     else
                %                         movex = (mod((2*eyeData(5,eye_iter)/3),theta))*cos(orient_rad);
                %                         movey = (mod((2*eyeData(5,eye_iter)/3),theta))*sin(orient_rad);
                %                     end
                if (frame == (last_frame+1)) && (eye_iter < size(eyeData,2))
                    eye_iter = eye_iter+1;
                    first_frame = last_frame;
                    last_frame = min(length(this.t), first_frame+eyeData(3, eye_iter));
                    if orient == 45
                        %                         movex = (rem(eyeData(4,eye_iter),theta))*cos(orient_rad);
                        %                         movey = (                     rem(eyeData(4,eye_iter),theta))*sin(orient_rad);
                        movex = FindProjection(eyeData(4,eye_iter)*base2, 0);
                        movey = FindProjection(eyeData(4,eye_iter)*base2, inf);
                    else
                        %                         movex = (rem(eyeData(5,eye_iter),theta))*cos(orient_rad);
                        %                         movey = (rem(eyeData(5,eye_iter),theta))*sin(orient_rad);
                        movex = FindProjection(eyeData(5,eye_iter)*base1, 0);
                        movey = FindProjection(eyeData(5,eye_iter)*base1, inf);
                    end
                end
                %
                % temp = (reshape(opto_stim(frame,:),size(move_aperture))).*move_aperture;
                temp = reshape(opto_stim(frame,:),size(this.input.aperture));
                %temp2 = imtranslate(temp, [-1*movetemp2x,movey -1*movey], 'FillValues',this.display.gray-0.002);
                temp2 = imtranslate(temp, [movex, movey]);
                eye_stim(frame, :) = temp2(:);

                % move_aperture = exp(-((this.input.params.x+movex).^2+(this.input.params.y+movey).^2)/(2*this.stim.sigma^2));
            end
        end


    end
end
