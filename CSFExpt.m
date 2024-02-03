classdef CSFExpt
    %% Public Properties
    properties (Access = public)
        % Configuration flags
        DEBUG = 0;
        preload = 1; waitframes = 1;
        DARK = 1;

        % Stimulus variable(s)
        % centre of the stimulus (in degrees of VA - usually set this to 0)
        % width/height of the stimulus prior to gaussian envelope
        % gaussian envelope over the grating (in degree of visual angle)
        stim = struct('condition', '', 'dur', 800, 'center', [0,0], 'size', [10,10], 'sigma', 1.25, 'numTrialsPerFreq', 0);
        input = {};

        % Display variables: Room 510 - Display++
        % display.width: measure the physical width of your screen in cm;
        %display.fRate: Frame rate of the display - 120Hz for ours
        % dispay.dist: in cm; measure your viewing distance
        display = struct('resolution', [1920 1080], ...
            'width', 71, 'dist', 80, ...
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
            this = SetScreenSettings(this);
            setKbSettings(this);
            this = SetFbSoundSettings(this, 0.02, 0.04, 0.01);
            this = DefineInputParams(this, config.inputType);
            %             this = RunExpt(this);
        end

        %% Routine to run the experiment
        function this = RunExpt(this)
            try

                ListenChar(2);
                ApplyKbFilter;

                % create a data structure to store all levels
                data = dictionary();
                redoTrials = []; % list of trials repeated

                for i= 1: length(this.FF)
                    qcsf = SetQCSF(this);
                    close all
                    % initialize redoTrial flag to avoid skipping frames
                    redoTrial = 0;
                    breakExperiment = 0;
                    for trial = 1:this.stim.numTrialsPerFreq
                        tic
                        %% Go through trial
                        % Play trial start sound
                        PlayFbSounds(this, '0');
                        Screen('FillRect',this.theScreen,this.display.gray)
                        DrawFormattedText(this.theScreen, '.', 'center', 'center', 0);
                        Screen('Flip',this.theScreen)

                        %%
                        if strcmp(this.type{1,1},'opto') || strcmp(this.type{1,1}, 'eye')
                            WaitSecs(0.005);
                        else
                            WaitSecs(0.5);
                        end
                        qcsf.data.trial = trial;
                        [qcsf,VF,contrast] = runQCSF(qcsf, 'pretrial');
                        %% Determine whether contrast stim is -45 or 45 deg (or 0/90, or whatever you end up picking)
                        % CR = ceil(2*rand);
                        CR = randi(2);
                        if CR == 1
                            orient = -45; % in degress
                        else
                            orient = 45;
                        end


                        if (strcmp(this.fixFreq,'temporal'))
                            [stimulus] = GenerateGrating(this, VF, this.FF(i), contrast, orient);
                            if strcmp(this.type{1,1},'baseline')
                                Tex = MakeTex(this, stimulus, this.FF(i), contrast);
                            elseif strcmp(this.type{1,1},'opto')
                                Tex = ApplyOptoFilter(this, stimulus, this.FF(i), contrast);
                            end


                        elseif(strcmp(this.fixFreq,'spatial'))
                            [stimulus] = GenerateGrating(this, this.FF(i), VF, contrast,orient);
                            if strcmp(this.type{1,1},'baseline')
                                Tex = MakeTex(this, stimulus, VF, contrast);
                            elseif strcmp(this.type{1,1},'opto')
                                Tex = ApplyOptoFilter(this, stimulus, VF, contrast); % returns the texture if opto
                            end
                        end


                        if this.preload
                            Screen('PreloadTextures', this.theScreen, Tex);
                            % disp('Textures are preloaded');
                        end


                        % Present qCSF
                        StartT = GetSecs; %Measure start time of session
                        k = 0;
                        frame_time = zeros(length(this.t),1);
                        while (GetSecs() - StartT) < this.stim.dur/1000
                            k = ceil((GetSecs() - StartT) .* this.display.fRate); % current frame based on current time
                            if k <= (this.stim.numFrames)
                                Screen('DrawTextures', this.theScreen, Tex(k));
                                Screen('Flip', this.theScreen); % show to participant
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
                        Screen('FillRect', this.theScreen, this.display.gray);
                        Screen('Flip', this.theScreen);
                        toc

                        tic
                        [keyPressed] = wait4key({'LeftArrow','RightArrow', 'ESCAPE', 'q'}, -3);
                        sprintf(keyPressed);
                        [returnCode, response] = PlayFbSounds(this,keyPressed, CR);
                        if returnCode == 0 && breakExperiment == 1
                            fprintf("Experiment quit by the user; Pressed %s key on trial = %d \n", keyPressed, trial);
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
                        Screen('Close');

                    end

                    output = struct('qcsf', qcsf, 'stim', this.stim,'display', this.display, 'sID', ...
                        this.sID, 'time', this.time, 'redo_trials', redoTrials);
                    if strcmp(this.fixFreq, 'temporal')
                        data(["TF_" + this.FF(i)]) = output;
                    elseif strcmp(this.fixFreq, 'spatial')
                        data(["SF_" + this.FF(i)]) = output;
                    end

                    if breakExperiment
                        break;
                    end
                end

                this.time.endExp = datestr(now);
                diary('off');
                Screen('CloseAll');
                ListenChar(0);
                ShowCursor;
                BitsPlusPlus('Close');

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
                BitsPlusPlus('Close');

                outputdir = [this.homedir filesep 'output' filesep this.sID];
                if ~(isfolder(outputdir))
                    mkdir(outputdir)
                end
                save([append(outputdir, filesep, this.sID, '_qCSF_', this.stim.condition, '_')  datestr(now, 'yy-mm-dd_HH-MM') '_error' '.mat'], 'data');
                rethrow(ME);

            end
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
                this.stim.ramp = [linspace(0,1,tails) ones(1, this.stim.numFrames-2*tails) linspace(1, 0, tails)];

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
                        this.stim.FF =  [7 10];
                        this.stim.numTrialsPerFreq = 3;
                        this.stim.seed = rng(230).Seed;
                        this.FF = this.stim.FF;
                    else
                        this.stim.FF = exp(linspace(log(0.5), log(30),5));
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
                PsychDefaultSetup(2);
                PsychImaging('PrepareConfiguration');
                PsychImaging('AddTask', 'General', 'FloatingPoint32Bit');
                PsychImaging('AddTask', 'FinalFormatting', 'DisplayColorCorrection', 'ClampOnly');
                PsychImaging('AddTask', 'General', 'EnableBits++Mono++Output');
                BitsPlusPlus('OpenBits#');
                
                if this.DARK
                    if strcmp(this.type{1,1},'opto') || strcmp(this.type{1,1}, 'eye')
                        this.display.gray = this.p.delta;
                    else
                        this.display.gray = 0.202;
                    end
                else
                    this.display.gray = WhiteIndex(this.display.screen)/2;
                end


                [this.theScreen, theScreenArea] = PsychImaging('OpenWindow',this.display.screen,this.display.gray, [0 0 screenRes(3) screenRes(4)]);

                % Generally recommended for graphics
                Screen('BlendFunction', this.theScreen, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

                ifi = Screen('GetFlipInterval', this.theScreen);
                MaxTime = this.stim.dur/1000; % Set maximum time for all stimuli to be presented in seconds
                %                 numFrames = round(MaxTime/ifi); % Number of frames per trial
                this.t = 1/this.display.fRate :1/this.display.fRate :MaxTime;  % asssume 120Hz
                this.stim.t = this.t;
                this.stim.dt = 1/this.display.fRate;
                % throw up an errr if the frame rate is far from 120.
                if abs(this.display.fRate - 1/ifi)>1
                    error('Frame rate is %5.2f Hz',1/ifi);
                end

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
                fprintf("Experiment quit by the user; Pressed %s key on trial = %d \n", keyPressed, trial);
                returnCode = 0;
                response = NaN;
            else
                fprintf('Experiment broken by the user; Pressed an invalid key on trial = %d \n', trial);
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

            % ramp
            ramp = cos(orient*pi/180)*this.input.params.x + sin(orient*pi/180)*this.input.params.y;
            grating = contrast*cos(2*pi*ramp*SF).*this.input.aperture;
            tt = sin(2*pi*this.t*TF)'.*this.stim.ramp'; % Use the outer product to get 'stimulus' which is the counterphase grating movie (very fast).
            stimulus = tt*grating(:)';
        end

        %% Setting the optogenetic parameters for a specified opto-protein model
        function [this] = SetOptoProperties(this, model)
            %% optogenetic model settings
            if strcmp(model, '4_BGAG_12_460_SNAPmGluR')
                this.p = struct('fitModel', model, ...
                    'baseline', 0, ...
                    'tau_on', 0.05, ...
                    'tau_off', 0.3, ...
                    'ampFac', 20, ...
                    'bFac', 0.7, ...
                    'tau_b', 0.5, ...
                    'delta', 0.202, ...
                    'offset', 0.9899, ...
                    'opto_black', 1);
            end
            %             this.p.fitModel = '4_BGAG_12_460_SNAPmGluR'; % Holt photoswitch
            %             p.baseline = 0;
            %             p.tau_on = .05;    % 'linear' filter onset time constant
            %             p.tau_off = .3;    % 'linear' filter offset time constant
            %             p.ampFac = 20;    % scale factor for stimulus influence on linear filter
            %             p.bFac = .7; %.05; % scale factor for stimulus influence on 'b'
            %             p.tau_b = .5;      % time constant for recovery of 'b' toward initial baseline
            %             % p.scaleFactor = 4; % scale factor based on getScalingFactor_DISPLAYPP.m run at TF = 0.5
            %
            %             % original gratings were scaled between 0 - 1
            %             % opto response at it's max contrast (low TF) goes between -0.9899 and 3.9596
            %             % rescale them, based on a low TF grating
            %             p.delta = 0.202;
            %             p.offset = 0.9899;
            % this is needed because as the TF increases the
            % amplitude decreasing effect of the optogenetic protein
            % increases. See Watson temporal sensitivity
        end

        %% Function to apply optogenetic filter on the given stimulus input
        function [Tex] = ApplyOptoFilter(this, stimulus, TF, contrast)
            result = this.p.delta*(zeros(size(stimulus))+this.p.offset);
            % Find the varying pixels.
            gv = ((1/(contrast*TF))*var(stimulus)>0.000000);
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

                    % For real time computation, compute each frame by uncommenting
                    % this line:
                    result(j,gv) = this.p.delta*(opto_out(j,:)+this.p.offset);
                end
                img = reshape(result(j,:),size(this.input.aperture));
                Tex(j) = Screen('MakeTexture', this.theScreen, img(:,:));
            end
        end

        function [Tex] = MakeTex(this, stimulus, TF, contrast)
            gv = ((1/(contrast*TF))*var(stimulus)>0); %0.000001);
            % Pull out the subset of the stimulus movie containing only varying pixels.
            G = stimulus(:,gv);
            G = this.display.gray*((G+1));
            gratings = (this.display.gray-0.002)*ones(length(this.t), size(this.input.aperture,1), size(this.input.aperture,2));
            for frame=1:length(this.t)
                gratings(frame,gv) = G(frame,:);
                img = reshape(gratings(frame,:),size(this.input.aperture));
                Tex(frame) = Screen('MakeTexture', this.theScreen, img(:,:));
            end
        end
    end
end

