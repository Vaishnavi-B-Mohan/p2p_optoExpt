clear all; close all; clc;

addpath(genpath('/qCSF_v1/'));
addpath(genpath('../UWToolbox'));

sID = input('Please enter subject initials? ', 's');
sID = upper(sID);
cond_prompt = input(['Please choose the experiment condition: \n' ...
    '0. Demo experiment \n' ...
    '1. for baseline with fixed TF: baseline_v1 \n' ...
    '2. for baseline with fixed SF: baseline_v2 \n' ...
    '3. for optofilter with fixed TF: opto_v1 \n' ...
    '4. for optofilter with fixed SF: opto_v2 \n' ...
    '5. for opto with eye movements with fixed TF: eye_v1 \n' ...
    '6. for opto with eye movements with fixed SF: eye_v2 \n']);

switch cond_prompt
    case 0
        cond_prompt = 'demo'; % Doesn't work yet. Don't bother
    case 1
        cond_prompt = 'baseline_v1'; % Baseline condition with fixed temporal frequencies
    case 2
        cond_prompt = 'baseline_v2'; % Baseline condition with fixed spatial frequencies
    case 3
        cond_prompt = 'opto_v1'; % Opto condition with fixed 5 temporal frequencies
    case 4
        cond_prompt = 'opto_v2'; % Opto condition with fixed spatial frequencies
    case 5
        cond_prompt = 'eye_v1'; % Opto (with eye-movements) with fixed temporal frequencies
    case 6
        cond_prompt = 'eye_v2'; % Opto (with eye-movements) with fixed spatial frequencies
    otherwise
        sprintf("Please enter a valid condition. Thanks!")
        return;                    

end

config.DEBUG = 1; % Set this to 1 to debug the functionality with only a few trials
config.DEBUG_WTS = 1; % Set this to 1 to debug without the psychtoolbox screen being displayed
config.fRate = 60;
config.jitter_flag = 0; % 0:= Regular eye-movements; 1:= Rapid changes to fixation; 2:= No eye movements
config.DEBUG_STIMCHANGE = 1; % Set this to 1 to observe the sinusoids and associated band power graphs
config.DEBUG_RECORDMOVIES = 0;
config.DEBUG_PLOTFFT = 0;
config.photoswitch = 'ChRmine'; % ["4_BGAG_12_460_SNAPmGluR", "ChR2", "ReaChR", "ChrimsonR", "CsChrimson", "bReaChES", "ChRmine"];
config.inputType = 'grating' ;
expt = CSFExpt(cond_prompt, sID, config);
RunExpt(expt); 

sprintf("Thank you. Goodbye!")
close all;
