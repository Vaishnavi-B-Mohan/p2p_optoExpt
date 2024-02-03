
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
        cond_prompt = 'demo';
    case 1
        cond_prompt = 'baseline_v1';
    case 2
        cond_prompt = 'baseline_v2';
    case 3
        cond_prompt = 'opto_v1';
    case 4
        cond_prompt = 'opto_v2';
    case 5
        cond_prompt = 'eye_v1';
    case 6
        cond_prompt = 'eye_v2';
    otherwise
        sprintf("Please enter a valid condition. Thanks!")
        return;

end

config.DEBUG = 1;
config.fRate = 60;
config.photoswitch = '4_BGAG_12_460_SNAPmGluR';
config.inputType = 'grating' ;
expt = CSFExpt(cond_prompt, sID, config);
RunExpt(expt);

sprintf("Lalala")

