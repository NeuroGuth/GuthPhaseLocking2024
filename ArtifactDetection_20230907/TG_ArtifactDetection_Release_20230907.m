%==========================================================================
% This script detects artifacts and interictal epileptiform discharges
% based on the function LK_ArtifactDetection_20220319.m
% (https://github.com/NeuroLuke/KunzNatureNeuroscience2024/blob/main/OpenField/Functions/LK_ArtifactDetection_20220319.m)
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.micro             = 'D:\TreasureHunt\MicroDownsampled_20210910'; % data path
paths.save              = 'D:\TreasureHunt\ArtifactDetection_20230907'; % save folder

% own functions
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% set configurations for artifact detection (base ond amplitudes, gradients and power)
param.time2Exclude     = 1; % additional time to exclude around artifacts (in seconds) (Staresina et al., Nat Neurosci, 2015)
param.amp.idx          = 1; % index for labeling
param.amp.threshFac    = 4; % factor to calculate amplitude threshold
param.gra.idx          = 2; % index for labeling
param.gra.threshFac    = 4; % factor to calculate gradient threshold
param.pow.idx          = 3; % index for labeling
param.pow.threshFac    = 4; % factor to calculate power threshold
param.pow.lowerFreq    = 1; % lower frequency for calculating the power spectrum
param.pow.upperFreq    = 60; % upper frequency for calculating the power spectrum
param.pow.numFreqs     = 30; % number of frequencies for calculating the power spectrum

%% subjects
subjects    = {...
    'TH01'; ...
    'TH02'; ...
    'TH03'; ...
    'TH04'; ...
    'TH05'; ...
    'TH06'; ...
    'TH07'; ...
    'TH08'; ...
    'TH09'; ...
    'TH10'; ...
    'TH11'; ...
    'TH12'; ...
    'TH13'; ...
    'TH14'; ...
    'TH15'; ...
    'TH16'; ...
    'TH17'; ...
    'TH18'; ...
    };

%% save settings
save(fullfile(paths.save, 'settings'));

%% loop through subjects
for iSub = 1:length(subjects)
    
    % get sessions
    sessions = dir(fullfile(paths.micro, subjects{iSub}, 'session*'));
    
    % loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % display session information
        fprintf('\n==================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        % get available microwires
        chanDir = TG_GetChanDir_20210812(paths.micro, subjects{iSub}, sessions(iSess).name);
        
        % loop through channels
        for iWire = 1:size(chanDir, 1)

            % print channel name
            fprintf('\tChannel name: %s.\n', chanDir(iWire).name);
            
            % load raw data
            data                = load(fullfile(paths.micro, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, 'datacutFt2000Hz.mat'));
            
            % artifact detection
            detectedArtifacts   = LK_ArtifactDetection_20220319(param, data);
            
            % save artifact detection results
            mkdir(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name));
            TG_save(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, 'detectedArtifacts.mat'), detectedArtifacts);
        end
    end
end

% shut down parallel pool
delete(gcp);
