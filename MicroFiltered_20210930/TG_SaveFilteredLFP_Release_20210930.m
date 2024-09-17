%==========================================================================
% This script filters downsampled LFP signals and saves them.
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.micro = 'D:\TreasureHunt\MicroDownsampled_20210910'; % data path
paths.save  = 'D:\TreasureHunt\MicroFiltered_20210930'; % save folder

% own functions
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% frequency band and file name
filterBand = [3, 10];
saveName   = strcat('datacutFt2000HzBP_', regexprep(num2str(filterBand), '\s+', '_'), '.mat');

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
        parfor iWire = 1:size(chanDir, 1)

            % print channel name
            fprintf('\tChannel name: %s.\n', chanDir(iWire).name);
            
            % load data
            data           = load(fullfile(paths.micro, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, 'datacutFt2000Hz.mat'));
            
            % set configurations for bandstop filter
            cfg            = [];
            cfg.bsfilter   = 'yes';
            cfg.bsfiltype  = 'fir';
            cfg.bsfreq     = [48, 52; 98, 102; 148, 152];
            
            % bandstop filter data
            bsData         = ft_preprocessing(cfg, data);
            
            % set configurations for bandpass filter
            cfg            = [];
            cfg.bpfilter   = 'yes';
            cfg.bpfilttype = 'fir';
            cfg.bpfreq     = filterBand;
            cfg.demean     = 'yes';
            
            % bandpass filter data
            filteredData   = ft_preprocessing(cfg, bsData);
            
            % save filtered data
            mkdir(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name));
            TG_save(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, saveName), filteredData);
        end
    end
end

% shut down parallel pool
delete(gcp);
