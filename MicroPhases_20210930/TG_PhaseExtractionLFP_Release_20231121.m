%==========================================================================
% This script extracts phases of filtered LFP signal using the
% generalized phase approach.
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.filtered  = 'D:\TreasureHunt\MicroFiltered_20210930'; % data path
paths.save      = 'D:\TreasureHunt\MicroPhases_20210930'; % save folder

% functions
addpath(genpath('D:\External\Functions'));
addpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions');

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% parameters
filterBand      = [1, 10]; % frequency band
lp              = 1; % cutoff for negative frequency detection [Hz]
phaseMethod     = 'generalized'; % phase extraction method

% load and save name
loadName        = strcat('datacutFt2000HzBP_', regexprep(num2str(filterBand), '\s+', '_'), '.mat');
saveName        = strcat('datacutFt2000HzBP_', phaseMethod, '_', regexprep(num2str(filterBand), '\s+', '_'), '.mat');

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
    sessions = dir(fullfile(paths.filtered, subjects{iSub}, 'session*'));
    
    % loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % display session information
        fprintf('\n==================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        % get available microwires
        chanDir = TG_GetChanDir_20210812(paths.filtered, subjects{iSub}, sessions(iSess).name);
        
        % loop through channels
        parfor iWire = 1:size(chanDir, 1)

            % print channel name
            fprintf('\tChannel name: %s.\n', chanDir(iWire).name);
            
            % load data
            data            = load(fullfile(paths.filtered, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, loadName));
            
            % preallocate phase data structure
            phaseData       = data;

            % generalized phase of data
            try
                generalizedData = TG_generalized_phase_vector(data.trial{1, 1}', data.fsample, lp);
                phaseData.trial = {transpose(generalizedData)};
            catch
                phaseData.trial = [];
            end

            % save filtered data
            mkdir(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name));
            TG_save(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, saveName), phaseData);
        end
    end
end

% shut down parallel pool
delete(gcp);