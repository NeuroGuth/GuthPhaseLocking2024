%==========================================================================
% This script extracts and saves the instantaneous oscillatory and 
% aperiodic component of the downsampled LFP using SPRiNT (Wilson,
% Castanheira and Baillet, eLife, 2022)
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.load   = 'D:\TreasureHunt\MicroDownsampled_20210910'; % data path
paths.save   = 'D:\TreasureHunt\MicroSPRiNT_20230713'; % save folder

% own functions
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% file name to load
loadName     = 'datacutFt2000Hz.mat';
saveName     = 'datacutSprint.mat';

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
    sessions = dir(fullfile(paths.load, subjects{iSub}, 'session*'));
    
    % loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % display session information
        fprintf('\n==================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        % get available microwires
        chanDir = TG_GetChanDir_20210812(paths.load, subjects{iSub}, sessions(iSess).name);
        
        % loop through channels
        parfor (iWire = 1:size(chanDir, 1), 16)
            
            % print channel name
            fprintf('\tChannel name: %s.\n', chanDir(iWire).name);
            
            % load data
            data                    = load(fullfile(paths.load, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, loadName));
            
            %% SPRiNT processing of data
            
            % general settings
            opt                     = [];
            opt.sfreq               = data.fsample;         % input sampling rate
            opt.WinLength           = 1;                    % STFT window length in seconds; default = 1
            opt.WinOverlap          = 50;                   % overlap between sliding windows (in %)
            opt.WinAverage          = 5;                    % number of sliding windows averaged by time point
            
            % specparam opts
            opt.freq_range          = [1, 40];
            opt.peak_width_limits   = [2, 6];
            opt.max_peaks           = 3;
            opt.min_peak_height     = 0.5;
            opt.aperiodic_mode      = 'fixed';      % alternative: knee
            opt.peak_threshold      = 2.0;          % 2 std dev: parameter for interface simplification
            
            % Matlab-only options
            opt.peak_type           = 'gaussian';   % alternative: cauchy
            opt.proximity_threshold = 1;            % to allow fitting peaks closer to the edge of the spectrum, default = 2
            opt.guess_weight        = 'none';
            opt.thresh_after        = true;         % threshold after fitting, always selected for Matlab
            opt.hOT                 = 1;            % using optimization toolbox
            opt.rmoutliers          = 'yes';
            opt.maxfreq             = 2.5;
            opt.maxtime             = 6;
            opt.minnear             = 3;
            
            % define variables
            freqs                   = 0 : (1 / opt.WinLength) : (opt.sfreq / 2);

            % SPRiNT analysis
            try
                % compute short-time Fourier transform
                [TF, ts]            = SPRiNT_stft(data.trial{:}, opt);
                sprintData          = SPRiNT_specparam_matlab(TF, freqs, opt, ts); % parameterize STFTs
                
                % save structure
                mkdir(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name));
                TG_save(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, saveName), sprintData);
            catch
                continue;
            end
        end
    end
end

% shut down parallel pool
delete(gcp);
