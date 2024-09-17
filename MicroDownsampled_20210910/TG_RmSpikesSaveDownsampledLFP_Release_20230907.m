%==========================================================================
% This script removes spikes and then downsamples and saves the LFP signals
% in the fieldtrip structure.
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.micro              	= 'D:\TreasureHunt\Micro_20210111'; % data path
paths.SUData              	= 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.save                  = 'D:\TreasureHunt\MicroDownsampled_20210910'; % save folder  

% add functions and toolboxes
addpath(genpath('D:\External\Functions'));
addpath('D:\External\Toolboxes\Wave_clus\Batch_files'); % for bandpass filter
addpath('D:\External\Toolboxes\Fieldtrip');

% default fieldtrip settings
ft_defaults;

% set configurations for bandpass-filter
bpcfg                       = [];
bpcfg.detect_order          = 4; % filter oder
bpcfg.detect_fmin           = 100; % minimum frequency
bpcfg.detect_fmax           = 3000; % maximum frequency

% set configurations for spike removal
rmcfg                       = [];
rmcfg.msBefore              = 1; % distance in ms to be included in mean spike
rmcfg.msAfter               = 2; % as in Zanos et al., 2011, JNeurophysiol

% set configurations for ft_resampledata
resampcfg                   = [];
resampcfg.resamplefs        = 2000;
resampcfg.detrend           = 'no';
resampcfg.demean            = 'no';
resampcfg.feedback          = 'no';

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
        
        % load sampling rates for all wires (for parallel computing)
        wireSr = nan(size(chanDir));
        for iSr = 1:size(chanDir, 1)
            thisSr         = load(fullfile(paths.micro, subjects{iSub}, sessions(iSess).name, chanDir(iSr).name, 'datacut.mat'), 'sr');
            wireSr(iSr, 1) = thisSr.sr;
        end
        
        % define sampling rate
        if size(unique(wireSr), 1) == 1
            bpcfg.sr = mode(wireSr); % sampling rate
        else
            error(strcat('sampling rates differ between wires for', 32, subjects{iSub}, 32, sessions(iSess).name));
        end

        % loop through channels
        parfor (iWire = 1:size(chanDir, 1), 16)
            
            % print channel name
            fprintf('\tChannel name: %s.\n', chanDir(iWire).name);
            
            % load raw data
            data            = load(fullfile(paths.micro, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, 'datacut.mat'));
            
            % band-pass filter data
            bpData          = spike_detection_filter(data.data, bpcfg);

            % copy data for spike removal
            rmData          = data.data;
            
            %% spike removal

            % path for this channel's unit data
            mwChanPath      = fullfile(paths.SUData, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name);

            % try spike removal
            try
                % load wave-clus output (for spike-depiction etc.)
                t = load(fullfile(mwChanPath, 'times_datacut.mat'));

                % get different clusters
                clusters        = unique(t.cluster_class(:, 1));

                % remove cluster with unassigned and artifact spikes
                clusters        = clusters(clusters > 0);

                % remove spikes, if clusters exist
                if ~isempty(clusters)

                    % spike peak times
                    peakTimes       = t.cluster_class(:, 2);

                    % spike peak samples
                    peakSamples     = round(peakTimes / 1000 * data.sr);

                    % loop through different clusters
                    for iClus = 1:size(clusters, 1)

                        % get spikes of this clus
                        clusIdx             = (t.cluster_class(:, 1) == clusters(iClus));

                        % get peak samples
                        clusPeakSamples     = peakSamples(clusIdx);

                        % samples before and after peak
                        addSamples          = -rmcfg.msBefore / 1000 * data.sr : rmcfg.msAfter / 1000 * data.sr;

                        % get all cluster spike samples
                        allSamples          = clusPeakSamples + addSamples;

                        % calculate mean spike
                        meanSpike           = mean(bpData(allSamples));

                        % tapering
                        meanSpike           = meanSpike - mean([meanSpike(1); meanSpike(end)]);
                        meanSpike           = TG_taper(meanSpike, 1/8);

                        % replace samples
                        repSamples          = repmat(meanSpike, sum(clusIdx), 1);

                        % remove spikes
                        rmData(allSamples)  = rmData(allSamples) - repSamples;
                    end
                end
            catch
                fprintf('No spikes removed for %s, session %d, channel %d', subjects{iSub}, iSess, iWire);
            end

            %% downsample data
             
            % convert data into fieldtrip format
            ftData                 = [];
            ftData.fsample         = data.sr;
            ftData.label           = {chanDir(iWire).name};
            ftData.trial           = {rmData};
            ftData.time            = {(1:size(rmData, 2)) / data.sr};
            
            % downsample data
            downsampled_data       = ft_resampledata(resampcfg, ftData);
            
            % save downsampled data
            mkdir(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name));
            TG_save(fullfile(paths.save, subjects{iSub}, sessions(iSess).name, chanDir(iWire).name, 'datacutFt2000Hz.mat'), downsampled_data);
        end
    end
end

% shut down parallel pool
delete(gcp);
