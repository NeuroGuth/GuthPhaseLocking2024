%==========================================================================
% This script aligns the times from the behavioral logfile 'trialInfo.mat'
% with the corresponding sample points of the LFP microwire signal. It
% detects the start and end sampling points of encoding, object recall
% and location recall.
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.data                 = 'D:\TreasureHunt\DataComplete_20200521\';
paths.micro                = 'D:\TreasureHunt\Micro_20210111\';
paths.behlog               = 'D:\TreasureHunt\Beh_20210111\';

% set parameters
params                     = [];
params.encDuration         = 1.5; % encoding segment duration in s
params.objDuration         = 4; % object recall segment duration in s
params.mwTrigTreshold      = 30000; % trigger threshold is 30,000
params.mwTreshTH09         = 5000;  % for this subject, trigger threshold is 5,000
params.condition           = {'encoding'; 'objectRecall'; 'locationRecall'};

% subjects
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

%% loop through subjects
parfor iSub = 1:length(subjects)
    
    % get sessions
    sessions = dir(fullfile(paths.data, subjects{iSub}, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % display session information
        fprintf('\n===================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        %% define trial

        % process behavioral trigger information
        EEGLogfile                    = dir(fullfile(sessions(iSess).folder, sessions(iSess).name, 'Beh', '*EEGLog.txt'));
        fileID                        = fopen(fullfile(EEGLogfile.folder, EEGLogfile.name));
        behtrig                       = textscan(fileID, '%f\t%d\t%s');
        fclose(fileID);
        behTrigTimepoints             = behtrig{1}(strcmp(behtrig{3}, 'ON')) ./ 1000; % [sec]
        
        % process microwire trigger information
        trSignal                      = load(fullfile(paths.micro, subjects{iSub}, sessions(iSess).name, 'ainp2', 'datacut.mat'));
        mwSamplingRate                = trSignal.sr;
        mwTrigData                    = trSignal.data;
        mwTrigSamples                 = 1:length(mwTrigData); % all sampling points
        mwTrigSamplepoints            = mwTrigSamples(diff(mwTrigData < params.mwTrigTreshold) > 0); % identify microwire sampling points of trigger start (when the photodiode was switched on)
        if strcmp(subjects{iSub, 1}, 'TH09')
            mwTrigSamplepoints        = mwTrigSamples(diff(mwTrigData > params.mwTreshTH09) > 0);
        end
        
        % sanity check
        [rho, pval]         = corr(diff(behTrigTimepoints), diff(mwTrigSamplepoints ./ mwSamplingRate)', 'type', 'spearman');
        if rho < 0.999
            error('Correlation between inter-trigger intervals not sufficient.');
        end
        
        % load behavioral logfile
        trialInfo           = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo           = trialInfo.trialInfo;
        
        % preallocate results structure
        segmentInfo         = [];
        segmentInfo.fsample = mwSamplingRate;
        
        % loop through encoding, object recall and location recall
        startTime           = [];
        for iCond = 1:size(params.condition, 1)

            %% get behavioral start times

            % encoding
            if strcmp(params.condition{iCond}, 'encoding')
                startTime = transpose(trialInfo.PLAYER_CHEST_ROTATION_ENDED(:, 1)); % get behavioral time points of encoding start

            % object recall
            elseif strcmp(params.condition{iCond}, 'objectRecall')
                startTime = transpose(trialInfo.RECORDING_STARTED(:, 1)); % get behavioral time points of object recall start

            % location recall
            elseif strcmp(params.condition{iCond}, 'locationRecall')
                startTime = transpose(trialInfo.timeCue4LocRecall(:, 1));  % get behavioral time points of location recall start
            end

            %% identify microwire start samples

            % find microwire sample points during segment starts by identifying the two trigger points closest to each segment start
            trDistancesStart                        = behTrigTimepoints - startTime;
            trDistancesStart(trDistancesStart >= 0) = nan;
            [minValuesStart, closestIndexStart]     = max(trDistancesStart); % index of triggers directly before events start and durations from triggers to corresponding events
            leftTrigValuesStart                     = transpose(behTrigTimepoints(closestIndexStart)); % behavioral times of triggers directly before events
            rightTrigValuesStart                    = transpose(behTrigTimepoints(closestIndexStart + 1)); % behavioral times of triggers directly after events
            ratioBehTrigStart                       = abs(minValuesStart) ./ (rightTrigValuesStart - leftTrigValuesStart); % each event's location between preceding and subsequent trigger
            
            % identification of corresponding microwire sampling points in microwires using 'ratioBehTrig'
            mwLeftTrigSamplesStart                  = mwTrigSamplepoints(closestIndexStart);
            mwRightTrigSamplesStart                 = mwTrigSamplepoints(closestIndexStart + 1);
            mwSampleDifferenceStart                 = mwRightTrigSamplesStart - mwLeftTrigSamplesStart;

            % get start samples
            startMwSamples                          = round(mwLeftTrigSamplesStart + (ratioBehTrigStart .* mwSampleDifferenceStart));

            %% create results structure

            % encoding
            if strcmp(params.condition{iCond}, 'encoding')
                
                % create encoding start and end samples
                encSampleDur            = params.encDuration * mwSamplingRate;
                segmentInfo.encoding    = [startMwSamples', startMwSamples' + encSampleDur];

            % object recall
            elseif strcmp(params.condition{iCond}, 'objectRecall')
                
                % create object recall start and end samples
                objSampleDur            = params.objDuration * mwSamplingRate;
                segmentInfo.objRecall   = [startMwSamples', startMwSamples' + objSampleDur];

            % location recall
            elseif strcmp(params.condition{iCond}, 'locationRecall')

                % end time
                endTime                             = transpose(trialInfo.timeLocRecall(:, 1));

                % find microwire sample points during segment ends by identifying the two trigger points closest to each segment end
                trDistancesEnd                      = behTrigTimepoints - endTime;
                trDistancesEnd(trDistancesEnd >= 0) = nan;
                [minValuesEnd, closestIndexEnd]     = max(trDistancesEnd); % index of triggers directly before events end and durations from triggers to corresponding events
                leftTrigValuesEnd                   = transpose(behTrigTimepoints(closestIndexEnd)); % behavioral times of triggers directly before events
                rightTrigValuesEnd                  = transpose(behTrigTimepoints(closestIndexEnd + 1)); % behavioral times of triggers direclty after events
                ratioBehTrigEnd                     = abs(minValuesEnd) ./ (rightTrigValuesEnd - leftTrigValuesEnd); % each event's location between preceding and subsequent trigger

                % identification of corresponding microwire sampling points in microwires using 'ratioBehTrig'
                mwLeftTrigSamplesEnd                = mwTrigSamplepoints(closestIndexEnd);
                mwRightTrigSamplesEnd               = mwTrigSamplepoints(closestIndexEnd + 1);
                mwSampleDifferenceEnd               = mwRightTrigSamplesEnd - mwLeftTrigSamplesEnd;

                % get end samples
                endMwSamples                        = round(mwLeftTrigSamplesEnd + (ratioBehTrigEnd .* mwSampleDifferenceEnd));

                % segment info
                segmentInfo.locRecall               = [startMwSamples', endMwSamples'];
            end
        end
        
        % save results
        TG_save(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, strcat('segmentInfo_20230208.mat')), segmentInfo);
    end
end

