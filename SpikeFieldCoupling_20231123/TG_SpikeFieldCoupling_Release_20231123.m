%==========================================================================
% This script calculates the spike-triggered average and spike-field
% coherence for each unit
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% set random seed
randSeedNum         = 444;
rng(randSeedNum, 'twister');
randSeedNum         = randi(100000, 100000, 1); % for randomseed

% paths
paths.behlog        = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.SUData        = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.LFPData       = 'D:\TreasureHunt\MicroDownsampled_20210910'; % LFP data
paths.artifact      = 'D:\TreasureHunt\ArtifactDetection_20230907'; % IED/artifact detection
paths.subInfo       = 'D:\TreasureHunt\SubjectInformation_20210111'; % subject information
paths.phaseRes      = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % results from phase analysis
paths.save          = 'D:\TreasureHunt\SpikeFieldCoupling_20231123'; % save folder
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% add functions
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

%% set configurations

% parameters
param               = [];
param.c4aName       = {'Cluster4Analysis_LK_20200416_113912', 'Cluster4Analysis_TG_20210315_020120'};
param.timeBorder    = 0.5; % time before and after spike for spike-field coherence in seconds
param.sfcFs         = 250; % sampling rate for SFC
param.nSubSmp       = 100; % number of subsamples

% settings for spectral analysis
specCfg             = [];
specCfg.output      = 'pow';
specCfg.channel     = 'all';
specCfg.method      = 'mtmfft';
specCfg.taper       = 'dpss'; % multitaper
specCfg.tapsmofrq   = 4; % fw: width of frequency smoothing in Hz
specCfg.foi         = 1:1:100; % frequencies of interest

% configurations for random locations in the Treasure Hunt arena to normalize drop error
RandLoccfg          = [];
RandLoccfg.maxR     = 50;
RandLoccfg.minR     = 0;
RandLoccfg.N        = 1000000;
RandLoccfg.centerX  = 370;
RandLoccfg.centerY  = 360;
randLocs            = TG_RandomPointsInCircle(RandLoccfg);

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

%% preallocate
allRes         = cell(length(subjects), 1);

% for sanity checking
allExByWC      = cell(length(subjects), 1); % excluded by wave-clus
allExByVisInsp = cell(length(subjects), 1); % excluded due to visual inspection

%% loop through subjects
parfor (iSub = 1:length(subjects), 9)
% for iSub = 1:length(subjects)

    % get sessions
    sessions         = dir(fullfile(paths.SUData, subjects{iSub}, 'session*'));

    % preallocate
    subRes           = [];

    % set random number generator
    rng(randSeedNum(iSub, 1), 'twister');

    % for sanity checking
    exByWC           = []; % excluded by wave-clus
    exByVisInsp      = []; % excluded due to visual inspection

    %% get each channel's brain region

    % load information about brain region
    micro2regionFile = dir(fullfile(paths.subInfo, subjects{iSub}, 'Microwires', 'Micro2Region.txt'));
    fileID           = fopen(fullfile(micro2regionFile.folder, micro2regionFile.name));
    micro2Region     = textscan(fileID, '%s %s');
    fclose(fileID);

    % load information about left or right hemisphere
    micro2macroFile  = dir(fullfile(paths.subInfo, subjects{iSub}, 'Microwires', 'Micro2Macro.txt'));
    fileID           = fopen(fullfile(micro2macroFile.folder, micro2macroFile.name));
    micro2macro      = textscan(fileID, '%*s %s');
    fclose(fileID);

    % information about hemisphere
    hemisphereShort  = cellfun(@(x) x(end), micro2macro{:}, 'UniformOutput', false);
    hemisphere       = replace(hemisphereShort, {'R'; 'L'; 'd'}, {'right'; 'left'; 'NotImplanted'});

    % channel list
    chanNumbers      = cellfun(@str2num, micro2Region{:, 1}, 'UniformOutput', false);
    chanList         = cat(2, chanNumbers{:, 1});

    % create list of each channel's brain region
    brainRegionsList       = cell(size(chanList, 2), 2);
    brainRegionsList(:, 2) = num2cell(chanList');
    regionNames            = [];
    for iMacro = 1:size(chanNumbers, 1)
        regionNames = cat(1, regionNames, repmat(strcat(micro2Region{1, 2}(iMacro, 1), '_', hemisphere{iMacro, 1}), numel(chanNumbers{iMacro, 1}), 1));
    end
    brainRegionsList(:, 1) = regionNames;

    %% loop through sessions
    for iSess = 1:size(sessions, 1)

        % define variables
        encSegmentsOrig = [];
        bGoodMem        = [];
        encSpikeIdx     = [];
        smpEncSucc      = [];
        smpEncFail      = [];
        smpRecSucc      = [];
        smpRecFail      = [];
        encSuccSFC      = [];
        encFailSFC      = [];
        recSuccSFC      = [];
        recFailSFC      = [];
        encSuccSTP      = [];

        % original session index
        sessIdx         = split(sessions(iSess).name, '_');
        sessIdx         = str2double(sessIdx{2});

        % display session information
        fprintf('\n==================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);

        % get available data
        mwChanDir       = TG_GetChanDir_20210812(paths.SUData, subjects{iSub}, sessions(iSess).name);  % unit data

        % brain region sanity check
        if ~isequal(size(mwChanDir, 1), size(brainRegionsList, 1))
            error('different number of microwire channels than channels in the brain region file');
        end

        %% chest-wise memory performance

        % load behavioral logfile and segment info file
        trialInfo       = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo       = trialInfo.trialInfo;
        segmentInfo     = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'segmentInfo_20230208.mat'));

        % each chest's Treasure Hunt trial index
        encTrialIdx     = transpose(rude(trialInfo.numChestsPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numChestsPerTrial))));
        objTrialIdx     = transpose(rude(trialInfo.numObjRecallPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numObjRecallPerTrial))));
        locTrialIdx     = transpose(rude(trialInfo.numLocRecallPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numLocRecallPerTrial))));

        % index for encoding encSegments
        encSegmentIdx   = trialInfo.chestIdx;

        % prelloacte - object recall
        encObjRec       = cell(size(encTrialIdx, 1), 1);
        objSegmentIdx   = nan(size(trialInfo.audioResponse, 1), 1);

        % prelloacte - location recall
        encLocRec       = nan(size(encTrialIdx, 1), 1);
        locRec          = nan(size(locTrialIdx, 1), 1);
        locSegmentIdx   = nan(size(locTrialIdx, 1), 1);

        % loop through chests
        for iChest = 1:size(encTrialIdx, 1)

            % object recall
            if strcmp(trialInfo.recallType{encTrialIdx(iChest)}, 'OBJECT')
                thisChestLoc = trialInfo.chestLoc(iChest, :);
                idxObjRec    = abs((trialInfo.cueingLoc(:, 1) - thisChestLoc(1))) < 0.0001 ...
                    & abs(trialInfo.cueingLoc(:, 3) - thisChestLoc(3)) < 0.0001 ...
                    & encTrialIdx(iChest) == objTrialIdx;
                encObjRec(iChest, 1)        = trialInfo.audioResponse(idxObjRec);
                objSegmentIdx(idxObjRec, 1) = iChest;

                % location recall
            elseif strcmp(trialInfo.recallType{encTrialIdx(iChest)}, 'LOCATION')
                thisChestLabel = trialInfo.TREASURE_LABEL_EN(iChest, 1);
                idxLocRec      = strcmp(trialInfo.cueingObject, thisChestLabel) & encTrialIdx(iChest) == locTrialIdx;
                if sum(idxLocRec == 1)
                    D                           = pdist2(trialInfo.CORRECT_TEST_POSITION(idxLocRec, [1, 3]), trialInfo.CHOSEN_TEST_POSITION(idxLocRec, [1, 3]));
                    DSurro                      = pdist2(trialInfo.CORRECT_TEST_POSITION(idxLocRec, [1, 3]), randLocs);
                    encLocRec(iChest, 1)        = sum(D < DSurro) / numel(DSurro);
                    locRec(idxLocRec, 1)        = sum(D < DSurro) / numel(DSurro);
                    locSegmentIdx(idxLocRec, 1) = iChest;
                else
                    warning('No location recall for %s, session %d, chest %d', subjects{iSub}, iSess, iChest);
                end
            end
        end

        % memory performance
        bGoodMemEnc     = strcmp(encObjRec, 'CORRECT') | encLocRec > median(encLocRec, 'omitnan'); % or: > 0.9; % logical index for good memory performance during encoding
        bGoodMemEncAbs  = strcmp(encObjRec, 'CORRECT') | encLocRec > 0.9;
        bGoodMemObj     = strcmp(trialInfo.audioResponse, 'CORRECT'); % logical index for successful object recall
        bGoodMemLoc     = locRec > median(locRec, 'omitnan'); % or '> 0.9'; % logical index for good location recall performance
        bGoodMemLocAbs  = locRec > 0.9;

        % concatenate object and location recall
        bGoodMemRec     = cat(1, bGoodMemObj, bGoodMemLoc);
        recSegmentIdx   = cat(1, objSegmentIdx, locSegmentIdx);

        % classify object (1) and location recall (2) encSegments
        encObjOrLoc                = zeros(size(encSegmentIdx, 1), 1);
        encObjOrLoc(objSegmentIdx(~isnan(objSegmentIdx))) = 1;
        encObjOrLoc(locSegmentIdx) = 2;
        recObjOrLoc                = [ones(size(bGoodMemObj, 1), 1); ones(size(bGoodMemLoc, 1), 1) * 2];

        %% loop through channels
        for iWire = 1:size(mwChanDir, 1)

            % print channel name
            fprintf('\tChannel name: %s.\n', mwChanDir(iWire).name);

            % path for this channel's unit data
            mwChanPath          = fullfile(paths.SUData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name);

            % load LFP data
            LFPdata             = load(fullfile(paths.LFPData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'datacutFt2000Hz.mat'));

            % load data for artifact removal
            artifactData        = load(fullfile(paths.artifact, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'detectedArtifacts.mat'), 'bArtifact');

            % load wave-clus output (for spike-depiction etc.)
            try
                t = load(fullfile(mwChanPath, 'times_datacut.mat'));
            catch
                exByWC = cat(1, exByWC, [iSub, iSess, iWire]); % bookkeeping
                fprintf('No wave_clus for %s, session %d, channel %d', subjects{iSub}, iSess, iWire);
                continue;
            end

            % load unit classification file
            try
                c4a = load(fullfile(mwChanPath, param.c4aName{1, 1}));
            catch
                c4a = load(fullfile(mwChanPath, param.c4aName{1, 2}));
            end
            c4a = c4a.Cluster4Analysis;

            %% loop through each cluster
            for iClus = 1:max(t.cluster_class(:, 1))

                % print cluster number
                fprintf('\t\tCluster: %d.\n', iClus);

                % continue if you decided that this cluster is not nice
                if strcmp(c4a(cell2mat({c4a.Cluster}) == iClus).Decision, 'no')
                    exByVisInsp = cat(1, exByVisInsp, [iSub, iSess, iWire, iClus]);
                    fprintf('You decided not to analyze this cluster.\n');
                    continue;
                end

                % get cluster
                thisClusterOrig         = t.cluster_class(t.cluster_class(:, 1) == iClus, :); % cluster-number, time in ms
                thisClusSamplesOrig     = round(thisClusterOrig(:, 2) * LFPdata.fsample / 1000); % convert from ms to LFP samples

                % IEDs/artifact removal
                thisClusArtifactIdx     = artifactData.bArtifact(thisClusSamplesOrig);
                thisClusSamples         = thisClusSamplesOrig(~thisClusArtifactIdx);

                %% spike-related information

                % adjust sampling rate of segment info file to LFP data
                encSegments     = segmentInfo.encoding / (segmentInfo.fsample / LFPdata.fsample);
                recSegments     = cat(1, segmentInfo.objRecall, segmentInfo.locRecall) / (segmentInfo.fsample / LFPdata.fsample);

                % loop through encoding segments
                encSpikeIdx     = cell(size(encSegments, 1), 1);
                encSpikeSamples = cell(size(encSegments, 1), 1);
                for iSeg = 1:size(encSegments, 1)

                    % get information
                    encSpikeIdx{iSeg, 1}      = find(thisClusSamples >= encSegments(iSeg, 1) & thisClusSamples <= encSegments(iSeg, 2));
                    encSpikeSamples{iSeg, 1}  = thisClusSamples(encSpikeIdx{iSeg, 1});
                end

                % loop through recall segments
                recSpikeIdx     = cell(size(recSegments, 1), 1);
                recSpikeSamples = cell(size(recSegments, 1), 1);
                for iSeg = 1:size(recSegments, 1)

                    % get information
                    recSpikeIdx{iSeg, 1}      = find(thisClusSamples >= recSegments(iSeg, 1) & thisClusSamples <= recSegments(iSeg, 2));
                    recSpikeSamples{iSeg, 1}  = thisClusSamples(recSpikeIdx{iSeg, 1});
                end

                % number of spikes during encoding and recall
                encNumSpikes        = cell2mat(cellfun(@(x) size(x, 1), encSpikeIdx, 'UniformOutput', 0));
                recNumSpikes        = cell2mat(cellfun(@(x) size(x, 1), recSpikeIdx, 'UniformOutput', 0));

                % memory performance index per spike
                bGoodMemEncSpikes   = repelem(bGoodMemEnc, encNumSpikes);
                bGoodMemRecSpikes   = repelem(bGoodMemRec, recNumSpikes);

                %% spike-field coherence

                % spike samples - encoding
                allEncSpikeIdx      = cat(1, encSpikeIdx{:});
                allEncSpikeSmp      = thisClusSamples(allEncSpikeIdx);

                % spike samples - retrieval
                allRecSpikeIdx      = cat(1, recSpikeIdx{:});
                allRecSpikeSmp      = thisClusSamples(allRecSpikeIdx);

                % continue, if there are no spikes during encoding or
                % retrieval (unit would be excluded later anyways)
                if isempty(allEncSpikeIdx) || isempty(allRecSpikeIdx)
                    unitRes = [];
                    subRes  = cat(1, subRes, unitRes);
                    continue;
                end

                % offset of first spike sample point
                smpOffsetST         = (param.timeBorder * LFPdata.fsample);

                % start and end samples - encoding
                encStartSmpST       = allEncSpikeSmp - smpOffsetST;
                encEndSmpST         = allEncSpikeSmp + smpOffsetST;

                % start and end samples - recall
                recStartSmpST       = allRecSpikeSmp - smpOffsetST;
                recEndSmpST         = allRecSpikeSmp + smpOffsetST;

                % extract spike-triggered LFP - encoding
                encCfg              = [];
                encCfg.trl          = [encStartSmpST, encEndSmpST, repmat(-smpOffsetST, size(allEncSpikeSmp))];
                encST               = ft_redefinetrial(encCfg, LFPdata);

                % extract spike-triggered LFP - recall
                recCfg              = [];
                recCfg.trl          = [recStartSmpST, recEndSmpST, repmat(-smpOffsetST, size(allRecSpikeSmp))];
                recST               = ft_redefinetrial(recCfg, LFPdata);

                % downsample
                rsmpCfg             = [];
                rsmpCfg.resamplefs  = param.sfcFs;
                encST               = ft_resampledata(rsmpCfg, encST);
                recST               = ft_resampledata(rsmpCfg, recST);

                % divide into successful and unsuccessful spikes - encoding
                encSuccST           = encST;
                encSuccST.trial     = encSuccST.trial(1, bGoodMemEncSpikes);
                encSuccST.time      = encSuccST.time(1, bGoodMemEncSpikes);
                encFailST           = encST;
                encFailST.trial     = encFailST.trial(1, ~bGoodMemEncSpikes);
                encFailST.time      = encFailST.time(1, ~bGoodMemEncSpikes);

                % divide into successful and unsuccessful spikes - recall
                recSuccST           = recST;
                recSuccST.trial     = recSuccST.trial(1, bGoodMemRecSpikes);
                recSuccST.time      = recSuccST.time(1, bGoodMemRecSpikes);
                recFailST           = recST;
                recFailST.trial     = recFailST.trial(1, ~bGoodMemRecSpikes);
                recFailST.time      = recFailST.time(1, ~bGoodMemRecSpikes);

                % number of spikes - encoding
                numEncSucc          = size(encSuccST.trial, 2);
                numEncFail          = size(encFailST.trial, 2);

                % number of spikes - recall
                numRecSucc          = size(recSuccST.trial, 2);
                numRecFail          = size(recFailST.trial, 2);

                % continue, if there are no spikes
                if numEncSucc == 0 || numEncFail == 0 || numRecSucc == 0 || numRecFail == 0
                    unitRes = [];
                    subRes  = cat(1, subRes, unitRes);
                    continue;
                end

                %% subsampling
                encSuccSFC          = nan(param.nSubSmp, size(specCfg.foi, 2));
                encFailSFC          = nan(param.nSubSmp, size(specCfg.foi, 2));
                recSuccSFC          = nan(param.nSubSmp, size(specCfg.foi, 2));
                recFailSFC          = nan(param.nSubSmp, size(specCfg.foi, 2));
                for iSmp = 1:param.nSubSmp

                    %% encoding
                    if numEncSucc > numEncFail

                        % draw a subsample from the successful spikes
                        smpEncFail          = encFailST;
                        smpEncSucc          = encSuccST;
                        smpEncSucc.trial    = datasample(encSuccST.trial, size(encFailST.trial, 2), 'Replace', false);
                        smpEncSucc.time     = encFailST.time;
                    elseif numEncSucc < numEncFail

                        % draw a subsample from the unsuccessful spikes
                        smpEncSucc          = encSuccST;
                        smpEncFail          = encFailST;
                        smpEncFail.trial    = datasample(encFailST.trial, size(encSuccST.trial, 2), 'Replace', false);
                        smpEncFail.time     = encSuccST.time;
                    elseif numEncSucc == numEncFail

                        % no subsampling needed
                        smpEncSucc          = encSuccST;
                        smpEncFail          = encFailST;
                    end

                    %% recall
                    if numRecSucc > numRecFail

                        % draw a subsample from the successful spikes
                        smpRecFail          = recFailST;
                        smpRecSucc          = recSuccST;
                        smpRecSucc.trial    = datasample(recSuccST.trial, size(recFailST.trial, 2), 'Replace', false);
                        smpRecSucc.time     = recFailST.time;
                    elseif numRecSucc < numRecFail

                        % draw a subsample from the unsuccessful spikes
                        smpRecSucc          = recSuccST;
                        smpRecFail          = recFailST;
                        smpRecFail.trial    = datasample(recFailST.trial, size(recSuccST.trial, 2), 'Replace', false);
                        smpRecFail.time     = recSuccST.time;

                    elseif numRecSucc == numRecFail

                        % no subsampling needed
                        smpRecSucc          = recSuccST;
                        smpRecFail          = recFailST;
                    end

                    %%  SFC

                    % spike-triggered average - encoding
                    avgCfg              = [];
                    encSuccSTA          = ft_timelockanalysis(avgCfg, smpEncSucc);
                    encFailSTA          = ft_timelockanalysis(avgCfg, smpEncFail);

                    % spike-triggered average - recall
                    avgCfg              = [];
                    recSuccSTA          = ft_timelockanalysis(avgCfg, smpRecSucc);
                    recFailSTA          = ft_timelockanalysis(avgCfg, smpRecFail);

                    % spike-triggered power - encoding
                    encSuccSTP          = ft_freqanalysis(specCfg, smpEncSucc);
                    encFailSTP          = ft_freqanalysis(specCfg, smpEncFail);

                    % spike-triggered power - recall
                    recSuccSTP          = ft_freqanalysis(specCfg, smpRecSucc);
                    recFailSTP          = ft_freqanalysis(specCfg, smpRecFail);

                    % power of spike-triggered average - encoding
                    encSuccSTApow       = ft_freqanalysis(specCfg, encSuccSTA);
                    encFailSTApow       = ft_freqanalysis(specCfg, encFailSTA);

                    % power of spike-triggered average - recall
                    recSuccSTApow       = ft_freqanalysis(specCfg, recSuccSTA);
                    recFailSTApow       = ft_freqanalysis(specCfg, recFailSTA);

                    % calculate spike-field coherence (SFC) - encoding
                    encSuccSFC(iSmp, :) = encSuccSTApow.powspctrm ./ encSuccSTP.powspctrm;
                    encFailSFC(iSmp, :) = encFailSTApow.powspctrm ./ encFailSTP.powspctrm;

                    % calculate spike-field coherrece (SFC) - recall
                    recSuccSFC(iSmp, :) = recSuccSTApow.powspctrm ./ recSuccSTP.powspctrm;
                    recFailSFC(iSmp, :) = recFailSTApow.powspctrm ./ recFailSTP.powspctrm;
                end

                %% collect information for this unit
                
                % basics
                unitRes                  = [];
                unitRes.idx              = [iSub, sessIdx, iWire, iClus];
                unitRes.brainRegion      = brainRegionsList{iWire, 1};
                unitRes.foi              = encSuccSTP.freq;

                % encoding
                unitRes.bGoodMemEnc      = bGoodMemEnc;
                unitRes.encSuccSFC       = mean(encSuccSFC);
                unitRes.encFailSFC       = mean(encFailSFC);

                % recall
                unitRes.bGoodMemRec      = bGoodMemRec;
                unitRes.recSuccSFC       = mean(recSuccSFC);
                unitRes.recFailSFC       = mean(recFailSFC);

                % collapse across units
                subRes = cat(1, subRes, unitRes);
            end
        end
    end

    % collect results across subjects
    allRes{iSub, 1}         = subRes;
    allExByWC{iSub, 1}      = exByWC;
    allExByVisInsp{iSub, 1} = exByVisInsp;
end

% unnest result cells
allRes                  = cat(1, allRes{:});
allExByWC               = cat(1, allExByWC{:});
allExByVisInsp          = cat(1, allExByVisInsp{:});
allIdx                  = cat(1, allRes.idx);

% load phase analysis results for exlcusion index
phaseRes                = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));
phaseResIdx             = cat(1, phaseRes.phaseRes.idx);
incIdx                  = ismember(allIdx, phaseResIdx, 'rows');

% save important output
save(fullfile(paths.save, 'spikeFieldResults'), 'allRes', 'allExByWC', 'allExByVisInsp', 'allIdx', '-v7.3');

% include only cells which are also included in phase analysis
incRes                  = allRes(incIdx);

% sanity check
if size(incRes, 1) ~= size(phaseRes.phaseRes, 1)
    error('different structure sizes');
end

% collect SFC results - encoding
encSucc                 = cat(1, incRes.encSuccSFC);
encFail                 = cat(1, incRes.encFailSFC);

% t-tests
[~, encP, ~, encT]      = ttest(encSucc, encFail);
encP                    = encP * size(encSucc, 2) * 2; % Bonferroni correction
encH                    = any(encP < 0.05);

% collect SFC results - retrieval
recSucc                 = cat(1, incRes.recSuccSFC);
recFail                 = cat(1, incRes.recFailSFC);

% t-tests
[~, recP, ~, recT]      = ttest(recSucc, recFail);
recP                    = recP * size(recSucc, 2) * 2; % Bonferroni correction
recH                    = any(recP < 0.05);

% plot - encoding
encFigSFC = figure;
TG_ShadeSEM_20210714(incRes(1).foi, encSucc * 100, [0.1, 0.6, 0.1], 0.4);
hold on;
TG_ShadeSEM_20210714(incRes(1).foi, encFail * 100, [0.6, 0.1, 0.1], 0.4);
xlabel('Frequency (Hz)');
xlim([0, 100]);
set(gca, 'XScale', 'log', 'tickDir', 'out');
box off;
ylabel('SFC (%)');
ylim([0, 4]);
title('Encoding');

% print figure
print(encFigSFC, fullfile(paths.save, 'encodingSFC'), '-dsvg', '-r300');

% plot - retrieval
recFigSFC = figure;
TG_ShadeSEM_20210714(incRes(1).foi, recSucc * 100, [0.1, 0.6, 0.1], 0.4);
hold on;
TG_ShadeSEM_20210714(incRes(1).foi, recFail * 100, [0.6, 0.1, 0.1], 0.4);
xlabel('Frequency (Hz)');
xlim([0, 100]);
set(gca, 'XScale', 'log', 'tickDir', 'out');
box off;
ylabel('SFC (%)');
ylim([0, 3]);
title('Retrieval');

% print figure
print(recFigSFC, fullfile(paths.save, 'recallSFC'), '-dsvg', '-r300');
