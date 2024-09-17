%==========================================================================
% This script performs power analysis for LFP data. It compares power
% during successful and unsuccessful encoding and recall segments
% using cluster-based permutation testing.
%
% Tim Guth, 2024
%==========================================================================

%% settings
clc; close all; clear;
rng(444);

% paths
paths.behlog    = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.SUData    = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.LFPData   = 'D:\TreasureHunt\MicroDownsampled_20210910'; % LFP data
paths.subInfo   = 'D:\TreasureHunt\SubjectInformation_20210111'; % subject information
paths.save      = 'D:\TreasureHunt\FreqanalysisSegments_20220110\Freqanalysis_20231130'; % save folder
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% own functions
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

%% set configurations

% parameters
param.freqMin        = 1; % minimum frequency to analyze
param.freqMax        = 40; % maximum frequency to analyze
param.freqNum        = 40; % number of frequencies of interest
param.buffer         = 5; % borders before and after segment in seconds
param.segDur         = 4; % duration of each segment to cut in seconds
param.timeRes        = 0.01; % length of time bins for TFR analysis
param.condition      = {'encoding'; 'objectRecall'; 'locationRecall'};
param.c4aName        = {'Cluster4Analysis_LK_20200416_113912', 'Cluster4Analysis_TG_20210315_020120'};

% for random locations in the Treasure Hunt arena to normalize drop error
RandLoccfg           = [];
RandLoccfg.maxR      = 50;
RandLoccfg.minR      = 0;
RandLoccfg.N         = 1000000;
RandLoccfg.centerX   = 370;
RandLoccfg.centerY   = 360;
randLocs             = TG_RandomPointsInCircle(RandLoccfg);

% for ft_freqanalysis
TFRcfg               = [];
TFRcfg.method        = 'wavelet';
TFRcfg.pad           = 'nextpow2';
TFRcfg.output        = 'pow';
TFRcfg.keeptrials    = 'yes';
TFRcfg.width         = 5; % number of wavelet cycles
TFRcfg.foi           = round(logspace(log10(param.freqMin), log10(param.freqMax), param.freqNum), 2); % see Miller et al. 2018 Nature Communications
TFRcfg.toi           = -param.buffer:param.timeRes:param.segDur + param.buffer;

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

%% preallocations

% main results
allRes = [];

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % get sessions
    sessions         = dir(fullfile(paths.LFPData, subjects{iSub}, 'session*'));
    
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
    
    % channel list
    chanNumbers      = cellfun(@str2num, micro2Region{:, 1}, 'UniformOutput', false);
    chanList         = cat(2, chanNumbers{:, 1});
    
    % information about hemisphere
    hemisphereShort  = cellfun(@(x) x(end), micro2macro{:}, 'UniformOutput', false);
    hemisphere       = replace(hemisphereShort, {'R'; 'L'; 'd'}, {'right'; 'left'; 'NotImplanted'});
    
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
        
        % original session index
        sessIdx = split(sessions(iSess).name, '_');
        sessIdx = str2double(sessIdx{2});
        
        % display session information
        fprintf('\n=============================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);

        % get available data
        mwChanDir  = TG_GetChanDir_20210812(paths.SUData, subjects{iSub}, sessions(iSess).name);  % unit data
        LFPChanDir = TG_GetChanDir_20210812(paths.LFPData, subjects{iSub}, sessions(iSess).name); % LFP data
        
        % brain region sanity check
        if ~isequal(size(LFPChanDir, 1), size(brainRegionsList, 1))
            error('different number of microwire channels than channels in the brain region file');
        end

        %% chest-wise memory performance

        % load behavioral logfile and segment info file
        trialInfo     = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo     = trialInfo.trialInfo;
        segmentInfo   = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'segmentInfo_20230208.mat'));
        
        % each chest's Treasure Hunt trial index
        encTrialIdx   = transpose(rude(trialInfo.numChestsPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numChestsPerTrial))));
        objTrialIdx   = transpose(rude(trialInfo.numObjRecallPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numObjRecallPerTrial))));
        locTrialIdx   = transpose(rude(trialInfo.numLocRecallPerTrial, trialInfo.trialIdx(1:numel(trialInfo.numLocRecallPerTrial))));
        
        % index for encoding segments
        encSegmentIdx = trialInfo.chestIdx;
        
        % prelloacte - object recall
        encObjRec     = cell(size(encTrialIdx, 1), 1);
        objSegmentIdx = nan(size(trialInfo.audioResponse, 1), 1);
        
        % prelloacte - location recall
        encLocRec     = nan(size(encTrialIdx, 1), 1);
        locRec        = nan(size(locTrialIdx, 1), 1);
        locSegmentIdx = nan(size(locTrialIdx, 1), 1);
        
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
        bGoodMemRecAbs  = cat(1, bGoodMemObj, bGoodMemLocAbs);
        recTrialIdx     = cat(1, objTrialIdx, locTrialIdx);
        recSegmentIdx   = cat(1, objSegmentIdx, locSegmentIdx);
        
        % classify object (1) and location recall (2) segments
        encObjOrLoc                = zeros(size(encSegmentIdx, 1), 1);
        encObjOrLoc(objSegmentIdx(~isnan(objSegmentIdx))) = 1;
        encObjOrLoc(locSegmentIdx) = 2;
        recObjOrLoc                = [ones(size(bGoodMemObj, 1), 1); ones(size(bGoodMemLoc, 1), 1) * 2];
        
        % preallocate
        allChanRes = cell(size(LFPChanDir, 1), 1);
        
        %% loop through channels
        parfor iWire = 1:size(LFPChanDir, 1)

            % brain region sanity check
            if ~isequal(str2double(LFPChanDir(iWire).name(5:end)), brainRegionsList{iWire, 2})
                error('mismatch between microwire channel and channel in the brain region file');
            end
            
            % load single-unit classification file
            pathSU = fullfile(paths.SUData, subjects{iSub}, sessions(iSess).name, LFPChanDir(iWire).name);
            try
                c4a = load(fullfile(pathSU, param.c4aName{1, 1}));
            catch
                c4a = load(fullfile(pathSU, param.c4aName{1, 2}));
            end
            c4a = c4a.Cluster4Analysis;
            
            % continue if there are no good clusters found in this channel
            if ~exist(fullfile(pathSU, 'times_datacut.mat'), 'file') || ~any(contains({c4a.Decision}, 'yes'))
                fprintf('No good clusters for this channel.\n');
                continue;
            end
            
            %% load data and perform time-frequency analysis

            % load data
            LFPData         = load(fullfile(paths.LFPData, subjects{iSub}, sessions(iSess).name, LFPChanDir(iWire).name, 'datacutFt2000Hz.mat'));
            
            % loop through encoding and recall
            segmentsOrig    = [];
            chanRes         = [];
            for iCond = 1:size(param.condition, 1)

                % include only encoding or recall periods
                if contains(param.condition{iCond}, 'encoding')
                    segmentsOrig    = segmentInfo.encoding;
                elseif contains(param.condition{iCond}, 'objectRecall')
                    segmentsOrig    = segmentInfo.objRecall;
                elseif contains(param.condition{iCond}, 'locationRecall')
                    segmentsOrig    = segmentInfo.locRecall;
                end

                % cut data
                cfg                     = [];
                cfg.trl                 = round(segmentsOrig / (segmentInfo.fsample / LFPData.fsample)); % adjust sampling rate of segmentInfo to data
                if contains(param.condition{iCond}, 'locationRecall')
                    cfg.trl(:, 1)            = cfg.trl(:, 2) - (param.segDur * LFPData.fsample) - (param.buffer * LFPData.fsample);
                else
                    cfg.trl(:, 1)            = cfg.trl(:, 1) - (param.buffer * LFPData.fsample);
                end
                cfg.trl(:, 2)            = cfg.trl(:, 2) + (param.buffer * LFPData.fsample);
                cfg.trl(:, 3)            = -(param.buffer * LFPData.fsample);
                allSegmentData           = ft_redefinetrial(cfg, LFPData);
                
                % time-frequency analysis
                TFR_allData              = ft_freqanalysis(TFRcfg, allSegmentData);
                
                % log-transformation
                TFR_allData.powspctrm    = log10(TFR_allData.powspctrm);
                
                % reshape and normalize TFR
                squeezedData             = permute(squeeze(TFR_allData.powspctrm), [3, 1, 2]);
                reshData1                = reshape(squeezedData, size(squeezedData, 1) * size(squeezedData, 2), size(squeezedData, 3));
                normData                 = normalize(reshData1, 1);
                reshData2                = reshape(normData, size(squeezedData, 1), size(squeezedData, 2), size(squeezedData, 3));
                permData                 = permute(reshData2, [2, 4, 3, 1]);
                TFR_allData.powspctrm    = permData; % z-transformed powerspectrum
                
                % change axis values to 'trick' fieldtrip to plot all frequencies with same height
                TFR_allData.freq         = 1:numel(TFR_allData.freq); % linearly increasing vector to 'trick' fieldtrip to plot all frequencies with same height
                
                % collect results across channels
                chanRes.idx              = [iSub, sessIdx, iWire];
                chanRes.subject          = subjects{iSub};
                chanRes.session          = sessions(iSess).name;
                chanRes.channel          = LFPChanDir(iWire).name;
                chanRes.brainRegion      = brainRegionsList{iWire, 1};
                
                % collect all results for each condition
                if strcmp(param.condition{iCond}, 'encoding')
                    chanRes.encBothPowspctrm = squeeze(mean(TFR_allData.powspctrm));
                    chanRes.encSuccPowspctrm = squeeze(mean(TFR_allData.powspctrm(bGoodMemEnc, :, :, :)));
                    chanRes.encFailPowspctrm = squeeze(mean(TFR_allData.powspctrm(~bGoodMemEnc, :, :, :)));
                elseif strcmp(param.condition{iCond}, 'objectRecall')
                    chanRes.objBothPowspctrm = squeeze(mean(TFR_allData.powspctrm));
                    chanRes.objSuccPowspctrm = squeeze(mean(TFR_allData.powspctrm(bGoodMemObj, :, :, :)));
                    chanRes.objFailPowspctrm = squeeze(mean(TFR_allData.powspctrm(~bGoodMemObj, :, :, :)));
                elseif strcmp(param.condition{iCond}, 'locationRecall')
                    chanRes.locBothPowspctrm = squeeze(mean(TFR_allData.powspctrm));
                    chanRes.locSuccPowspctrm = squeeze(mean(TFR_allData.powspctrm(bGoodMemLoc, :, :, :)));
                    chanRes.locFailPowspctrm = squeeze(mean(TFR_allData.powspctrm(~bGoodMemLoc, :, :, :)));
                end
            end
            
            % collect results across channels
            allChanRes{iWire, 1} = chanRes;
            
        end
        
        % collect all results across subjects
        allRes = cat(2, allRes, [allChanRes{:}]);
    end
end

% shut down parallel pool
delete(gcp);

% save average results for each channel
save(fullfile(paths.save, 'results'), '-v7.3');

%% additional analysis and figures

%==========================================================================
% in this part of the script power effects during successful and
% failed encoding and object/location recall are compared
%==========================================================================

% load previously saved files
r = load(fullfile(paths.save, 'results'));

% configuration for random number generator
rng(444);

% configurations for ft_freqstatistics
FScfg                     = [];
FScfg.channel             = {'merged_channel'};
FScfg.parameter           = 'powspctrm';
FScfg.method              = 'montecarlo';
FScfg.statistic           = 'ft_statfun_depsamplesT';
FScfg.correctm            = 'cluster';
FScfg.clusterstatistic    = 'maxsum';
FScfg.correcttail         = 'alpha';
FScfg.neighbours          = [];
FScfg.tail                = 0;
FScfg.clustertail         = 0;
FScfg.clusteralpha        = 0.05;
FScfg.alpha               = 0.05;
FScfg.numrandomization    = 1001;
FScfg.ivar                = 1;
FScfg.uvar                = 2;

% configurations for ft_singleplot
Plotcfg                   = [];
Plotcfg.figure            = 'gcf';
Plotcfg.interactive       = 'no';
Plotcfg.xlim              = [-1, 5];
Plotcfg.zlim              = [-0.05, 0.05];
newYLabels                = 2 .^ (ceil((log10(min(r.TFRcfg.foi)) / log10(2))):floor((log10(max(r.TFRcfg.foi)) / log10(2))));  % create new y tick values
Plotcfg.newYLabels        = num2cell(newYLabels);
Plotcfg.yticks            = (log10(newYLabels) - log10(r.param.freqMin) ...
                            + ((log10(r.param.freqMax) - log10(r.param.freqMin)) / (r.param.freqNum - 1))) ...
                            * ((r.param.freqNum - 1) / (log10(r.param.freqMax) - log10(r.param.freqMin))); % create new y ticks positions
Plotcfg.parameter         = 'powspctrm';

% preallocate
biggestClusStat           = cell(size(r.param.condition, 1), 2);

%% loop through encoding, object recall and location recall
for iCond = 1:size(r.param.condition, 1)
    
    % condition-specific settings/variables
    if strcmp(r.param.condition{iCond}, 'encoding')
        succPow             = cat(3, r.allRes.encSuccPowspctrm);
        failPow             = cat(3, r.allRes.encFailPowspctrm);
        FScfg.latency       = [0, 1.5];
    elseif strcmp(r.param.condition{iCond}, 'objectRecall')
        succPow             = cat(3, r.allRes.objSuccPowspctrm);
        failPow             = cat(3, r.allRes.objFailPowspctrm);
        FScfg.latency       = [0, 4];
    elseif strcmp(r.param.condition{iCond}, 'locationRecall')
        succPow             = cat(3, r.allRes.locSuccPowspctrm);
        failPow             = cat(3, r.allRes.locFailPowspctrm);
        FScfg.latency       = [-1, 4];
    end

    % print number of channels
    fprintf('Total number of channels: %d.\n', size(r.allRes, 2));

    % latency start and end time index
    [~, latStartIdx]    = min(abs(r.TFRcfg.toi - FScfg.latency(1, 1)));
    [~, latEndIdx]      = min(abs(r.TFRcfg.toi - FScfg.latency(1, 2)));

    % define variables
    allSessSucc         = {};
    allSessFail         = {};
    numSess             = 0; % counter variable
    
    % loop through subjects
    for iSub = 1:size(r.subjects, 1)
        
        % get sessions
        sessions = dir(fullfile(r.paths.LFPData, r.subjects{iSub}, 'session*'));
        
        % loop through sessions
        for iSess = 1:size(sessions, 1)
            
            % get index for this session
            oneSessIdx               = strcmp({r.allRes.subject}, r.subjects{iSub}) & strcmp({r.allRes.session}, sessions(iSess).name);
            
            % remove session if there are no channels
            if oneSessIdx == 0
                continue;
            end
            
            % number of sessions
            numSess                  = numSess + 1;
            
            % get powerspectrum for successful and failed trials
            oneSessSucc              = mean(succPow(:, :, oneSessIdx), 3, 'omitnan');
            oneSessFail              = mean(failPow(:, :, oneSessIdx), 3, 'omitnan');
            
            % create fieldtrip structure for this session - successful
            oneSessSuccTFR           = [];
            oneSessSuccTFR.label     = {'merged_channel'};
            oneSessSuccTFR.dimord    = 'chan_freq_time';
            oneSessSuccTFR.freq      = 1:size(oneSessSucc, 1);
            oneSessSuccTFR.time      = r.TFRcfg.toi;
            oneSessSuccTFR.powspctrm = permute(oneSessSucc, [3, 1, 2]);

            % create fieldtrip structure for this session - fail
            oneSessFailTFR           = [];
            oneSessFailTFR.label     = {'merged_channel'};
            oneSessFailTFR.dimord    = 'chan_freq_time';
            oneSessFailTFR.freq      = 1:size(oneSessFail, 1);
            oneSessFailTFR.time      = r.TFRcfg.toi;
            oneSessFailTFR.powspctrm = permute(oneSessFail, [3, 1, 2]);
            
            % collect results across sessions
            allSessSucc{1, numSess}  = oneSessSuccTFR;
            allSessFail{1, numSess}  = oneSessFailTFR;
        end
    end
    
    %% grand average over sessions
    cfg         = [];
    cfg.nanmean = 'yes';
    grAvgSucc   = ft_freqgrandaverage(cfg, allSessSucc{:});
    grAvgFail   = ft_freqgrandaverage(cfg, allSessFail{:});
    
    %% ft_freqstatistics
    
    % design matrix for ft_freqstatistics
    FScfg.design             = [ones(1, size(allSessSucc, 2)), ones(1, size(allSessFail, 2)) * 2; ...
        1:size(allSessSucc, 2), 1:size(allSessFail, 2)];
    
    % apply ft_freqstatistics (permutation test)
    statSuccVsFail           = ft_freqstatistics(FScfg, allSessSucc{:}, allSessFail{:});
    
    % subtract powerspectra
    grAvgSuccVsFail                 = grAvgSucc;
    grAvgSuccVsFail.powspctrm       = grAvgSucc.powspctrm - grAvgFail.powspctrm;
    
    % collect time-frequency clusters with smallest p-values for each region
    biggestClusStat{iCond, 1}       = [statSuccVsFail.posclusters(1).prob, statSuccVsFail.posclusters(1).clusterstat, statSuccVsFail.posclusters(1).stddev, statSuccVsFail.posclusters(1).cirange];
    biggestClusStat{iCond, 2}       = [statSuccVsFail.negclusters(1).prob, statSuccVsFail.negclusters(1).clusterstat, statSuccVsFail.negclusters(1).stddev, statSuccVsFail.negclusters(1).cirange];
    
    % create TFR figure - difference
    TFR_Diff                = figure;
    ft_singleplotTFR(Plotcfg, grAvgSuccVsFail);
    yticks(Plotcfg.yticks); % change positions of y ticks
    set(gca, 'YTickLabel', Plotcfg.newYLabels); % change values of y ticks
    title(strcat('(', num2str(size(allSessSucc, 2)), {' '}, 'sessions)'), 'FontSize', 15);
    set(gca, 'fontsize', 15);
    if strcmp(r.param.condition{iCond}, 'encoding')
        xline([0, 1.5], ':', 'LineWidth', 2);
        xlim([-1, 2.5]);
    elseif strcmp(r.param.condition{iCond}, 'objectRecall')
        xline([0, 4], ':', 'LineWidth', 2);
    elseif strcmp(r.param.condition{iCond}, 'locationRecall')
        xline(4, ':', 'LineWidth', 2);

        % change x-axis labels
        xticks                  = get(gca, 'XTick');
        xticklabels             = get(gca, 'XTickLabel');
        xticklabels_numeric     = str2double(xticklabels);
        new_xticklabels_numeric = xticklabels_numeric - 4;
        new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
        set(gca, 'XTickLabel', new_xticklabels);
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    % save figure
    print(TFR_Diff, fullfile(r.paths.save, strcat('TFR_diff_', r.param.condition{iCond})), '-dtiff');
    
    % create TFR figure - successful
    TFR_Succ     = figure;
    ft_singleplotTFR(Plotcfg, grAvgSucc);
    yticks(Plotcfg.yticks); % change positions of y ticks
    set(gca, 'YTickLabel', Plotcfg.newYLabels); % change values of y ticks
    title(strcat('(', num2str(size(allSessSucc, 2)), {' '}, 'sessions)'), 'FontSize', 15);
    set(gca, 'fontsize', 15);
    if strcmp(r.param.condition{iCond}, 'encoding')
        xline([0, 1.5], ':', 'LineWidth', 2);
        xlim([-1, 2.5]);
    elseif strcmp(r.param.condition{iCond}, 'objectRecall')
        xline([0, 4], ':', 'LineWidth', 2);
    elseif strcmp(r.param.condition{iCond}, 'locationRecall')
        xline(4, ':', 'LineWidth', 2);

        % change x-axis labels
        xticks                  = get(gca, 'XTick');
        xticklabels             = get(gca, 'XTickLabel');
        xticklabels_numeric     = str2double(xticklabels);
        new_xticklabels_numeric = xticklabels_numeric - 4;
        new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
        set(gca, 'XTickLabel', new_xticklabels);
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % save figure
    print(TFR_Succ, fullfile(r.paths.save, strcat('TFR_succ_', r.param.condition{iCond})), '-dtiff');
    
    % create TFR figure - fail
    TFR_Fail     = figure;
    ft_singleplotTFR(Plotcfg, grAvgFail);
    yticks(Plotcfg.yticks); % change positions of y ticks
    set(gca, 'YTickLabel', Plotcfg.newYLabels); % change values of y ticks
    title(strcat('(', num2str(size(allSessFail, 2)), {' '}, 'sessions)'), 'FontSize', 15);
    set(gca, 'fontsize', 15);
    if strcmp(r.param.condition{iCond}, 'encoding')
        xline([0, 1.5], ':', 'LineWidth', 2);
        xlim([-1, 2.5]);
    elseif strcmp(r.param.condition{iCond}, 'objectRecall')
        xline([0, 4], ':', 'LineWidth', 2);
    elseif strcmp(r.param.condition{iCond}, 'locationRecall')
        xline(4, ':', 'LineWidth', 2);

        % change x-axis labels
        xticks                  = get(gca, 'XTick');
        xticklabels             = get(gca, 'XTickLabel');
        xticklabels_numeric     = str2double(xticklabels);
        new_xticklabels_numeric = xticklabels_numeric - 4;
        new_xticklabels         = cellstr(num2str(new_xticklabels_numeric));
        set(gca, 'XTickLabel', new_xticklabels);
    end
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    
    % save figure
    print(TFR_Fail, fullfile(r.paths.save, strcat('TFR_fail_', r.param.condition{iCond})), '-dtiff');
end
