%==========================================================================
% This script performs phase locking analysis of units to
% the phase angle of LFP data during encoding and recall segments.
% It also calculates a power spectrum for each segment and the
% spike-triggered average for each unit during encoding and recall.
%
% Tim Guth, 2024
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
paths.phaseData     = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.sprintData    = 'D:\TreasureHunt\MicroSPRiNT_20230713'; % data from SPRiNT algorithm
paths.subInfo       = 'D:\TreasureHunt\SubjectInformation_20210111'; % subject information
paths.save          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % save folder
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
param                 = [];
param.polarHistEdges  = linspace(-pi, pi, 21);
param.polarHistLabels = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
param.c4aName         = {'Cluster4Analysis_LK_20200416_113912', 'Cluster4Analysis_TG_20210315_020120'};
param.filterBand      = [1, 10]; % filter band for phase detection
param.phaseMethod     = 'generalized'; % generalized phase approach
param.segBlDur        = 1.5; % baseline duration in seconds before segments
param.timeBorderSTA   = 1; % time window spike-triggered average in seconds
param.nSur            = 1001;
param.condition       = {'all'; 'baseline'; 'enc_both'; 'enc_succ'; 'enc_fail'; 'rec_both'; 'rec_succ'; 'rec_fail'};

% configurations for random locations in the Treasure Hunt arena to normalize drop error
RandLoccfg            = [];
RandLoccfg.maxR       = 50;
RandLoccfg.minR       = 0;
RandLoccfg.N          = 1000000;
RandLoccfg.centerX    = 370;
RandLoccfg.centerY    = 360;
randLocs              = TG_RandomPointsInCircle(RandLoccfg);

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
allNoSprint    = cell(length(subjects), 1); % excluded because there is no SPRiNT file for this wire

%% loop through subjects
parfor (iSub = 1:length(subjects), 6)
% for iSub = 1:length(subjects)

    % set random number generator
    rng(randSeedNum(iSub, 1), 'twister');

    % set variables
    bGoodMem         = [];
    bGoodMemSeg      = [];
    subRes           = [];

    % for sanity checking
    exByWC           = []; % excluded by wave-clus
    exByVisInsp      = []; % excluded due to visual inspection
    noSprint         = []; % no SPRiNT file existant

    % get sessions
    sessions         = dir(fullfile(paths.SUData, subjects{iSub}, 'session*'));

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

        % original session index
        sessIdx = split(sessions(iSess).name, '_');
        sessIdx = str2double(sessIdx{2});

        % display session information
        fprintf('\n==================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);

        % get available data
        mwChanDir       = TG_GetChanDir_20210812(paths.SUData, subjects{iSub}, sessions(iSess).name);  % unit data
        phaseChanDir    = TG_GetChanDir_20210812(paths.phaseData, subjects{iSub}, sessions(iSess).name); % phase data

        % channel names sanity check
        if any(~strcmp({mwChanDir.name}, {phaseChanDir.name}))
            error('different channels for unit data and phase data');
        end

        % brain region sanity check
        if ~isequal(size(phaseChanDir, 1), size(brainRegionsList, 1))
            error('different number of microwire channels than channels in the brain region file');
        end

        %% chest-wise memory performance

        % load behavioral logfile and segment info file
        trialInfo     = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo     = trialInfo.trialInfo;
        segmentInfo   = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'segmentInfo_20230208.mat'));

        % each chest's TreasureHunt trial index
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

        %% loop through channels
        for iWire = 1:size(phaseChanDir, 1)

            % print channel name
            fprintf('\tChannel name: %s.\n', mwChanDir(iWire).name);

            % path for this channel's unit data
            mwChanPath          = fullfile(paths.SUData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name);

            %% LFP data

            % load LFP data
            LFPdata             = load(fullfile(paths.LFPData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'datacutFt2000Hz.mat'));

            %% artifact data

            % load data for artifact removal
            artifactData        = load(fullfile(paths.artifact, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'detectedArtifacts.mat'), 'bArtifact');
            bArtifact           = artifactData.bArtifact;

            % artifact ratio
            artifactRatio       = sum(bArtifact) / size(bArtifact, 2);

            %% generalized phase data

            % load phase data
            phaseDataName       = strcat('datacutFt2000HzBP_', param.phaseMethod, '_', regexprep(num2str(param.filterBand), '\s+', '_'), '.mat');
            phaseData           = load(fullfile(paths.phaseData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, phaseDataName));

            % instantaneous complex value
            if ~isempty(phaseData.trial)
                instComplex     = phaseData.trial{1, 1};
            else
                continue;
            end

            % instantaneous phase
            instPhase           = angle(instComplex);

            % instantaneous power
            instPower           = abs(instComplex) .^ 2;

            % mean log power across session
            avgLogPower         = mean(log10(instPower));

            %% power index for this wire

            % get power quantiles
            powerMedian         = median(instPower(~bArtifact));

            % creater power index
            powerIdx            = double(instPower > powerMedian);

            %% SPRiNT data

            % load SPRiNT data
            try
                sprintData      = load(fullfile(paths.sprintData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'datacutSprint.mat'));
            catch
                noSprint        = cat(1, noSprint, [iSub, iSess, iWire]); % bookkeeping
                fprintf('No sprint file for %s, session %d, channel %d', subjects{iSub}, iSess, iWire);
                continue;
            end

            % extract SPRiNT time and slope
            sprintTimeOrig      = cat(2, sprintData.SPRiNT.channel.aperiodics.time);
            sprintSlopeOrig     = cat(2, sprintData.SPRiNT.channel.aperiodics.exponent);

            % SPRiNT time in samples
            sprintSamples       = sprintTimeOrig * phaseData.fsample;

            % extract sprint times with theta peaks
            centerFrequency     = cat(1, sprintData.SPRiNT.channel.peaks.center_frequency); % center frequency of peaks
            thetaPeakIdx        = centerFrequency >= param.filterBand(1, 1) & centerFrequency <= param.filterBand(1, 2); % extract only peaks in theta band

            % get time of all peaks
            peakTime            = cat(1, sprintData.SPRiNT.channel.peaks.time);

            % get only times of theta peaks
            timesWithTheta      = unique(peakTime(thetaPeakIdx));

            % create index for theta oscillations in sprint signal
            [~, thetaIdxOrig]     	        = min(abs(sprintTimeOrig - timesWithTheta), [], 2);
            bThetaOrig                      = zeros(size(sprintTimeOrig));
            bThetaOrig(1, thetaIdxOrig)     = ones;

            % overall time coverage of SPRiNT windows
            sprintWindow        = sprintData.SPRiNT.options.WinLength + (sprintData.SPRiNT.options.WinAverage - 1) * (sprintData.SPRiNT.options.WinLength * (1 - (sprintData.SPRiNT.options.WinOverlap / 100)));

            % samples of overall SPRiNT windows
            sprintWindowsSmp    = sprintSamples' + ((-(sprintWindow / 2 * phaseData.fsample) + 1):1:((sprintWindow / 2 * phaseData.fsample)));

            % find artifacts 
            slopeBArtifact      = any(artifactData.bArtifact(sprintWindowsSmp), 2)';

            % remove artifact data
            sprintSlopeOrig(slopeBArtifact) = NaN;
            bThetaOrig(slopeBArtifact)      = NaN;

            % upsample SPRiNT results to match phase data sampling rate
            sprintFsample       = 1 / mode(diff(sprintTimeOrig));
            upsamplingRatio     = phaseData.fsample / sprintFsample;
            sprintSlope         = repelem(sprintSlopeOrig, upsamplingRatio);
            bTheta              = repelem(bThetaOrig, upsamplingRatio);

            % add  NaNs to start and end borders of SPRiNT data for padding
            missingSprint       = size(instComplex, 2) - size(sprintSlope, 2);
            sprintSlope         = cat(2, repelem(NaN, missingSprint / 2), sprintSlope, repelem(NaN, missingSprint / 2));
            bTheta              = cat(2, repelem(NaN, missingSprint / 2), bTheta, repelem(NaN, missingSprint / 2));

            %% slope index for this wire

            % get slope quantiles
            slopeMedian         = median(sprintSlopeOrig, 'omitnan');

            % create slope index
            slopeIdx            = double(sprintSlope > slopeMedian);

            % exclude NaNs in slope data
            slopeIdx(isnan(sprintSlope))    = NaN;

            %% sanity checks

            % check whether LFP data and phase data have the same sampling rate
            if phaseData.fsample ~= LFPdata.fsample
                error('phase data has different sampling rate than LFP data');
            end

            % sanity check for artifact detection data
            if size(bArtifact, 2) ~= size(instComplex, 2)
                error('phase data does not match artifact detection data');
            end

            %% wave-clus output

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

            %% segment indices

            % endoding
            encSegments         = segmentInfo.encoding / (segmentInfo.fsample / phaseData.fsample);
            encSegments(:, 1)   = ceil(encSegments(:, 1));
            encSegments(:, 2)   = floor(encSegments(:, 2));
            encSegmentDataIdx   = zeros(size(instPhase));

            % loop through encoding segments
            for iEnc = 1:size(encSegments, 1)
                encSegmentDataIdx(encSegments(iEnc, 1):encSegments(iEnc, 2)) = iEnc;
            end

            % recall
            segmentInfo.recall  = cat(1, segmentInfo.objRecall, segmentInfo.locRecall);
            recSegments         = segmentInfo.recall / (segmentInfo.fsample / phaseData.fsample);
            recSegments(:, 1)   = ceil(recSegments(:, 1));
            recSegments(:, 2)   = floor(recSegments(:, 2));
            recSegmentDataIdx   = zeros(size(instPhase));

            % loop through recall segments
            for iRec = 1:size(recSegments, 1)
                recSegmentDataIdx(recSegments(iRec, 1):recSegments(iRec, 2)) = iRec;
            end

            % baseline
            blSegmentDataIdx    = bwlabel(~encSegmentDataIdx & ~recSegmentDataIdx);

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

                %% cluster properties

                % get cluster times and samples
                thisClusOrig                        = t.cluster_class(t.cluster_class(:, 1) == iClus, :); % cluster-number, time in ms
                thisClusSamplesOrig                 = round(thisClusOrig(:, 2) * phaseData.fsample / 1000); % convert from ms to phase angle samples

                % index for this cluster
                bThisClusOrig                       = false(size(instPhase));
                bThisClusOrig(thisClusSamplesOrig)  = true;

                %% check whether cluster is from a single unit

                % signal-to-noise ratio
                thisClusSpikes          = t.spikes(t.cluster_class(:, 1) == iClus, :);
                thisClusMeanSpike       = mean(thisClusSpikes);
                thisClusAmplitude       = abs(thisClusMeanSpike(:, 20));
                thisClusNoise           = std(thisClusSpikes(:, 1));
                thisClusSNR             = thisClusAmplitude / thisClusNoise;

                % inter-spike interval
                thisClusISI             = diff(thisClusOrig(:, 2));
                thisClusIsiThreeMs      = thisClusISI < 3;
                thisClusRatioThreeMs    = sum(thisClusIsiThreeMs) / size(thisClusIsiThreeMs, 1);

                % determine whether it is a single unit
                thisClusSingleUnit      = false;
                if thisClusSNR > 3 && thisClusRatioThreeMs < 0.01
                    thisClusSingleUnit  = true;
                end

                %% remove artifacts

                % IEDs/artifact removal
                thisClusArtifactIdx         = bArtifact(bThisClusOrig)';
                thisClus                    = thisClusOrig(~thisClusArtifactIdx, :);
                thisClusSamples             = thisClusSamplesOrig(~thisClusArtifactIdx);

                % index for this cluster
                bThisClus                   = false(size(instPhase));
                bThisClus(thisClusSamples)  = true;

                %% spike attributes

                % encoding and retrieval index
                thisClusEncRecIdx           = zeros(size(thisClusSamples));
                thisClusEncRecIdx(encSegmentDataIdx(bThisClus) > 0) = 1;
                thisClusEncRecIdx(recSegmentDataIdx(bThisClus) > 0) = 2;

                %% phase angles

                % prealloacte
                condDataOrigIdx         = [];
                segSpikeArtifact        = [];
                segDur                  = [];
                segFr                   = [];
                segBlFr                 = [];
                segPhase                = [];
                segSpikeTime            = [];
                segSpikePower           = [];
                segThetaIdx             = [];
                segPowerIdx             = [];
                segSlopeIdx             = [];
                segSpikeSlope           = [];
                STA                     = [];

                % all
                allDur                  = [];
                allFr                   = [];
                allComplex              = [];
                allSpikeArtifact        = [];
                allPowerIdx             = [];
                allSlopeIdx             = [];
                allThetaIdx             = [];
                allSpikeSlope           = [];
                allZval                 = [];
                allMrv                  = [];
                allRank                 = [];

                % baseline
                blDur                   = [];
                blFr                    = [];
                blSpikeArtifact         = [];
                blSpikePower            = [];
                blZval                  = [];
                blMrv                   = [];
                blRank                  = [];

                % encoding
                encStartStop            = [];
                encDur                  = [];
                encFr                   = [];
                encBlFr                 = [];
                encPhase                = [];
                encSpikeTime            = [];
                encSpikeArtifact        = [];
                encSpikePower           = [];
                encPowerIdx             = [];
                encSlopeIdx             = [];
                encThetaIdx             = [];
                encSpikeSlope           = [];
                encBothZval             = [];
                encBothMrv              = [];
                encBothRank             = [];
                encSuccRank             = [];
                encFailRank             = [];
                encBothSTA              = [];
                encSuccSTA              = [];
                encFailSTA              = [];

                % recall
                recStartStop            = [];
                recDur                  = [];
                recFr                   = [];
                recBlFr                 = [];
                recPhase                = [];
                recSpikeTime            = [];
                recSpikeArtifact        = [];
                recSpikePower           = [];
                recPowerIdx             = [];
                recSlopeIdx             = [];
                recThetaIdx             = [];
                recSpikeSlope           = [];
                recBothZval             = [];
                recBothMrv              = [];
                recBothRank             = [];
                recSuccRank             = [];
                recFailRank             = [];
                recBothSTA              = [];
                recSuccSTA              = [];
                recFailSTA              = [];

                % loop through conditions
                for iCond = 1:size(param.condition, 1)

                    % segment index
                    if contains(param.condition{iCond}, 'all')
                        condDataOrigIdx = true(size(instPhase));
                        bGoodMemSeg     = true;
                    elseif contains(param.condition{iCond}, 'baseline')
                        condDataOrigIdx = blSegmentDataIdx;
                        bGoodMemSeg     = true(max(blSegmentDataIdx), 1);
                    elseif contains(param.condition{iCond}, 'enc')
                        condDataOrigIdx = encSegmentDataIdx;
                        bGoodMem        = bGoodMemEnc;
                    elseif contains(param.condition{iCond}, 'rec')
                        condDataOrigIdx = recSegmentDataIdx;
                        bGoodMem        = bGoodMemRec;
                    end

                    % memory performance index
                    if contains(param.condition{iCond}, 'both')
                        bGoodMemSeg = true(size(bGoodMem, 1), 1);
                    elseif contains(param.condition{iCond}, 'succ')
                        bGoodMemSeg = bGoodMem;
                    elseif contains(param.condition{iCond}, 'fail')
                        bGoodMemSeg = ~bGoodMem;
                    end

                    % only include segments with specific memory performance
                    thisCondDataOrigIdx         = condDataOrigIdx;
                    thisCondDataOrigIdx(~any(thisCondDataOrigIdx == find(bGoodMemSeg))) = 0;

                    % exclude artifacts
                    thisCondDataIdx             = thisCondDataOrigIdx;
                    thisCondDataIdx(bArtifact)  = 0;

                    % extract cluster data for this condition
                    thisCondClus                = bThisClus(thisCondDataIdx > 0);

                    % extract data for this condition
                    thisCondPhase               = instPhase(thisCondDataIdx > 0);

                    % empirical spike phases of this condition
                    catPhases                   = thisCondPhase(thisCondClus)';

                    % get empirical mean resultant vector
                    [segMrv, ~]                 = circ_axialmean(catPhases, [], 1);

                    % get empirical Rayleigh z-value
                    if ~isempty(catPhases)
                        [~, segZval]            = circ_rtest(catPhases);
                    else
                        segZval                 = NaN;
                    end

                    % get surrogate MRV
                    surSegMrv       = NaN(1, param.nSur);
                    surSegZval      = NaN(1, param.nSur);
                    for iSur = 1:param.nSur

                        %% create surrogates

                        % random number for circular shift
                        surShiftSmpl        = randi(size(thisCondPhase, 2), 1, 1);

                        % circularly shift phases of this condition
                        surThisCondPhase    = circshift(thisCondPhase, surShiftSmpl);

                        % surrogate phases
                        surCatPhases        = surThisCondPhase(thisCondClus)';

                        % get surrogate Rayleigh z-values
                        try
                            [surSegMrv(:, iSur), ~]     = circ_axialmean(surCatPhases, [], 1);
                            [~, surSegZval(:, iSur)]    = circ_rtest(surCatPhases);

                        catch
                            surSegMrv(:, iSur)          = NaN;
                            surSegZval(:, iSur)         = NaN;
                        end
                    end

                    % get rank of empirical MRV in surrogate distribution
                    segRank             = sum(segZval > surSegZval, 2) / size(surSegZval, 2);

                    %% get additional segment-related data

                    % get information about segments
                    thisSegments        = unique(thisCondDataOrigIdx(thisCondDataOrigIdx > 0))';

                    % preallocate
                    segSpikeSamplesOrig = cell(size(thisSegments, 1), 1);
                    segSpikeSamples     = cell(size(thisSegments, 1), 1);
                    segComplex          = cell(size(thisSegments, 1), 1);
                    segPhase            = cell(size(thisSegments, 1), 1);
                    segSpikeTime        = cell(size(thisSegments, 1), 1);
                    segSpikePower       = cell(size(thisSegments, 1), 1);
                    segSpikeSlope       = cell(size(thisSegments, 1), 1);
                    segSpikeArtifact    = cell(size(thisSegments, 1), 1);
                    segPowerIdx         = cell(size(thisSegments, 1), 1);
                    segSlopeIdx         = cell(size(thisSegments, 1), 1);
                    segThetaIdx         = cell(size(thisSegments, 1), 1);
                    segStartStop        = nan(size(thisSegments, 1), 2);
                    segDur              = nan(size(thisSegments, 1), 1);
                    segFr               = nan(size(thisSegments, 1), 1);
                    segBlFr             = nan(size(thisSegments, 1), 1);

                    % loop through segments
                    for iSeg = 1:size(thisSegments, 1)

                        % segment index
                        bThisSegment                    = thisCondDataOrigIdx == thisSegments(iSeg);

                        % get spike indices
                        segSpikeSamplesOrig{iSeg, 1}    = find(bThisClusOrig & bThisSegment)';
                        segSpikeSamples{iSeg, 1}        = find(bThisClus & bThisSegment)';

                        % get complex signal results
                        segComplex{iSeg, 1}             = instComplex(segSpikeSamples{iSeg, 1})';

                        % get phases
                        segPhase{iSeg, 1}               = instPhase(segSpikeSamples{iSeg, 1})';

                        % get spike times within segment
                        segSpikeTime{iSeg, 1}           = (segSpikeSamples{iSeg, 1} - find(bThisSegment, 1, 'first')) / phaseData.fsample * 1000; % convert to ms

                        % get artifact ratio
                        segSpikeArtifact{iSeg, 1}       = bArtifact(segSpikeSamplesOrig{iSeg, 1})';

                        % get power
                        segSpikePower{iSeg, 1}          = instPower(segSpikeSamples{iSeg, 1})';

                        % get slope
                        segSpikeSlope{iSeg, 1}          = sprintSlope(segSpikeSamples{iSeg, 1})';

                        % get power index
                        segPowerIdx{iSeg, 1}            = powerIdx(segSpikeSamples{iSeg, 1})';

                        % get slope index
                        segSlopeIdx{iSeg, 1}            = slopeIdx(segSpikeSamples{iSeg, 1})';

                        % get theta index
                        segThetaIdx{iSeg, 1}            = bTheta(segSpikeSamples{iSeg, 1})';

                        %% firing rate

                        % segment duration
                        segSmplNum              = find(bThisSegment, 1, 'last') - find(bThisSegment, 1, 'first'); % start and end sample difference
                        segStartStop(iSeg, :)   = round([find(bThisSegment, 1, 'first'), find(bThisSegment, 1, 'last')] / phaseData.fsample * 1000);
                        segDur(iSeg, 1)         = segSmplNum / phaseData.fsample; % divide by sampling rate

                        % firing rate during segment
                        segFr(iSeg, 1)          = size(segSpikeSamplesOrig{iSeg, 1}, 1) / segDur(iSeg, 1);

                        % firing rate during baseline segment before
                        if contains(param.condition{iCond}, {'enc', 'rec'})
                            segBlSmpl           = param.segBlDur * phaseData.fsample;
                            segBlSpikeIdx       = find(thisClusSamplesOrig >= (find(bThisSegment, 1, 'first') - segBlSmpl) & thisClusSamplesOrig <= find(bThisSegment, 1, 'first'));
                            segBlFr(iSeg, 1)    = size(segBlSpikeIdx, 1) / param.segBlDur;
                        else
                            segBlFr(iSeg, 1)    = NaN;
                        end
                    end

                    %% spike-triggered average

                    % only during encoding and recall
                    if contains(param.condition{iCond}, {'enc', 'rec'})

                        % get spike samples
                        allSegSpikeSamples = cat(1, segSpikeSamples{:});

                        % get sample window
                        sampleWindowSTA    = -param.timeBorderSTA * LFPdata.fsample:1:param.timeBorderSTA * LFPdata.fsample;
                        allSamplesSTA      = allSegSpikeSamples + sampleWindowSTA;

                        % get spike-triggered LFP data
                        allLFPdataSTA      = LFPdata.trial{1, 1}(allSamplesSTA);

                        % get STA
                        STA                = mean(allLFPdataSTA, 1);
                    end

                    %% collect results

                    % get results for each condition
                    if strcmp(param.condition{iCond}, 'all')
                        allDur              = segDur;
                        allFr               = segFr;
                        allComplex          = segComplex;
                        allSpikeArtifact    = segSpikeArtifact;
                        allPowerIdx         = segPowerIdx;
                        allThetaIdx         = segThetaIdx;
                        allSpikeSlope       = segSpikeSlope;
                        allSlopeIdx         = segSlopeIdx;
                        allZval             = segZval;
                        allMrv              = segMrv;
                        allRank             = segRank;
                    elseif strcmp(param.condition{iCond}, 'baseline')
                        blDur               = segDur;
                        blFr                = segFr;
                        blStartStop         = segStartStop;
                        blPhase             = segPhase;
                        blSpikeArtifact     = segSpikeArtifact;
                        blSpikePower        = segSpikePower;
                        blZval              = segZval;
                        blMrv               = segMrv;
                        blRank              = segRank;
                    elseif strcmp(param.condition{iCond}, 'enc_both')
                        encDur              = segDur;
                        encFr               = segFr;
                        encBlFr             = segBlFr;
                        encStartStop        = segStartStop;
                        encPhase            = segPhase;
                        encSpikeTime        = segSpikeTime;
                        encSpikeArtifact    = segSpikeArtifact;
                        encSpikePower       = segSpikePower;
                        encPowerIdx         = segPowerIdx;
                        encThetaIdx         = segThetaIdx;
                        encSpikeSlope       = segSpikeSlope;
                        encSlopeIdx         = segSlopeIdx;
                        encBothZval         = segZval;
                        encBothMrv          = segMrv;
                        encBothRank         = segRank;
                        encBothSTA          = STA;
                    elseif strcmp(param.condition{iCond}, 'enc_succ')
                        encSuccRank         = segRank;
                        encSuccSTA          = STA;
                    elseif strcmp(param.condition{iCond}, 'enc_fail')
                        encFailRank         = segRank;
                        encFailSTA          = STA;
                    elseif strcmp(param.condition{iCond}, 'rec_both')
                        recDur              = segDur;
                        recFr               = segFr;
                        recBlFr             = segBlFr;
                        recStartStop        = segStartStop;
                        recPhase            = segPhase;
                        recSpikeTime        = segSpikeTime;
                        recSpikeArtifact    = segSpikeArtifact;
                        recSpikePower       = segSpikePower;
                        recPowerIdx         = segPowerIdx;
                        recThetaIdx         = segThetaIdx;
                        recSpikeSlope       = segSpikeSlope;
                        recSlopeIdx         = segSlopeIdx;
                        recBothZval         = segZval;
                        recBothMrv          = segMrv;
                        recBothRank         = segRank;
                        recBothSTA          = STA;
                    elseif strcmp(param.condition{iCond}, 'rec_succ')
                        recSuccRank         = segRank;
                        recSuccSTA          = STA;
                    elseif strcmp(param.condition{iCond}, 'rec_fail')
                        recFailRank         = segRank;
                        recFailSTA          = STA;
                    end
                end

                %% collect information for this unit

                % basics
                unitRes                     = [];
                unitRes.idx                 = [iSub, sessIdx, iWire, iClus];
                unitRes.brainRegion         = brainRegionsList{iWire, 1};
                unitRes.labelEnglish        = trialInfo.TREASURE_LABEL_EN;
                unitRes.labelGerman         = trialInfo.TREASURE_LABEL_GER;
                unitRes.bSingleUnit         = thisClusSingleUnit;

                % all
                unitRes.allArtRatio         = artifactRatio;
                unitRes.allComplex          = allComplex{:};
                unitRes.allEncRecIdx        = thisClusEncRecIdx;
                unitRes.allSpikeArtifact    = allSpikeArtifact{:};
                unitRes.allFr               = allFr;
                unitRes.allZval             = allZval;
                unitRes.allMrv              = allMrv;
                unitRes.allRank             = allRank;
                unitRes.allSpikeSlope       = allSpikeSlope{:};
                unitRes.allPowerIdx         = allPowerIdx{:};
                unitRes.allSlopeIdx         = allSlopeIdx{:};
                unitRes.allThetaIdx         = allThetaIdx{:};
                unitRes.allTime             = thisClus(:, 2);

                % baseline
                unitRes.blSpikeArtifact     = blSpikeArtifact;
                unitRes.blSpikePower        = blSpikePower;
                unitRes.blDur               = blDur;
                unitRes.blFr                = blFr;
                unitRes.blZval              = blZval;
                unitRes.blMrv               = blMrv;
                unitRes.blRank              = blRank;

                % encoding
                unitRes.encStartStop        = encStartStop;
                unitRes.encDur              = encDur;
                unitRes.encFr               = encFr;
                unitRes.encBlFr             = encBlFr;
                unitRes.bGoodMemEnc         = bGoodMemEnc;
                unitRes.bGoodMemEncAbs      = bGoodMemEncAbs;
                unitRes.encSpikeArtifact    = encSpikeArtifact;
                unitRes.encTrialIdx         = encTrialIdx;
                unitRes.encSegmentIdx       = encSegmentIdx;
                unitRes.encObjOrLoc         = encObjOrLoc;
                unitRes.encPhase            = encPhase;
                unitRes.encSpikeTime        = encSpikeTime;
                unitRes.encSpikePower       = encSpikePower;
                unitRes.encSpikeSlope       = encSpikeSlope;
                unitRes.encPowerIdx         = encPowerIdx;
                unitRes.encSlopeIdx         = encSlopeIdx;
                unitRes.encThetaIdx         = encThetaIdx;
                unitRes.encBothZval         = encBothZval;
                unitRes.encBothMrv          = encBothMrv;
                unitRes.encBothRank         = encBothRank;
                unitRes.encSuccRank         = encSuccRank;
                unitRes.encFailRank         = encFailRank;
                unitRes.encBothSTA          = encBothSTA;
                unitRes.encSuccSTA          = encSuccSTA;
                unitRes.encFailSTA          = encFailSTA;

                % recall
                unitRes.recStartStop        = recStartStop;
                unitRes.recDur              = recDur;
                unitRes.recFr               = recFr;
                unitRes.recBlFr             = recBlFr;
                unitRes.bGoodMemRec         = bGoodMemRec;
                unitRes.bGoodMemRecAbs      = bGoodMemRecAbs;
                unitRes.recSpikeArtifact    = recSpikeArtifact;
                unitRes.recTrialIdx         = recTrialIdx;
                unitRes.recSegmentIdx       = recSegmentIdx;
                unitRes.recObjOrLoc         = recObjOrLoc;
                unitRes.recPhase            = recPhase;
                unitRes.recSpikeTime        = recSpikeTime;
                unitRes.recSpikePower       = recSpikePower;
                unitRes.recSpikeSlope       = recSpikeSlope;
                unitRes.recPowerIdx         = recPowerIdx;
                unitRes.recSlopeIdx         = recSlopeIdx;
                unitRes.recThetaIdx         = recThetaIdx;
                unitRes.recBothZval         = recBothZval;
                unitRes.recBothMrv          = recBothMrv;
                unitRes.recBothRank         = recBothRank;
                unitRes.recSuccRank         = recSuccRank;
                unitRes.recFailRank         = recFailRank;
                unitRes.recBothSTA          = recBothSTA;
                unitRes.recSuccSTA          = recSuccSTA;
                unitRes.recFailSTA          = recFailSTA;

                % collapse across units
                subRes                      = cat(1, subRes, unitRes);
            end
        end
    end

    % collect results across subjects
    allRes{iSub, 1}         = subRes;
    allExByWC{iSub, 1}      = exByWC;
    allExByVisInsp{iSub, 1} = exByVisInsp;
    allNoSprint{iSub, 1}    = noSprint;
end

% unnest result cells
allRes                  = cat(1, allRes{:});
allExByWC               = cat(1, allExByWC{:});
allExByVisInsp          = cat(1, allExByVisInsp{:});
allNoSprint             = cat(1, allNoSprint{:});
allIdx                  = cat(1, allRes.idx);

% save important output
save(fullfile(paths.save, 'allResults'), 'allRes', 'allExByWC', 'allExByVisInsp', 'allNoSprint', 'allIdx', '-v7.3');

%==========================================================================
% additional analyses at the unit level
%==========================================================================

% load previously saved phase analysis results
r = load(fullfile(paths.save, 'allResults.mat'));
s = load(fullfile(paths.save, 'settings.mat'));
fprintf('Total number of cells: %d.\n', size(r.allRes, 1));

% settings
rng(444, 'twister');

% define variables
conditions              = {'both'; 'succ'; 'fail'}; % memory performance
parforRes               = [];
parforRes.param         = s.param;
parforRes.paths         = s.paths;
parforRes.subjects      = s.subjects;
parforRes.randSeedNum   = s.randSeedNum;
parforRes.idx           = {r.allRes.idx}';
parforRes.brainRegion   = {r.allRes.brainRegion}';

% cfg for STA bandpass filter
STAcfg                  = [];
STAcfg.bpfilter         = 'yes';
STAcfg.bpfilttype       = 'fir';
STAcfg.bpfreq           = s.param.filterBand;
STAcfg.demean           = 'yes';

%% define variables

% time axis for STA
timeAxisSTA      = linspace(-s.param.timeBorderSTA, s.param.timeBorderSTA, size(r.allRes(1).encBothSTA, 2));

% encoding
allbGoodMemEnc   = {r.allRes.bGoodMemEnc}';
allEncPhase      = {r.allRes.encPhase}';
allEncPowerIdx   = {r.allRes.encPowerIdx}';

% recall
allbGoodMemRec   = {r.allRes.bGoodMemRec}';
allRecPhase      = {r.allRes.recPhase}';
allRecPowerIdx   = {r.allRes.recPowerIdx}';

% preallocate
bExclude         = false(size(r.allRes, 1), 1);

%% loop through different memory conditions
for iCond = 1:size(conditions, 1)

    % include all trials
    if strcmp(conditions{iCond}, 'both')

        % encoding
        encRanks        = [r.allRes.encBothRank]';
        encSTA          = cat(1, r.allRes.encBothSTA);
        encPhases       = allEncPhase;
        encPowerIndices = allEncPowerIdx;

        % recall
        recRanks        = [r.allRes.recBothRank]';
        recSTA          = cat(1, r.allRes.recBothSTA);
        recPhases       = allRecPhase;
        recPowerIndices = allRecPowerIdx;

        % include only successful trials
    elseif strcmp(conditions{iCond}, 'succ')

        % encoding
        encRanks        = [r.allRes.encSuccRank]';
        encSTA          = cat(1, r.allRes.encSuccSTA);
        encPhases       = cellfun(@(x, y) x(y, :), allEncPhase, allbGoodMemEnc, 'UniformOutput', false);
        encPowerIndices = cellfun(@(x, y) x(y, :), allEncPowerIdx, allbGoodMemEnc, 'UniformOutput', false);

        % recall
        recRanks        = [r.allRes.recSuccRank]';
        recSTA          = cat(1, r.allRes.recSuccSTA);
        recPhases       = cellfun(@(x, y) x(y, :), allRecPhase, allbGoodMemRec, 'UniformOutput', false);
        recPowerIndices = cellfun(@(x, y) x(y, :), allRecPowerIdx, allbGoodMemRec, 'UniformOutput', false);

        % include only unsuccessful trials
    elseif strcmp(conditions{iCond}, 'fail')

        % encoding
        encRanks        = [r.allRes.encFailRank]';
        encSTA          = cat(1, r.allRes.encFailSTA);
        encPhases       = cellfun(@(x, y) x(~y, :), allEncPhase, allbGoodMemEnc, 'UniformOutput', false);
        encPowerIndices = cellfun(@(x, y) x(~y, :), allEncPowerIdx, allbGoodMemEnc, 'UniformOutput', false);

        % recall
        recRanks        = [r.allRes.recFailRank]';
        recSTA          = cat(1, r.allRes.recFailSTA);
        recPhases       = cellfun(@(x, y) x(~y, :), allRecPhase, allbGoodMemRec, 'UniformOutput', false);
        recPowerIndices = cellfun(@(x, y) x(~y, :), allRecPowerIdx, allbGoodMemRec, 'UniformOutput', false);
    end

    %% preallocate

    % encoding
    allEncMrv      = cell(size(r.allRes, 1), 1);
    allEncMean     = cell(size(r.allRes, 1), 1);
    allEncRank     = cell(size(r.allRes, 1), 1);

    % recall
    allRecMrv     = cell(size(r.allRes, 1), 1);
    allRecMean    = cell(size(r.allRes, 1), 1);
    allRecRank    = cell(size(r.allRes, 1), 1);

    % rank of Watson-William permutation test
    allWwRank     = cell(size(r.allRes, 1), 1);

    %% loop through cells
    parfor iCell = 1:size(r.allRes, 1)

        % set random number generator
        rng(parforRes.randSeedNum(iCell, 1), 'twister');

        % print channel name
        fprintf('\tCell number: %d.\n', iCell);

        % get phase angles of this cell
        cellEncPhase    = cat(1, encPhases{iCell, 1}{:});
        cellRecPhase    = cat(1, recPhases{iCell, 1}{:});

        % include only units with >= 25 spikes during encoding and during recall
        % and with spikes in >= 20 percent of all encoding and recall segments
        if size(cellEncPhase, 1) < 25 || ...
                size(cellRecPhase, 1) < 25 || ...
                (sum(~cellfun(@isempty, encPhases{iCell, 1})) / size(encPhases{iCell, 1}, 1)) < 0.2 || ...
                (sum(~cellfun(@isempty, recPhases{iCell, 1})) / size(recPhases{iCell, 1}, 1)) < 0.2

            % encoding
            allEncMrv{iCell, 1}     = NaN;
            allEncMean{iCell, 1}    = NaN;
            allEncRank{iCell, 1}    = NaN;

            % recall
            allRecMrv{iCell, 1}     = NaN;
            allRecMean{iCell, 1}    = NaN;
            allRecRank{iCell, 1}    = NaN;

            % rank of Watson-Williams test
            allWwRank{iCell, 1}     = NaN;

            % exclude cell
            bExclude(iCell, 1)      = true;
            continue;
        end

        %% bandpass filter STA data

        % encoding STA in fieldtrip structure
        encSTAFt            = [];
        encSTAFt.trial      = {encSTA(iCell, :)};
        encSTAFt.time       = {timeAxisSTA};
        encSTAFt.label      = {'STA'};
        encSTAFilt          = ft_preprocessing(STAcfg, encSTAFt);

        % recall STA in fieldtrip structure
        recSTAFt            = [];
        recSTAFt.trial      = {recSTA(iCell, :)};
        recSTAFt.time       = {timeAxisSTA};
        recSTAFt.label      = {'STA'};
        recSTAFilt          = ft_preprocessing(STAcfg, recSTAFt);

        %% perform Watson-Williams-permutation test

        % concatenate all encoding and recall phases
        allPhases           = [encPhases{iCell, 1}; recPhases{iCell, 1}];

        % class variable defining encoding (1) and recall (2)
        class               = [ones(size(encPhases{iCell, 1}, 1), 1); ones(size(recPhases{iCell, 1}, 1), 1) * 2];

        % empirical absolute mean angles and mean resultant vectors
        [cellEncMrv, cellEncMean] = circ_axialmean(cat(1, allPhases{class == 1}));
        [cellRecMrv, cellRecMean] = circ_axialmean(cat(1, allPhases{class == 2}));

        % empirical Watson-Williams test (see Van Rullen et al. 2016)
        [~, wwTable]        = circ_wwtest(cat(1, allPhases{class == 1}), cat(1, allPhases{class == 2}));
        wwF                 = wwTable{2, 5};

        % surrogate tests
        surWwF    = zeros(parforRes.param.nSur, 1);
        for iSur = 1:parforRes.param.nSur

            % shuffle class for surrogate dataset
            surClass            = datasample(class, numel(class), 'replace', false);

            % surrogate Watson-Williams test
            [~, surWwTable]     = circ_wwtest(cat(1, allPhases{surClass == 1}), cat(1, allPhases{surClass == 2}));
            surWwF(iSur, 1)     = surWwTable{2, 5};
        end

        % rank of empirical WW-test f-value in surrogate dataset
        wwRank      = sum(wwF > surWwF) / parforRes.param.nSur;

        % get p-values
        encP        = 1 - encRanks(iCell, 1);
        recP        = 1 - recRanks(iCell, 1);
        wwP         = 1 - wwRank;

        % correct low p-values of permutation tests (see Phipson and Smyth, 2010)
        encP(encP < (1 / parforRes.param.nSur))     = 1 / parforRes.param.nSur;
        recP(recP < (1 / parforRes.param.nSur))     = 1 / parforRes.param.nSur;
        wwP(wwP < (1 / parforRes.param.nSur))       = 1 / parforRes.param.nSur;

        %% collect data

        % encoding
        allEncMrv{iCell, 1}     = cellEncMrv;
        allEncMean{iCell, 1}    = cellEncMean;
        allEncRank{iCell, 1}    = encRanks(iCell, 1);

        % recall
        allRecMrv{iCell, 1}     = cellRecMrv;
        allRecMean{iCell, 1}    = cellRecMean;
        allRecRank{iCell, 1}    = recRanks(iCell, 1);

        % rank of Watson-Williams test
        allWwRank{iCell, 1}     = wwRank;

        %% create figure with phase angle plot for encoding and recall

        % create figure
        close all;
        angleplot      = figure('units', 'centimeters', 'position', [2, 2, 16, 19]);
        colororder({'k', 'b'});
        ax1            = axes('units', 'centimeters', 'position', [0, 0, 16, 18]);
        thisCellRegion = strsplit(parforRes.brainRegion{iCell, 1}, '_');
        title(strcat(conditions{iCond}, '{ Unit Phases - }', thisCellRegion{1, 1}, '{ }', thisCellRegion{1, 2}), 'FontSize', 15);
        set(ax1, 'Visible', 'off');
        set(findall(ax1, 'type', 'text'), 'visible', 'on');

        %% surrogate distribution plot

        % plot distribution
        perm         = axes('units', 'centimeters', 'position', [10, 13.5, 4, 3]);
        permDistr    = histogram(surWwF, 'FaceColor', [0.5, 0.5, 0.5]);
        hold on;
        empVal       = xline(wwF, 'r');
        title({strcat('{Permutation test (P = }', num2str(wwP, 3), ')')});
        xlabel('f-value Watson-Williams test');
        ylabel('Number of surrogates');
        box off;
        set(gca, 'tickDir', 'out');

        %% spike density plot

        % get data
        cellDir      = TG_GetChanDir_20210812(parforRes.paths.SUData, parforRes.subjects{parforRes.idx{iCell, 1}(1, 1)}, strcat('session_', num2str(parforRes.idx{iCell, 1}(1, 2))));
        cc           = load(fullfile(cellDir(parforRes.idx{iCell, 1}(1, 3)).folder, cellDir(parforRes.idx{iCell, 1}(1, 3)).name, 'times_datacut.mat'));
        cc.thisSpike = cc.spikes(cc.cluster_class(:, 1) == parforRes.idx{iCell, 1}(1, 4), :);

        % create axis
        spk          = axes('units', 'centimeters', 'position', [2, 13.5, 4, 3]);
        hold on;

        % spike density plot
        spikeTime    = ((0:size(cc.thisSpike, 2) - 1) ./ cc.par.sr) .* 1000; % [msec]
        outPlot      = TG_DensityPlot(spikeTime, cc.thisSpike);
        hold off;
        set(gca, ...
            'xlim', round([min(spikeTime), max(spikeTime)]), ...
            'xtick', round([min(spikeTime), max(spikeTime)]), ...
            'ylim', [outPlot.lbound, outPlot.ubound], ...
            'ytick', [outPlot.lbound, outPlot.ubound], ...
            'ticklength', [0, 0]);
        colormap(spk, 'parula');
        title('Spike density plot');
        box off;

        % enhance axes
        xlabel('ms', ...
            'FontUnits', 'centimeters', 'units', 'normalized', 'position', [0.5, -0.1, 0]);
        ylabel('\muV', ...
            'FontUnits', 'centimeters', 'units', 'normalized', 'position', [-0.12, 0.5, 0]);

        %% plot polar histogram for encoding
        enc1     = axes('units', 'centimeters', 'position', [2, 6.5, 4, 4]);
        encHistc = histcounts(cellEncPhase, parforRes.param.polarHistEdges);
        encHistg = polarhistogram('BinEdges', parforRes.param.polarHistEdges, 'BinCounts', encHistc, 'FaceColor', 'k');
        hold on;
        polarplot([cellEncMean, cellEncMean], [0, max(encHistc)], 'r', 'LineWidth', 2);
        set(gca, 'ThetaTickLabel', parforRes.param.polarHistLabels);
        rlim([0, max(encHistc)]);
        title({strcat('{Encoding (}', '{n = }', num2str(size(cellEncPhase, 1)), '{ spikes, }'), strcat('{P = }', num2str(encP, 3), ')')});

        %% plot polar histogram for recall
        rec1     = axes('units', 'centimeters', 'position', [10, 6.5, 4, 4]);
        recHistc = histcounts(cellRecPhase, parforRes.param.polarHistEdges);
        recHistg = polarhistogram('BinEdges', parforRes.param.polarHistEdges, 'BinCounts', recHistc, 'FaceColor', 'k');
        hold on;
        polarplot([cellRecMean, cellRecMean], [0, max(recHistc)], 'r', 'LineWidth', 2);
        set(gca, 'ThetaTickLabel', parforRes.param.polarHistLabels);
        rlim([0, max(recHistc)]);
        title({strcat('{Recall (}', '{n = }', num2str(size(cellRecPhase, 1)), '{ spikes, }'), strcat('{P = }', num2str(recP, 3), ')')});

        %% spike-triggered average for encoding

        % y-axis limits
        yAxisLim   = max(abs([encSTAFilt.trial{1, 1}, recSTAFilt.trial{1, 1}])) * 1.2;
        if yAxisLim < 20 && yAxisLim >= 10
            yAxisLim = 20;
        elseif yAxisLim < 10
            yAxisLim = 10;
        end

        % plot
        enc2       = axes('units', 'centimeters', 'position', [1, 1, 6, 4]);
        encSTAplot = plot(encSTAFilt.time{1, 1}, encSTAFilt.trial{1, 1}, 'Color', 'k');
        hold on;
        xline(0, 'LineStyle', ':', 'LineWidth', 1.5);
        xlim([-0.5, 0.5]);
        ylim([-yAxisLim, yAxisLim]);
        % ylim([-10, 10]);
        box off;
        set(gca, 'tickDir', 'out');
        title('Spike-triggered average');
        xlabel('Time [s]');
        ylabel('Voltage [µV]');

        %% spike-triggered average for recall
        rec2       = axes('units', 'centimeters', 'position', [9, 1, 6, 4]);
        recSTAplot = plot(recSTAFilt.time{1, 1}, recSTAFilt.trial{1, 1}, 'Color', 'k');
        hold on;
        xline(0, 'LineStyle', ':', 'LineWidth', 1.5);
        xlim([-0.5, 0.5]);
        ylim([-yAxisLim, yAxisLim]);
        % ylim([-10, 10]);
        box off;
        set(gca, 'tickDir', 'out');
        title('Spike-triggered average');
        xlabel('Time [s]');
        ylabel('Voltage [µV]');

        %% save figure
        unitName    = strjoin(string(parforRes.idx{iCell, 1}), '_');
        figureName  = strcat(unitName , '_', conditions{iCond}, '_PhaseAngleFigure');
        print(angleplot, fullfile(paths.save, 'UnitFigures', figureName), '-djpeg');

        %% save examples as svg

        % array with examples
        examples = [4, 0, 1, 1; ... % phase locking
            5, 0, 4, 3; ...
            11, 0, 3, 1; ...
            5, 0, 32, 2; ... % phase shifts
            5, 0, 46, 2; ...
            10, 1, 34, 2; ...
            10, 1, 40, 4];

        % save example as svg
        if ismember(parforRes.idx{iCell, 1}, examples, 'row') && ...
                strcmp(conditions{iCond}, 'both')
            set(angleplot, 'renderer', 'painters');
            saveas(angleplot, fullfile(paths.save, 'UnitFigures', strcat(figureName, '.svg')), 'svg');
        end
    end

    %% collect results

    % cells to include in results structure
    condResCell   = cat(2, parforRes.idx(:), parforRes.brainRegion(:), ...
        allEncMrv, allEncMean, allEncRank, allRecMrv, ...
        allRecMean, allRecRank, allWwRank);

    % fieldnames
    condResFields = {'idx', 'brainRegion', ...
        'allEncMrv', 'allEncMean', 'allEncRank', 'allRecMrv', ...
        'allRecMean', 'allRecRank', 'allWwRank'};

    % create results structure
    condRes       = cell2struct(condResCell, condResFields, 2);

    % collect information for each condition
    if strcmp(conditions{iCond}, 'both')
        bothRes = condRes;
    elseif strcmp(conditions{iCond}, 'succ')
        succRes = condRes;
    elseif strcmp(conditions{iCond}, 'fail')
        failRes = condRes;
    end
end

% shut down parallel pool
delete(gcp);

% exclude units with too few spikes
bothRes             = bothRes(~bExclude);
succRes             = succRes(~bExclude);
failRes             = failRes(~bExclude);
phaseRes            = r.allRes(~bExclude);

% remove big fields from phaseRes structure
phaseRes            = rmfield(phaseRes, { ...
    'allComplex', 'allSpikeSlope', 'allEncRecIdx', ...
    'allPowerIdx', 'allSlopeIdx', 'allThetaIdx', 'allTime', ...
    'encStartStop', 'recStartStop', ...
    'encBothSTA', 'encSuccSTA', 'encFailSTA', ...
    'recBothSTA', 'recSuccSTA', 'recFailSTA'});

%% detect object cells

% preallocate
objectCellRank  = nan(size(phaseRes));
bObjectCell     = false(size(phaseRes));

% loop through cells
for iCell = 1:size(phaseRes, 1)

    % rank of empirical t-value in surrogate values
    objectCellRank(iCell, 1)    = TG_PermTest1D_2PS_20230727(phaseRes(iCell).encFr, phaseRes(iCell).encBlFr, parforRes.param.nSur, s.randSeedNum);

    % logical whether object cell or not
    if objectCellRank(iCell, 1) > 0.95
        bObjectCell(iCell, 1)   = true;
    end
end

% remove old object cell fields from phase results structure
if any(contains(fieldnames(phaseRes), 'bObjectCell'))
    phaseRes    = rmfield(phaseRes, 'bObjectCell');
end
if any(contains(fieldnames(phaseRes), 'bObjectCell'))
    phaseRes    = rmfield(phaseRes, 'bObjectCell');
end

% create phase results structure with object cell index
phaseResCell    = struct2cell(phaseRes);
phaseResCellNew = cat(1, phaseResCell(1:2, :), num2cell(bObjectCell'), num2cell(objectCellRank'), phaseResCell(3:end, :));
oldFieldnames   = fieldnames(phaseRes);
newFieldnames   = cat(1, oldFieldnames(1:2), 'bObjectCell', 'objectCellRank', oldFieldnames(3:end));
phaseRes        = cell2struct(phaseResCellNew, newFieldnames);

%% save additional results
save(fullfile(paths.save, 'additionalResults'), 'bothRes', 'succRes', 'failRes', 'phaseRes', 'bExclude', '-v7.3');

%% save version without large fields
phaseRes = rmfield(phaseRes, {'labelEnglish', 'labelGerman'});
save(fullfile(paths.save, 'additionalResultsSmall'), 'bothRes', 'succRes', 'failRes', 'phaseRes', 'bExclude', '-v7.3');
