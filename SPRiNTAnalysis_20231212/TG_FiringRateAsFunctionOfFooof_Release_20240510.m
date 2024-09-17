%==========================================================================
% This script is for analyzing the relationship between fooof slope/offset
% and neuronal firing rates.
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % results folder
paths.artifact          = 'D:\TreasureHunt\ArtifactDetection_20230907'; % IED/artifact detection
paths.phaseData         = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.sprint            = 'D:\TreasureHunt\MicroSPRiNT_20230713'; % SPRiNT results
paths.SUData            = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.save              = 'D:\TreasureHunt\SPRiNTAnalysis_20231212'; % save folder

% parameters
param                   = [];
param.nSur              = 10001; % 10001
param.surHistEdges      = linspace(-10, 10, 201);
param.filterBand        = [1, 10];

% set random seed
randSeed    = 444;
rng(randSeed, 'twister');
randSeedNum = randi(100000, param.nSur, 1);

% own functions
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieldtrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

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

%% load phase results data
results             = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));
% allResults          = load(fullfile(paths.phaseRes, 'allResults.mat'));

% get data from all results
allIdx              = cat(1, results.phaseRes.idx); % unit index

% preallocate
topSlopeFr          = nan(size(allIdx, 1), 1);
lowSlopeFr          = nan(size(allIdx, 1), 1);
thetaFr             = nan(size(allIdx, 1), 1);
noThetaFr           = nan(size(allIdx, 1), 1);

%% loop through units
parfor iCell = 1:size(allIdx, 1)
    
    % display progress
    disp(iCell);
    
    %% get cell information

    % get index of this unit
    cellIdx             = allIdx(iCell, :);
    
    % get subject and session name
    cellSub             = subjects{cellIdx(1, 1)};
    cellSess            = strcat('session_', num2str(cellIdx(1, 2)));
    
    % get wire name
    mwChanDir           = TG_GetChanDir_20210812(paths.SUData, cellSub, cellSess);  % unit data
    cellWire            = mwChanDir(cellIdx(1, 3)).name;

    % get unit number
    cellUnit            = cellIdx(1, 4);
    
    %% load data

    % load sprint result of this wire
    sprintData          = load(fullfile(paths.sprint, cellSub, cellSess, cellWire, 'datacutSprint.mat'));

    % load data for artifact removal
    artifactData        = load(fullfile(paths.artifact, cellSub, cellSess, cellWire, 'detectedArtifacts.mat'), 'bArtifact');
    
    % load sampling rate of artifact index
    phaseData           = load(fullfile(paths.phaseData, cellSub, cellSess, cellWire, 'datacutFt2000HzBP_generalized_1_10.mat'));

    % load unit data
    t                   = load(fullfile(paths.SUData, cellSub, cellSess, cellWire, 'times_datacut.mat'));

    %% slope index for this wire

    % get SPRiNT time
    sprintTimeOrig      = cat(2, sprintData.SPRiNT.channel.aperiodics.time);

    % SPRiNT time in samples
    sprintSamples       = sprintTimeOrig * phaseData.fsample;

    % get SPRiNT slope
    sprintSlopeOrig     = cat(2, sprintData.SPRiNT.channel.aperiodics.exponent);

    % overall time coverage of SPRiNT windows
    sprintWindow        = sprintData.SPRiNT.options.WinLength + (sprintData.SPRiNT.options.WinAverage - 1) * (sprintData.SPRiNT.options.WinLength * (1 - (sprintData.SPRiNT.options.WinOverlap / 100)));

    % samples of overall SPRiNT windows
    sprintWindowsSmp    = sprintSamples' + ((-(sprintWindow / 2 * phaseData.fsample) + 1):1:((sprintWindow / 2 * phaseData.fsample)));

    % find artifacts
    sprintBArtifact     = any(artifactData.bArtifact(sprintWindowsSmp), 2)';

    % upsample slope and artifact index
    sprintTimeWindow    = mode(diff(sprintTimeOrig));
    sprintFsample       = 1 / sprintTimeWindow;
    upsamplingRatio     = phaseData.fsample / sprintFsample; % upsampling factor
    sprintSlope         = repelem(sprintSlopeOrig, upsamplingRatio);
    bArtifact           = repelem(sprintBArtifact, upsamplingRatio);

    % add NaNs to SPRiNT data to start and end borders for padding
    missingSprint       = size(artifactData.bArtifact, 2) - size(sprintSlope, 2);
    sprintSlope         = cat(2, repelem(NaN, missingSprint / 2), sprintSlope, repelem(NaN, missingSprint / 2));
    bArtifact           = cat(2, repelem(true, missingSprint / 2), bArtifact, repelem(true, missingSprint / 2));

    % get slope quantiles
    slopeMedian         = median(sprintSlopeOrig(~sprintBArtifact));

    % create slope index
    slopeIdx            = double(sprintSlope > slopeMedian);

    % do not include padded start and end of slope index
    slopeIdx(1, 1:missingSprint / 2)                = NaN;
    slopeIdx(1, end - missingSprint / 2 + 1: end)   = NaN;

    %% theta index for this wire

    % extract sprint times with theta peaks
    centerFrequency                 = cat(1, sprintData.SPRiNT.channel.peaks.center_frequency); % center frequency of peaks

     % extract only peaks in theta band
    thetaPeakIdx                    = centerFrequency >= param.filterBand(1, 1) & centerFrequency <= param.filterBand(1, 2);

    % get time of all peaks
    peakTime                        = cat(1, sprintData.SPRiNT.channel.peaks.time);

    % get only times of theta peaks
    timesWithTheta                  = unique(peakTime(thetaPeakIdx));

    % create index for theta oscillations in sprint signal
    [~, thetaIdxOrig]               = min(abs(sprintTimeOrig - timesWithTheta), [], 2);
    bThetaOrig                      = zeros(size(sprintTimeOrig));
    bThetaOrig(1, thetaIdxOrig)     = ones;
    thetaIdx                        = repelem(bThetaOrig, upsamplingRatio);

    % add NaNs to SPRiNT data to start and end borders for padding
    thetaIdx                        = cat(2, repelem(NaN, missingSprint / 2), thetaIdx, repelem(NaN, missingSprint / 2));

    %% unit data

    % get spike times of this cluster
    thisClusterOrig         = t.cluster_class(t.cluster_class(:, 1) == cellUnit, 2) / 1000; % convert to seconds
    thisClusSamplesOrig     = round(thisClusterOrig * phaseData.fsample);

    % exclude artifacts
    thisClusArtifactIdx     = bArtifact(thisClusSamplesOrig)';
    thisCluster             = thisClusterOrig(~thisClusArtifactIdx);
    thisClusSamples         = thisClusSamplesOrig(~thisClusArtifactIdx);

    % get slope and theta index for each spike
    thisClusSlopeIdx        = slopeIdx(thisClusSamples)';
    thisClusThetaIdx        = thetaIdx(thisClusSamples)';

    % exclude spikes during first and last time of signal, because there is
    % no sprint result for them
    sprintStart             = sprintTimeOrig(1) - sprintTimeWindow / 2;
    sprintEnd               = sprintTimeOrig(end) + sprintTimeWindow / 2;
    spikeExclusionIdx       = thisCluster > sprintStart & thisCluster < sprintEnd;
    spikeSlopeIdx           = thisClusSlopeIdx(spikeExclusionIdx);
    spikeThetaIdx           = thisClusThetaIdx(spikeExclusionIdx);
    
    %% calculate firing rates for different indices

    % calculate overall time for different slope indices
    durTopSlope             = sum(slopeIdx == 1 & bArtifact == 0) / phaseData.fsample;
    durLowSlope             = sum(slopeIdx == 0 & bArtifact == 0) / phaseData.fsample;

    % calculate overall time for different theta indices
    durTheta                = sum(thetaIdx == 1 & bArtifact == 0) / phaseData.fsample;
    durNoTheta              = sum(thetaIdx == 0 & bArtifact == 0) / phaseData.fsample;

    % calculate firing rates
    topSlopeFr(iCell, 1)    = sum(spikeSlopeIdx == 1) / durTopSlope;
    lowSlopeFr(iCell, 1)    = sum(spikeSlopeIdx == 0) / durLowSlope;

    % calculate firing rates
    thetaFr(iCell, 1)       = sum(spikeThetaIdx == 1) / durTheta;
    noThetaFr(iCell, 1)     = sum(spikeThetaIdx == 0) / durNoTheta;
end

%% histogram to compare firing rates

% slope
medianSlopeFr   = median([topSlopeFr, lowSlopeFr], 'omitnan');
slopeFrFig    	= figure('units', 'centimeters', 'Position', [10, 10, 6, 12]);
parallelcoords([topSlopeFr, lowSlopeFr], 'Color', [0.4, 0.4, 0.4, 0.6], ...
    'LineStyle', '-', 'LineWidth', 0.1);
hold on;
parallelcoords(medianSlopeFr, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 1);
set(gca, 'yscale', 'log');
set(gca, 'XTickLabel',{'high', 'low'});
xlim([0.75, 2.25]);
ylim([0.1, 30]);
set(gca, 'tickDir', 'out');
set(slopeFrFig, 'renderer', 'painters');
saveas(slopeFrFig, fullfile(paths.save, 'ALL_SlopeMedianFiringRate.svg'), 'svg');

% theta
medianThetaFr   = median([thetaFr, noThetaFr], 'omitnan');
thetaFrFig    	= figure('units', 'centimeters', 'Position', [10, 10, 6, 12]);
parallelcoords([thetaFr, noThetaFr], 'Color', [0.4, 0.4, 0.4, 0.6], ...
    'LineStyle', '-', 'LineWidth', 0.1);
hold on;
parallelcoords(medianThetaFr, 'Color', 'blue', 'LineStyle', '-', 'LineWidth', 1);
set(gca, 'yscale', 'log');
set(gca, 'XTickLabel',{'theta', 'no theta'});
xlim([0.75, 2.25]);
ylim([0.01, 30]);
set(gca, 'tickDir', 'out');
set(thetaFrFig, 'renderer', 'painters');
saveas(thetaFrFig, fullfile(paths.save, 'ALL_ThetaMedianFiringRate.svg'), 'svg');

%% surrogate tests
[slopePval, slopeT, surSlopeT] = TG_PermTest1D_2PS_20230727(topSlopeFr, lowSlopeFr, param.nSur, randSeedNum);
[thetaPval, thetaT, surThetaT] = TG_PermTest1D_2PS_20230727(thetaFr, noThetaFr, param.nSur, randSeedNum);

%% plot histogram of surrogate and empirical t-value

% slope
surSlopeFigure     = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
surSlopeTCounts    = histcounts(surSlopeT, param.surHistEdges);
surSlopeHist       = histogram('BinCounts', surSlopeTCounts, 'BinEdges', param.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(slopeT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-10, 10]);
xticks([-10, 0, 10]);
ylim([0, 500]);
yticks([0, 500]);
set(gca,'TickDir','out');
box off;
saveas(surSlopeFigure, fullfile(paths.save, 'ALL_SlopeSurDistrFiringRate.svg'));

% theta
surThetaFigure     = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
surThetaTCounts    = histcounts(surThetaT, param.surHistEdges);
surThetaHist       = histogram('BinCounts', surThetaTCounts, 'BinEdges', param.surHistEdges, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
xline(thetaT, 'color', [0, 0, 0], 'LineWidth', 1);
xlim([-10, 10]);
xticks([-10, 0, 10]);
ylim([0, 500]);
yticks([0, 500]);
set(gca,'TickDir','out');
box off;
saveas(surThetaFigure, fullfile(paths.save, 'ALL_ThetaSurDistrFiringRate.svg'));
