%==========================================================================
% This script creates examples of LFP data, corresponding spike trains, and
% corresponding SPRiNT data
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;
rng(444);

% paths
paths.behlog    = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.SUData    = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.LFPData   = 'D:\TreasureHunt\MicroDownsampled_20210910'; % LFP data
paths.phaseData = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.sprint    = 'D:\TreasureHunt\MicroSPRiNT_20230713'; % SPRiNT data
paths.subInfo   = 'D:\TreasureHunt\SubjectInformation_20210111'; % subject information
paths.save      = 'D:\TreasureHunt\SPRiNTAnalysis_20231212\SprintExamples'; % save folder
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% own functions
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% parameters
param            	= [];
param.c4aName       = {'Cluster4Analysis_LK_20200416_113912', 'Cluster4Analysis_TG_20210315_020120'};
param.timeBorders  	= 2; % time in seconds before and after trial to show
param.filterBand    = [1, 10]; % theta frequency band

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

%% get data

% get example subject, session, wire and cluster
% iSub = 11; iSess = 1; iWire = 3; iClus = 1; iSeg = 90; % phase locking
iSub = 5; iSess = 1; iWire = 2; iClus = 1; iSeg = 95; % slopes

% get sessions
sessions            = dir(fullfile(paths.SUData, subjects{iSub}, 'session*'));

% load behavioral segment info file
segmentInfo         = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'segmentInfo_20230208.mat'));

% get available data
mwChanDir           = TG_GetChanDir_20210812(paths.SUData, subjects{iSub}, sessions(iSess).name);  % unit data

% path for this channel's unit data
mwChanPath          = fullfile(paths.SUData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name);

% load LFP data
LFPdata             = load(fullfile(paths.LFPData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'datacutFt2000Hz.mat'));

% load phase data
phaseData           = load(fullfile(paths.phaseData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'datacutFt2000HzBP_generalized_1_10.mat'));

% load sprint data
sprintData          = load(fullfile(paths.sprint, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'datacutSprint.mat'));

%% extract relevant data

% SPRiNT
sprintFreqs         = sprintData.Freqs';
sprintLogFreqs      = log10(sprintFreqs);
sprintSlope         = cat(2, sprintData.SPRiNT.channel.aperiodics.exponent); % exponent (negative slope)
sprintTime          = cat(1, sprintData.SPRiNT.channel.data.time)'; % time in seconds
sprintPow           = cat(1, sprintData.SPRiNT.channel.data.power_spectrum)';
sprintFooof         = cat(1, sprintData.SPRiNT.channel.data.fooofed_spectrum)';
sprintAperiodic     = cat(1, sprintData.SPRiNT.channel.data.ap_fit)';
sprintPeak          = cat(1, sprintData.SPRiNT.channel.data.peak_fit)';

% slope index for this wire
slopeMedian         = median(sprintSlope); % get slope quantiles
slopeIdx            = sprintSlope > slopeMedian; % create slope index

% extract sprint times with theta peaks
centerFrequency     = cat(1, sprintData.SPRiNT.channel.peaks.center_frequency); % center frequency of peak
relevantPeakIdx     = centerFrequency >= param.filterBand(1, 1) & centerFrequency <= param.filterBand(1, 2); % extract only peaks in the theta band
peakTime            = cat(1, sprintData.SPRiNT.channel.peaks.time); % get time of all peaks
timesWithPeak       = unique(peakTime(relevantPeakIdx)); % get only times with at least one peak present

% create index for theta oscillations in sprint signal
[~, sprintTimeDiff]     	= min(abs(sprintTime - timesWithPeak), [], 2);
thetaIdx                   	= false(size(sprintTime));
thetaIdx(1, sprintTimeDiff)	= true;

% load wave-clus output (for spike-depiction etc.)
t = load(fullfile(mwChanPath, 'times_datacut.mat'));

% load unit classification file
try
    c4a = load(fullfile(mwChanPath, param.c4aName{1, 1}));
catch
    c4a = load(fullfile(mwChanPath, param.c4aName{1, 2}));
end
c4a = c4a.Cluster4Analysis;

% choose cluster with most spikes
allClusters         = t.cluster_class(:, 1);
relClusters         = allClusters(allClusters ~= 0);

% get data for only selected cluster
thisCluster         = t.cluster_class(t.cluster_class(:, 1) == iClus, 2) / 1000; % cluster-number, convert to seconds

% info about encoding segments
timeInfo            = segmentInfo.encoding / segmentInfo.fsample; % convert to seconds

% add time before and after trial
borderInfo          = timeInfo;
borderInfo(:, 1)    = timeInfo(:, 1) - param.timeBorders;
borderInfo(:, 2)    = timeInfo(:, 2) + param.timeBorders;

% cut data
cfg                 = [];
cfg.trl             = round(borderInfo * LFPdata.fsample);
cfg.trl(:, 3)       = -param.timeBorders * LFPdata.fsample; % adjust to sampling rate of data
allSegLFPData       = ft_redefinetrial(cfg, LFPdata);
allSegPhaseData     = ft_redefinetrial(cfg, phaseData);

%% get spikes and SPRiNT results of this segment

% get spikes of this segment
segSpikesIdx        = find(thisCluster > borderInfo(iSeg, 1) & thisCluster < borderInfo(iSeg, 2));
segSpikesTime       = thisCluster(segSpikesIdx);
spikeTimes          = segSpikesTime - timeInfo(iSeg, 1);

% get SPRiNT results of this segment
sprintSegIdx        = sprintTime > borderInfo(iSeg, 1) & sprintTime < borderInfo(iSeg, 2);
sprintSegTime       = sprintTime(sprintSegIdx) - timeInfo(iSeg, 1);
sprintSegPow       	= sprintPow(:, sprintSegIdx);
sprintSegFooof      = sprintFooof(:, sprintSegIdx);
sprintSegAperiodic  = sprintAperiodic(:, sprintSegIdx);
sprintSegSlope      = sprintSlope(:, sprintSegIdx);
sprintSegSlopeIdx   = slopeIdx(:, sprintSegIdx);
sprintSegThetaIdx   = thetaIdx(:, sprintSegIdx);

%% create figure
close all;
approach            = figure('units', 'normalized', 'outerposition', [0, 0, 1, 1]);
set(gcf, 'color', 'w');

% box settings
scaleFactor         = 100; % rescaling log-power results to make them fit to axis
boxBottom           = -10 * scaleFactor; % position of box bottom
boxHeight           = 5 * scaleFactor;
boxWidth            = mean(diff(sprintSegTime));

% virtual box axes
xAxisLimits         = [sprintData.Freqs(1, 1), sprintData.Freqs(1, end)];
yAxisLimits         = [0.1, ((-boxHeight - boxBottom) - scaleFactor) / scaleFactor];

%% mark encoding segment
rectangle('Position', [0, boxBottom, 1.5, 1600], 'FaceColor', [0.9, 0.9, 0.9], 'EdgeColor', [0.9, 0.9, 0.9]);
set(gca, 'Visible', 'off');
hold on;

%% LFP plot
lfpPlot             = plot(allSegLFPData.time{1, iSeg}, allSegLFPData.trial{1, iSeg}, 'Color', 'k', 'LineWidth', 1);

% scale bar
xbar = line([3, 3], [-250, -150], 'Color', 'k', 'LineWidth', 3);
ybar = line([1.8, 2.8], [-270, -270], 'Color', 'k', 'LineWidth', 3);
text(2.2, -300, '1 second');
text(3.1, -240, '100 \muV', 'rotation', 90);

%% plot spikes
if size(spikeTimes, 1) > 0
    spikePlot           = line([spikeTimes'; spikeTimes'],[-450, -400], 'Color', 'k', 'LineWidth', 1);
end

%% plot SPRiNT results (underneath the LFP and spiketrain)

% create virtual x-axes for sprint segments
virtXAxes           = sprintSegTime + (sprintLogFreqs - (sprintLogFreqs(end) / 2)) * (mean(diff(sprintSegTime)) / sprintLogFreqs(end)) * 0.8;

% plot box and slope level
for iSprint = 1:sum(sprintSegIdx)

    % plot box
    rectangle('Position', [sprintSegTime(1, iSprint) - (boxWidth * 0.5), boxBottom, boxWidth, boxHeight], 'FaceColor', [1, 1, 1], 'EdgeColor', [0.5, 0.5, 0.5]);

    % label with slope index
    if sprintSegSlopeIdx(1, iSprint) == 1
        text(virtXAxes(3, iSprint) + 0.02, boxBottom + 75, 'high');
    elseif sprintSegSlopeIdx(1, iSprint) == 0
        text(virtXAxes(4, iSprint), boxBottom + 75, 'low');
    end

    % label with theta index
    if sprintSegThetaIdx(1, iSprint) == 1
        text(virtXAxes(3, iSprint), boxBottom + 25, 'theta');
    elseif sprintSegThetaIdx(1, iSprint) == 0
        text(virtXAxes(2, iSprint), boxBottom + 25, 'no theta');
    end
end

% plot sprint results
powPlot             = plot(virtXAxes, log10(sprintSegPow) * scaleFactor + boxBottom + scaleFactor, 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);
fooofPlot           = plot(virtXAxes, log10(sprintSegFooof) * scaleFactor + boxBottom + scaleFactor, 'Color', 'r', 'LineWidth', 1);
aperiodicPlot       = plot(virtXAxes, log10(sprintSegAperiodic) * scaleFactor + boxBottom + scaleFactor, 'Color', 'b', 'LineStyle', ':', 'LineWidth', 2);

% plot one example
figure;
powPlot1            = plot(sprintFreqs, sprintSegPow(:, 1), 'Color', [0.4, 0.4, 0.4], 'LineWidth', 2);
hold on;
fooofPlot1          = plot(sprintFreqs, sprintSegFooof(:, 1), 'Color', 'r', 'LineWidth', 1);
aperiodicPlot1      = plot(sprintFreqs, sprintSegAperiodic(:, 1), 'Color', 'b', 'LineStyle', ':', 'LineWidth', 2);
apSlope             = mean(diff(log10(sprintSegAperiodic(:, 1))) ./ diff(log10(sprintFreqs)));
sprintSegShift      = sprintSegFooof(:, 1) .* (sprintSegAperiodic(:, 1) ./ sprintSegAperiodic(end, 1));
% shiftPlot           = plot(sprintFreqs, sprintSegShift, 'Color', 'r', 'LineStyle', ':', 'LineWidth', 2);
ylim([0.1, 10000]);
set(gca, 'XScale', 'log', 'YScale', 'log', 'TickDir', 'out', 'visible', 'on');
% box off;

% print figure
segIdx = regexprep(num2str([iSub, iSess, iWire, iClus, iSeg]), '\s+', '_');
print(approach, fullfile(paths.save, strcat(segIdx, '_SprintExample')), '-djpeg', '-r300');
