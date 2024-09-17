%==========================================================================
% This script analyzes the occurrence and frequency of oscillation peaks
% between 1 and 10 Hz in the power spectrum during the whole experiment.
%
% Tim Guth, 2024
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.sprint    = 'D:\TreasureHunt\MicroSPRiNT_20230713'; % folder with SPRiNT results
paths.artifact  = 'D:\TreasureHunt\ArtifactDetection_20230907'; % IED/artifact detection
paths.phaseData = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.SUData    = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % folder with single unit data
paths.phaseRes	= 'D:\TreasureHunt\PhaseAnalysis_20230921'; % phase analysis folder
paths.save    	= 'D:\TreasureHunt\SPRiNTAnalysis_20231212'; % save folder
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% add functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

%% set configurations
% parameters
param                 = [];
param.percEdges       = linspace(0, 1, 11);
param.binEdges        = linspace(1, 10, 19);
param.binCenters      = (param.binEdges(1 : end - 1) + param.binEdges(2 : end)) / 2;
param.wireEdges       = linspace(1, 10, 10);
param.filterBand      = [1, 10];
param.sprintName      = 'datacutSprint.mat';

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

%% load phase results to get list of included wires
phaseRes            = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

%% wire index
unitIdx             = cat(1, phaseRes.phaseRes.idx);
[wireIdx, selIdx]   = unique(unitIdx(:, 1:3), 'rows'); % all included wires

%% brain region

% brain region for all included units
brainRegIdx          = {phaseRes.phaseRes.brainRegion}';
brainRegIdxSplit     = split(brainRegIdx, '_');

% brain region for all wires
brainRegWireIdx      = brainRegIdxSplit(selIdx, 1);

% names of unique brain regions
uniqueBrainReg       = unique(brainRegIdxSplit(:, 1));

%% loop through wires

% define variables
thetaPerc            = nan(size(wireIdx(:, 1)));
thetaFreq            = cell(size(wireIdx(:, 1)));

% loop through wires
parfor iWire = 1:size(wireIdx, 1)
    
    % report progress
    disp(iWire);
    
    % get path of sprint data
    subName             = subjects{wireIdx(iWire, 1)};
    sessName            = strcat('session_', num2str(wireIdx(iWire, 2)));
    chanDir             = TG_GetChanDir_20210812(paths.SUData, subName, sessName);
    
    % load sprint data
    sprintData          = load(fullfile(paths.sprint, subName, sessName, chanDir(wireIdx(iWire, 3)).name, param.sprintName));
    
    % load data for artifact removal
    artifactData        = load(fullfile(paths.artifact, subName, sessName, chanDir(wireIdx(iWire, 3)).name, 'detectedArtifacts.mat'), 'bArtifact');
    
    % load sampling rate of artifact index
    phaseData           = load(fullfile(paths.phaseData, subName, sessName, chanDir(wireIdx(iWire, 3)).name, 'datacutFt2000HzBP_generalized_1_10.mat'), 'fsample');

    % SPRiNT times
    sprintTime          = cat(2, sprintData.SPRiNT.channel.aperiodics.time);

    % overall time coverage of SPRiNT windows
    sprintWindow        = sprintData.SPRiNT.options.WinLength + (sprintData.SPRiNT.options.WinAverage - 1) * (sprintData.SPRiNT.options.WinLength * (1 - (sprintData.SPRiNT.options.WinOverlap / 100)));

    % SPRiNT times in samples
    sprintSamples       = sprintTime * phaseData.fsample;

    % samples of overall SPRiNT windows
    sprintWindowsSmp    = sprintSamples' + ((-(sprintWindow / 2 * phaseData.fsample) + 1):1:((sprintWindow / 2 * phaseData.fsample)));

    % find artifacts
    sprintBArtifact     = any(artifactData.bArtifact(sprintWindowsSmp), 2)';
    artifactTime        = sprintTime(sprintBArtifact);

    % extract peak times
    peakTime            = cat(1, sprintData.SPRiNT.channel.peaks.time);
    peakFreq            = cat(1, sprintData.SPRiNT.channel.peaks.center_frequency);

    % remove artifacts
    if size(artifactTime, 2) > 0
        peaksToRemove   = min(abs(peakTime - artifactTime), [], 2) < 0.001;
        peakTime        = peakTime(~peaksToRemove);
        peakFreq        = peakFreq(~peaksToRemove);
    end

    % create theta index
    thetaIdx            = (peakFreq >= param.filterBand(1, 1) & peakFreq <= param.filterBand(1, 2));
    thetaTime           = unique(peakTime(thetaIdx));
    thetaPerc(iWire, 1)	= size(thetaTime, 1) / (size(sprintTime, 2) - size(artifactTime, 2));
    thetaFreq{iWire, 1}	= peakFreq(thetaIdx);
end

%% results for different brain regions

% add 'all' to loop
brainRegForLoop = cat(1, {'ALL'}, uniqueBrainReg);

% loop through brain regions
for iReg = 1:size(brainRegForLoop, 1)

    %  selection index for specific regions
    if strcmp(brainRegForLoop{iReg}, 'ALL')
        bSel = true(size(brainRegWireIdx));
    else
        bSel = strcmp(brainRegForLoop{iReg}, brainRegWireIdx);
    end
    
    % only include brain areas with data from at least 5 sessions
    subAndSess          = wireIdx(bSel, 1:2);
    uniqueSess          = unique(subAndSess, 'rows');
    if size(uniqueSess, 1) < 5
        continue;
    end

    % select results
    thisRegThetaPerc    = thetaPerc(bSel);
    thisRegThetaFreq    = thetaFreq(bSel);

    %% print results

    % print histogram of theta percentages
    prctFig     = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    percCounts  = histcounts(thisRegThetaPerc, param.percEdges);
    histogram('BinCounts', percCounts, 'BinEdges', param.percEdges * 100, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Theta presence (%)');
    ylabel('Number of wires');
    title(strcat(brainRegForLoop{iReg}, 32, num2str(sum(bSel)), 32, 'wires'));
    saveas(prctFig, fullfile(paths.save, strcat(brainRegForLoop{iReg}, '_RatioOscillatory.svg')));

    % print histogram of theta frequencies for all wires
    wireFig         = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    wireThetaFreq   = cell2mat(cellfun(@(x) mode(round(x)), thisRegThetaFreq, 'Uni', 0));
    wireThetaCounts = histcounts(wireThetaFreq, param.wireEdges);
    histogram('BinCounts', wireThetaCounts, 'BinEdges', param.wireEdges, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Frequency (Hz)');
    ylabel('Number of wires');
    title(strcat(brainRegForLoop{iReg}, 32, num2str(sum(bSel)), 32, 'wires'));
    saveas(wireFig, fullfile(paths.save, strcat(brainRegForLoop{iReg}, '_FrequenciesWires.svg')));

    % print histogram of theta frequencies for all peaks (normalized)
    binFig      = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
    thetaCounts = cell2mat(cellfun(@(x) histcounts(x, param.binEdges, 'Normalization', 'probability'), thisRegThetaFreq, 'Uni', 0));
    TG_ShadeSEM_20210714(param.binCenters, thetaCounts * 100, 'k', 0.5);
    box off;
    set(gca, 'TickDir', 'out');
    xlabel('Frequency (Hz)');
    ylabel('Probability (%)');
    xticks(0:2:10);
    ylim([0, 15]);
    grid on;
    title(strcat(brainRegForLoop{iReg}, 32, num2str(size(cat(1, thisRegThetaFreq{:}), 1)), 32, 'peaks'));
    saveas(binFig, fullfile(paths.save, strcat(brainRegForLoop{iReg}, '_FrequenciesOfPeaks.svg')));
end
