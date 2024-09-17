%==========================================================================
% This script creates examples of 1-10 Hz filtered data and
% corresponding spike trains
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.behlog    = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.SUData    = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.filtData  = 'D:\TreasureHunt\MicroFiltered_20210930'; % filtered data
paths.LFPData   = 'D:\TreasureHunt\MicroDownsampled_20210910'; % LFP data
paths.phaseData = 'D:\TreasureHunt\MicroPhases_20210930'; % phase data
paths.subInfo   = 'D:\TreasureHunt\SubjectInformation_20210111'; % subject information
paths.save      = 'D:\TreasureHunt\DataExamples_20221230'; % save folder

% add functions
addpath(genpath('D:\External\Functions'));
addpath(genpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% parameters
param                 = [];
param.filterBand      = [1, 10]; % filter ban for phase detection (for loading of correct file)
param.phaseMethod     = 'generalized';
param.c4aName         = {'Cluster4Analysis_LK_20200416_113912', 'Cluster4Analysis_TG_20210315_020120'};
param.timeBorders     = 2; % time in seconds before and after trial to show
param.polarHistLabels = {'0'; '90'; '\pm180'; '-90'};

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

% get example data
exampleSegments = ...
    [1, 1, 1, 1, 33; ...
    1, 3, 32, 2, 69; ...
    2, 1, 61, 1, 29; ...
    2, 2, 43, 1, 59; ...
    3, 2, 12, 1, 38; ...
    3, 3, 11, 1, 28; ...
    4, 1, 2, 1, 75; ...
    4, 2, 7, 1, 58; ...
    5, 1, 7, 1, 49; ...
    5, 2, 4, 2, 83; ...
    6, 1, 24, 3, 84; ...
    6, 2, 24, 1, 30; ...
    7, 1, 3, 1, 58; ...
    8, 1, 33, 1, 70; ...
    8, 2, 37, 1, 43; ...
    9, 1, 46, 1, 41; ...
    10, 1, 42, 2, 17; ...
    10, 2, 46, 1, 91; ...
    11, 1, 5, 1, 89; ...
    11, 1, 3, 1, 90; ... % example in main figure
    12, 1, 21, 1, 95; ...
    13, 1, 4, 1, 82; ...
    14, 1, 27, 1, 75; ...
    15, 1, 25, 1, 43; ...
    16, 1, 32, 1, 13; ...
    17, 1, 15, 1, 54; ...
    18, 1, 2, 1, 92; ...
    18, 2, 1, 1, 99];

% loop through examples
for iEx = 1:size(exampleSegments, 1)
    
    % define subject, session, wire, cluster and segment
    iSub            = exampleSegments(iEx, 1);
    iSess           = exampleSegments(iEx, 2);
    iWire           = exampleSegments(iEx, 3);
    iClus           = exampleSegments(iEx, 4);
    iSeg            = exampleSegments(iEx, 5);
    
    % get sessions
    sessions        = dir(fullfile(paths.SUData, subjects{iSub}, 'session*'));
    
    % load behavioral segment info file
    segmentInfo     = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'segmentInfo_20230208.mat'));
    
    % get available data
    mwChanDir       = TG_GetChanDir_20210812(paths.SUData, subjects{iSub}, sessions(iSess).name);  % unit data
    LFPChanDir      = TG_GetChanDir_20210812(paths.LFPData, subjects{iSub}, sessions(iSess).name); % LFP data
    FiltChanDir     = TG_GetChanDir_20210812(paths.filtData, subjects{iSub}, sessions(iSess).name); % filtered data
    phaseChanDir    = TG_GetChanDir_20210812(paths.phaseData, subjects{iSub}, sessions(iSess).name); % phase data
    
    % path for this channel's unit data
    mwChanPath      = fullfile(paths.SUData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name);
    
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
    
    % load LFP data
    LFPdata         = load(fullfile(paths.LFPData, subjects{iSub}, sessions(iSess).name, LFPChanDir(iWire).name, 'datacutFt2000Hz.mat'));
    
    % load filtered data
    filtDataName    = strcat('datacutFt2000HzBP_', regexprep(num2str(param.filterBand), '\s+', '_'), '.mat');
    filtData        = load(fullfile(paths.filtData, subjects{iSub}, sessions(iSess).name, LFPChanDir(iWire).name, filtDataName));
    
    % load phase data
    phaseDataName   = strcat('datacutFt2000HzBP_', param.phaseMethod, '_', regexprep(num2str(param.filterBand), '\s+', '_'), '.mat');
    phaseData       = load(fullfile(paths.phaseData, subjects{iSub}, sessions(iSess).name, phaseChanDir(iWire).name, phaseDataName));
    
    % load wave-clus output (for spike-depiction etc.)
    t = load(fullfile(mwChanPath, 'times_datacut.mat'));
    
    % load unit classification file
    try
        c4a = load(fullfile(mwChanPath, param.c4aName{1, 1}));
    catch
        c4a = load(fullfile(mwChanPath, param.c4aName{1, 2}));
    end
    c4a = c4a.Cluster4Analysis;
    
    % get phase data for only selected cluster
    thisCluster     = t.cluster_class(t.cluster_class(:, 1) == iClus, :); % cluster-number, time in ms
    thisClusSamples = round(thisCluster(:, 2) * segmentInfo.fsample / 1000); % convert from ms to phase angle samples
    
    % info about encoding segments
    segmentsOrig    = segmentInfo.encoding;
    
    % add time before and after trial
    segmentsOrig(:, 1) = segmentsOrig(:, 1) - (param.timeBorders * segmentInfo.fsample);
    segmentsOrig(:, 2) = segmentsOrig(:, 2) + (param.timeBorders * segmentInfo.fsample);
    
    % cut data
    cfg             = [];
    cfg.trl         = round(segmentsOrig / (segmentInfo.fsample / LFPdata.fsample)); % adjust sampling rate of segmentInfo to data
    cfg.trl(:, 3)   = -param.timeBorders * LFPdata.fsample; % adjust to sampling rate of data
    allSegLFPData   = ft_redefinetrial(cfg, LFPdata);
    allSegFiltData  = ft_redefinetrial(cfg, filtData);
    allSegPhaseData = ft_redefinetrial(cfg, phaseData);
    
    % get spikes of this segment
    segSpikesIdx     = find(thisClusSamples > segmentsOrig(iSeg, 1) & thisClusSamples < segmentsOrig(iSeg, 2));
    segSpikesSamples = thisClusSamples(segSpikesIdx);
    spikeTimes       = (segSpikesSamples - segmentsOrig(iSeg, 1)) / segmentInfo.fsample - param.timeBorders;
    
    %% create figure
    close all;
    approach = figure('units', 'centimeters', 'position', [5, 5, 20, 4.5]);
    ax1      = axes('units', 'centimeters', 'position', [1, 1, 20, 4.5]);
    set(ax1, 'Visible', 'off');
    hold on;
    set(gcf, 'color', 'w');
    
    % raw data
    exampleData = axes('units', 'centimeters', 'position', [0.5, 0, 20, 3.5], 'FontSize', 24);
    xLimit      = [-1.5, 3.6];
    yLimit      = [min(allSegLFPData.trial{1, iSeg}) - 50, max(allSegLFPData.trial{1, iSeg}) + 20];
    rectangle('Position', [0, yLimit(1, 1) + 1, 1.5, abs(yLimit(1, 1)) + yLimit(1, 2)], 'EdgeColor', 'none', 'FaceColor', [0.9, 0.9, 0.9]);
    hold on;
    plot(allSegLFPData.time{1, iSeg}, allSegLFPData.trial{1, iSeg}, 'Color', 'k', 'LineWidth', 1);
    set(exampleData, 'Visible', 'off', 'TickDir', 'out');
    box off;
    xlim(xLimit);
    ylim(yLimit);
    titleString = strcat({'sub_'}, {num2str(iSub)}, {', '}, {sessions(iSess).name}, {', '}, brainRegionsList{iWire, 1});
    titleString = strrep(titleString, 'sub_', 'subject ');
    titleString = strrep(titleString, '_', ' ');
    title(titleString, 'Interpreter', 'none', 'FontSize', 16, 'FontWeight', 'Normal');
    set(findall(exampleData, 'type', 'text'), 'visible', 'on');
    
    % phase
    cl              = cline(allSegFiltData.time{1, iSeg}, allSegFiltData.trial{1, iSeg}, [], angle(allSegPhaseData.trial{1, iSeg}));
    set(cl, 'LineWidth', 3.5);
    
    % colormap
    hmap(1:256, 1)  = linspace(0, 1, 256);
    hmap(:, [2, 3]) = 0.8; % brightness
    huemap          = hsv2rgb(hmap);
    colormap(huemap);
    
    % scale bar
    rectangle('Position', [3, yLimit(1, 1) + 1, 1, abs(yLimit(1, 1)) + yLimit(1, 2)], 'EdgeColor', 'none', 'FaceColor', [1, 1, 1]);
    xbar = line([3.2, 3.2], [-50, 50], 'Color', 'k', 'LineWidth', 3);
    
    % spiketrain
    spikeTrain      = axes('units', 'centimeters', 'position', [0.5, 0, 20, 0.3]);
    xlim(xLimit);
    set(spikeTrain, 'Visible', 'off');
    if ~isempty(spikeTimes)
        spiketrain  = line([spikeTimes'; spikeTimes'], [0; 1], 'Color', 'k', 'LineWidth', 1);
    end
    rectangle('Position', [3, 0, 1, 1], 'EdgeColor', 'none', 'FaceColor', [1, 1, 1]);

    % segment name
    segIdx = regexprep(num2str([iSub, iSess, iWire, iClus, iSeg]), '\s+', '_');

    % print figure
    if ismember([11, 1, 3, 1, 90], exampleSegments(iEx, :), 'rows')
        set(approach, 'renderer', 'painters');
        saveas(approach, fullfile(paths.save, strcat(segIdx, '_GeneralizedPhaseAndSpikeTrainExample.svg')), 'svg');
    else
        print(approach, fullfile(paths.save, strcat(segIdx, '_GeneralizedPhaseAndSpikeTrainExample')), '-dtiff', '-r300');
    end
end
