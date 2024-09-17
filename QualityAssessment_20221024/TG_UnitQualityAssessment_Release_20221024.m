%==========================================================================
% For recording quality assessment, this script calculates the number of
% units recorded on each wire; the ISI refractoriness for each unit;
% the mean firing rate for each unit;and the waveform peak
% signal-to-noise ratio (SNR) for each unit. It also counts the number of
% units per brain region.
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.SUData    = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % unit data from wave_clus
paths.LFPData   = 'D:\TreasureHunt\MicroDownsampled_20210910'; % LFP data
paths.subInfo   = 'D:\TreasureHunt\SubjectInformation_20210111'; % subject information
paths.save      = 'D:\TreasureHunt\QualityAssessment_20221024'; % save folder
saveName        = 'unitQualityAssesment_20221024';

% own functions
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% parameters
param                 = [];
param.c4aName         = {'Cluster4Analysis_LK_20200416_113912', 'Cluster4Analysis_TG_20210315_020120'};

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

%% preallocate
allRes         = cell(length(subjects), 1);

% for sanity checking
allExByWC      = cell(length(subjects), 1); % excluded by wave-clus
allExByVisInsp = cell(length(subjects), 1); % excluded due to visual inspection

%% loop through subjects
parfor iSub = 1:length(subjects)
    
    % get sessions
    sessions         = dir(fullfile(paths.SUData, subjects{iSub}, 'session*'));

    % preallocate
    subRes           = [];
    subWires         = [];

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
        
        % original session index
        sessIdx = split(sessions(iSess).name, '_');
        sessIdx = str2double(sessIdx{2});
        
        % display session information
        fprintf('\n==================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        % get available data
        mwChanDir    = TG_GetChanDir_20210812(paths.SUData, subjects{iSub}, sessions(iSess).name);  % unit data

        % brain region sanity check
        if ~isequal(size(mwChanDir, 1), size(brainRegionsList, 1))
            error('different number of microwire channels than channels in the brain region file');
        end
        
        %% loop through channels
        for iWire = 1:size(mwChanDir, 1)
            
            % print channel name
            fprintf('\tChannel name: %s.\n', mwChanDir(iWire).name);
            
            % path for this channel's unit data
            mwChanPath    = fullfile(paths.SUData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name);

            % load LFP data for fooof analysis
            LFPdata       = load(fullfile(paths.LFPData, subjects{iSub}, sessions(iSess).name, mwChanDir(iWire).name, 'datacutFt2000Hz.mat'));

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
            
            % preallocate number of units per wire
            numUnits = 0;
            
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
                
                % get data for only this cluster
                thisCluster     = t.cluster_class(t.cluster_class(:, 1) == iClus, :); % cluster-number, time in ms
                intSpIv         = diff(thisCluster(:, 2));
                
                % percentage of spikes with inter-spike interval < 3 ms
                percLowISI      = sum(intSpIv < 3) / size(intSpIv, 1);
                
                % mean firing rate of unit
                meanFR          = size(thisCluster, 1) / round(LFPdata.time{1, 1}(end), 3);
                
                % signal to noise ratio
                thisClusSpikes  = t.spikes(t.cluster_class(:, 1) == iClus, :);
                meanSpike       = mean(thisClusSpikes);
                peakAmp         = abs(min(meanSpike));
                stdClus         = mean(std(thisClusSpikes));
                sigToNoise      = peakAmp / stdClus;
                
                % count number of units per wire
                numUnits        = numUnits + 1;
                
                %% collect information for this unit

                % basics
                unitRes                  = [];
                unitRes.idx              = [iSub, sessIdx, iWire, iClus];
                unitRes.brainRegion      = brainRegionsList{iWire, 1};
                unitRes.percLowISI       = percLowISI;
                unitRes.meanFR           = meanFR;
                unitRes.sigToNoise       = sigToNoise;
                
                % collapse across units
                subRes = cat(1, subRes, unitRes);
            end
            
            % collapse number of units per wire across wires
            subWires = cat(1, subWires, numUnits);
        end
    end
    
    % collect results across subjects
    allRes{iSub, 1}         = subRes;
    allWires{iSub, 1}       = subWires;
    allExByWC{iSub, 1}      = exByWC;
    allExByVisInsp{iSub, 1} = exByVisInsp;
end

% unnest result cells
allRes         = cat(1, allRes{:});
allWires       = cat(1, allWires{:});
allWires       = allWires(~(allWires == 0));
allExByWC      = cat(1, allExByWC{:});
allExByVisInsp = cat(1, allExByVisInsp{:});

% save important output
save(fullfile(paths.save, saveName));

%% figures and information

% units per wire (only considering wires with at least one unit)
unitsFig     = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
unitsPerWire = histogram(allWires, 'FaceColor', [0.6, 0.6, 0.6]);
xlabel('Units per wire');
ylabel('Number of wires');
set(gca, 'TickDir', 'out');
box off;
meanUnitsPerWire = mean(allWires);
semUnitsPerWire  = std(allWires) / (sqrt(size(allWires, 1)));

% print figure
print(unitsFig, fullfile(paths.save, 'unitsPerWire'), '-dsvg', '-r300');

% percentage of ISI < 3 ms
isiFig        = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
lowIsiPerUnit = histogram([allRes.percLowISI] * 100, 'FaceColor', [0.6, 0.6, 0.6]);
xlabel('% ISI <3 ms');
ylabel('Number of units');
set(gca, 'TickDir', 'out');
box off;
meanISI = mean([allRes.percLowISI] * 100);
semISI  = std([allRes.percLowISI] * 100) / (sqrt(size([allRes.percLowISI], 2)));
ISI5    = sum([allRes.percLowISI] > 0.05);

% print figure
print(isiFig, fullfile(paths.save, 'lowISIPercentage'), '-dsvg', '-r300');

% mean FR
frFig           = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
meanFRperUnit   = histogram([allRes.meanFR], 'FaceColor', [0.6, 0.6, 0.6]);
xlabel('Mean FR (Hz)');
ylabel('Number of units');
set(gca, 'TickDir', 'out');
box off;
meanFRall       = mean([allRes.meanFR]);
semFRall        = std([allRes.meanFR]) / (sqrt(size([allRes.meanFR], 2)));

% print figure
print(frFig, fullfile(paths.save, 'meanFR'), '-dsvg', '-r300');


% signal-to-noise ratio
snrFig          = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
SNRperUnit      = histogram([allRes.sigToNoise], 'FaceColor', [0.6, 0.6, 0.6]);
xlabel('Peak SNR');
ylabel('Number of units');
set(gca, 'TickDir', 'out');
box off;
meanSNR         = mean([allRes.sigToNoise]);
semSNR          = std([allRes.sigToNoise]) / (sqrt(size([allRes.sigToNoise], 2)));

% print figure
print(snrFig, fullfile(paths.save, 'peakSNR'), '-dsvg', '-r300');

% number of units per region
brainRegionIdx        = {allRes.brainRegion}';
brainRegionIdxSplit   = split(brainRegionIdx, '_');
brainRegion           = unique(brainRegionIdxSplit(:, 1));

% preallocate
numPerRegion          = [brainRegion, num2cell(NaN(size(brainRegion, 1), 1))];

% loop through brain regions
for iReg = 1:size(brainRegion, 1)
    regIdx                = strcmp(brainRegionIdxSplit(:, 1), brainRegion{iReg});
    numPerRegion{iReg, 2} = sum(regIdx);
end
