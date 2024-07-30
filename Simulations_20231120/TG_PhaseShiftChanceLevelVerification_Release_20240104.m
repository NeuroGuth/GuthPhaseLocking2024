%==========================================================================
% This script empirically checks whether the a priori chosen chance level
% of 5% for phase shifting units in the script "TG_PhaseShiftAnalysis" is
% justified.
%
% Tim Guth, 2024
%==========================================================================

% start
clear; close all; clc;

% add paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % results folder
paths.save              = 'D:\TreasureHunt\Simulations_20231120';

% add functions
addpath(genpath('D:\External\Functions\'));

%% settings

% set random seed
randSeedNum             = 444;
rng(randSeedNum, 'twister');
randSeedNum             = randi(100000, 100000, 1); % for randomseed

% parameters
param                   = [];
param.nSur              = 1001;

% load result
allResults              = load(fullfile(paths.phaseRes, 'allResults.mat'));
additionalResults       = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

% extract phase results
allRes                  = allResults.allRes(~additionalResults.bExclude);
phaseRes                = additionalResults.phaseRes;

% extract structure fields of interest
allEncRecIdx            = cat(2, {allRes.allEncRecIdx});
allBGoodMemEnc          = cat(2, {phaseRes.bGoodMemEnc});
allBGoodMemRec          = cat(2, {phaseRes.bGoodMemRec});

% preallocate
encRecRank              = nan(size(phaseRes, 1), param.nSur);

%% loop through units and unfold the trial-wise information
parfor iCell = 1:size(phaseRes, 1)

    % display progress
    disp(iCell);

    % set variables
    cellEncRecRank      = nan(1, param.nSur);
    cellEncRecRankSucc  = nan(1, param.nSur);
    cellEncRecRankFail  = nan(1, param.nSur);

    % set random seed
    rng(randSeedNum(iCell, 1), 'twister');

    % run simulation n times
    for iSim = 1:param.nSur

        %% simulate random phase-locking neuron

        % empirical phases
        encPhasesOrig                       = phaseRes(iCell).encPhase;
        recPhasesOrig                       = phaseRes(iCell).recPhase;
        
        % concatenate phases
        encCatPhasesOrig                    = cat(1, encPhasesOrig{:});
        recCatPhasesOrig                    = cat(1, recPhasesOrig{:});

        % encoding-retrieval index
        encRecIdx                           = allEncRecIdx{iCell}(allEncRecIdx{iCell} ~= 0);

        % concatenate encoding and retrieval spikes
        encRecCatPhasesOrig                 = nan(size(encRecIdx));
        encRecCatPhasesOrig(encRecIdx == 1) = encCatPhasesOrig;
        encRecCatPhasesOrig(encRecIdx == 2) = recCatPhasesOrig;
        
        % circularly shift concatenated encoding and retrieval phases
        circEncRecCatPhasesOrig             = circshift(encRecCatPhasesOrig, randi(size(encRecCatPhasesOrig, 1)));

        % extract surrogate encoding and recall phases
        randEncPhases                       = circEncRecCatPhasesOrig(encRecIdx == 1);
        randRecPhases                       = circEncRecCatPhasesOrig(encRecIdx == 2);

        %% number of spikes

        % number of spikes per trial
        encNumSpikes    = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).encSpikePower, 'UniformOutput', 0));
        recNumSpikes    = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).recSpikePower, 'UniformOutput', 0));

        %% get encoding phases

        % recreate segment-wise results
        encPhases       = mat2cell(randEncPhases, encNumSpikes, 1);

        %% get recall phases

        % recreate segment-wise results
        recPhases       = mat2cell(randRecPhases, recNumSpikes, 1);

        % sort the recall phases so that they match the encoding segments
        recPhasesSorted = cell(size(encPhases, 1), 1);
        for iS = 1:size(phaseRes(iCell).encSegmentIdx, 1)
            logIdx                    = phaseRes(iCell).recSegmentIdx == phaseRes(iCell).encSegmentIdx(iS);
            recPhasesSorted(iS, 1)    = recPhases(logIdx);
        end

        %% Watson-Williams test for phase shifts

        % concatenate phases of two groups
        encRecPhases                = cat(1, encPhases, recPhasesSorted);

        % class variable
        classEncRec                 = [ones(size(encPhases, 1), 1); ones(size(recPhasesSorted, 1), 1) * 2];

        % empirical Watson-Williams test
        [~, wwEncRecTable]          = circ_wwtest(cat(1, encRecPhases{classEncRec == 1}), cat(1, encRecPhases{classEncRec == 2}));

        % F-values
        wwEncRec                    = wwEncRecTable{2, 5};

        % surrogate tests
        cellSurWwEncRec             = nan(param.nSur, 1);
        for iSur = 1:param.nSur

            % shuffle class for surrogate dataset
            surClassEncRec                  = datasample(classEncRec, numel(classEncRec), 'replace', false);

            % surrogate Watson-Williams test
            [~, surWwEncRecTable]           = circ_wwtest(cat(1, encRecPhases{surClassEncRec == 1}), cat(1, encRecPhases{surClassEncRec == 2}));

            % surrogate F-values
            cellSurWwEncRec(iSur, 1)        = surWwEncRecTable{2, 5};
        end

        % rank of empirical WW-test F-value in surrogate dataset
        cellEncRecRank(1, iSim)         = sum(wwEncRec > cellSurWwEncRec) / param.nSur;
    end

    % collect results across cells
    encRecRank(iCell, :)        = cellEncRecRank;
end

%% find cells with significant phase shift

% phase locking and phase difference
bShift      = encRecRank > 0.95;

% binomial tests
allPval     = myBinomTest(sum(bShift), size(bShift, 1), 0.05, 'one');

%% plot distribution of chance level
chanceLvl = figure;
histogram(sum(bShift) / size(bShift, 1), 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4], 'Orientation', 'horizontal');
yline(mean(sum(bShift) / size(bShift, 1)));
ylim([0, 0.15]);
ylabel('Chance level');
set(gca, 'tickDir', 'out');
box off;

% print figure
print(chanceLvl, fullfile(paths.save, 'ChanceLevelDistribution'), '-dsvg', '-r300');
