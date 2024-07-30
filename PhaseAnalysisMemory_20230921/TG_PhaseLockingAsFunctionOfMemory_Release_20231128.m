%==========================================================================
% This script applies surrogate tests for comparing phase locking between
% successful and unsuccessful encoding-recall pairs.
%
% Tim Guth, 2023
%==========================================================================

% start
clear; close all; clc;

% add paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % phase results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisMemory_20230921';

% add functions
addpath(genpath('D:\External\Functions\'));

% load result
results                 = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

% set random seed
randSeed                = 444;
rng(randSeed, 'twister');
randSeedNum             = randi(100000, 100000, 1); % for randomseed

% settings
param.nSur              = 10001;
param.polarHistLabels   = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
param.locRecThreshold   = 'median'; % 'median', or 'absolute'
param.recallType        = 'all'; % 'all', 'object', or 'location'
param.phaseLocking      = 'all'; % 'all', 'sig', or 'no'
param.unitType          = 'all'; % 'all', 'single', or 'multi'
param.cellType          = 'all'; % 'all', 'object', or 'noobject'
param.powLvl            = 'all'; % 'all', 'high', or 'low'
param.slopeLvl          = 'all'; % 'all', 'high', or 'low'
param.thetaLvl          = 'all'; % 'all', 'theta' or 'notheta'
param.resultName        = strcat(param.recallType, 'Recall', '_', ...
    param.phaseLocking, 'PL', '_', ...
    param.unitType, 'Units', '_', ...
    param.cellType, 'Cells', '_', ...
    param.powLvl, 'Powers', '_', ...
    param.slopeLvl, 'Slopes', '_', ...
    param.thetaLvl, 'Osc', '_');

% include only specified units (all, single units or multi units)
if strcmp(param.unitType, 'all') % all cells
    bUnit = true(1, size([results.phaseRes], 1));
elseif strcmp(param.unitType, 'single')
    bUnit = [results.phaseRes.bSingleUnit];
elseif strcmp(param.unitType, 'multi')
    bUnit = ~[results.phaseRes.bSingleUnit];
end

% include only specified units (all, object cells, non-object cells)
if strcmp(param.cellType, 'all') % all cells
    bCell = true(1, size([results.phaseRes], 1));
elseif strcmp(param.cellType, 'object') % cells with increased firing rate during encoding
    bCell = [results.phaseRes.bObjectCell];
elseif strcmp(param.cellType, 'noobject')
    bCell = ~[results.phaseRes.bObjectCell];
end

% include only specified units (all, significantly phase locking, non-PL cells)
if strcmp(param.phaseLocking, 'all') % all cells
    bSig = true(1, size([results.phaseRes], 1));
elseif strcmp(param.phaseLocking, 'sig')
    bSig = [results.phaseRes.encBothRank] > 0.95 & [results.phaseRes.recBothRank] > 0.95;
elseif strcmp(param.phaseLocking, 'no')
    bSig = ~([results.phaseRes.encBothRank] > 0.95 & [results.phaseRes.recBothRank] > 0.95);
end

% extract phase results
bInc                    = bUnit & bCell & bSig;
phaseRes                = results.phaseRes(bInc);

% preallocate mean phases for succ. and unsucc. encoding and recall
encZval                 = nan(size(phaseRes, 1), 1);
recZval                 = nan(size(phaseRes, 1), 1);
encSuccZval             = nan(size(phaseRes, 1), 1);
encFailZval             = nan(size(phaseRes, 1), 1);
recSuccZval             = nan(size(phaseRes, 1), 1);
recFailZval             = nan(size(phaseRes, 1), 1);
encMrv                  = nan(size(phaseRes, 1), 1);
recMrv                  = nan(size(phaseRes, 1), 1);
encSuccMrv              = nan(size(phaseRes, 1), 1);
encFailMrv              = nan(size(phaseRes, 1), 1);
recSuccMrv              = nan(size(phaseRes, 1), 1);
recFailMrv              = nan(size(phaseRes, 1), 1);
encAngle                = nan(size(phaseRes, 1), 1);
recAngle                = nan(size(phaseRes, 1), 1);
encSuccAngle            = nan(size(phaseRes, 1), 1);
encFailAngle            = nan(size(phaseRes, 1), 1);
recSuccAngle            = nan(size(phaseRes, 1), 1);
recFailAngle            = nan(size(phaseRes, 1), 1);
encCellRank             = nan(size(phaseRes, 1), 1);
recCellRank             = nan(size(phaseRes, 1), 1);
surEncSuccZval          = nan(size(phaseRes, 1), param.nSur);
surEncFailZval          = nan(size(phaseRes, 1), param.nSur);
surRecSuccZval          = nan(size(phaseRes, 1), param.nSur);
surRecFailZval          = nan(size(phaseRes, 1), param.nSur);
smpEncSuccZval          = nan(size(phaseRes, 1), param.nSur);
smpEncFailZval          = nan(size(phaseRes, 1), param.nSur);
smpRecSuccZval          = nan(size(phaseRes, 1), param.nSur);
smpRecFailZval          = nan(size(phaseRes, 1), param.nSur);
surEncSuccAngle         = nan(size(phaseRes, 1), param.nSur);
surEncFailAngle         = nan(size(phaseRes, 1), param.nSur);
surRecSuccAngle         = nan(size(phaseRes, 1), param.nSur);
surRecFailAngle         = nan(size(phaseRes, 1), param.nSur);

% loop through units and unfold the trial-wise information
parfor iCell = 1:size(phaseRes, 1)

    % initialize variables
    bGoodMemEnc     = [];
    bGoodMemRec     = [];
    encRecallIdx    = [];
    recRecallIdx    = [];
    encPowerIdx     = [];
    recPowerIdx     = [];
    encSlopeIdx     = [];
    recSlopeIdx     = [];
    encThetaIdx     = [];
    recThetaIdx     = [];
    recNotEmptyIdx  = [];
    smpEncPhaseSucc = [];
    smpEncPhaseFail = [];
    smpRecPhaseSucc = [];
    smpRecPhaseFail = [];

    % set random seed
    rng(randSeedNum(iCell, 1), 'twister');

    % display progress
    disp(iCell);

    %% spike-wise memory performance

    % number of spikes
    encNumSpikes            = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).encSpikePower, 'UniformOutput', 0));
    recNumSpikes            = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).recSpikePower, 'UniformOutput', 0));

    % memory performance index per spike
    if strcmp(param.locRecThreshold, 'median')
        bGoodMemEnc         = repelem(phaseRes(iCell).bGoodMemEnc, encNumSpikes);
        bGoodMemRec         = repelem(phaseRes(iCell).bGoodMemRec, recNumSpikes);
    elseif strcmp(param.locRecThreshold, 'absolute')
        bGoodMemEnc         = repelem(phaseRes(iCell).bGoodMemEncAbs, encNumSpikes);
        bGoodMemRec         = repelem(phaseRes(iCell).bGoodMemRecAbs, recNumSpikes);
    end

    %% extract phases, power and slope

    % phases
    allEncPhase         = cat(1, phaseRes(iCell).encPhase{:});
    allRecPhase         = cat(1, phaseRes(iCell).recPhase{:});

    % power
    encPower            = cat(1, phaseRes(iCell).encSpikePower{:});
    recPower            = cat(1, phaseRes(iCell).recSpikePower{:});

    % theta
    encTheta            = cat(1, phaseRes(iCell).encThetaIdx{:});
    recTheta            = cat(1, phaseRes(iCell).recThetaIdx{:});

    % slope
    encSlope            = cat(1, phaseRes(iCell).encSpikeSlope{:});
    recSlope            = cat(1, phaseRes(iCell).recSpikeSlope{:});

    %% inclusion index for spikes according to settings

    % extract only spikes from specific recall type
    if strcmp(param.recallType, 'all')
        encRecallIdx = ones(size(allEncPhase));
        recRecallIdx = ones(size(allRecPhase));
    elseif strcmp(param.recallType, 'object')
        encRecallIdx = repelem(phaseRes(iCell).encObjOrLoc, encNumSpikes) == 1;
        recRecallIdx = repelem(phaseRes(iCell).recObjOrLoc, recNumSpikes) == 1;
    elseif strcmp(param.recallType, 'location')
        encRecallIdx = repelem(phaseRes(iCell).encObjOrLoc, encNumSpikes) == 2;
        recRecallIdx = repelem(phaseRes(iCell).recObjOrLoc, recNumSpikes) == 2;
    end

    % extract only spikes with specific power
    if strcmp(param.powLvl, 'all')
        encPowerIdx     = ones(size(encPower));
        recPowerIdx     = ones(size(recPower));
    elseif strcmp(param.powLvl, 'high')
        encPowerIdx     = encPower > median(encPower);
        recPowerIdx     = recPower > median(recPower);
    elseif strcmp(param.powLvl, 'low')
        encPowerIdx     = encPower <= median(encPower);
        recPowerIdx     = recPower <= median(recPower);
    end

    % extract only spikes with specific slope
    if strcmp(param.slopeLvl, 'all')
        encSlopeIdx     = ones(size(encSlope));
        recSlopeIdx     = ones(size(recSlope));
    elseif strcmp(param.slopeLvl, 'high')
        encSlopeIdx     = encSlope > median(encSlope, 'omitnan');
        recSlopeIdx     = recSlope > median(recSlope, 'omitnan');
    elseif strcmp(param.slopeLvl, 'low')
        encSlopeIdx     = encSlope <= median(encSlope, 'omitnan');
        recSlopeIdx     = recSlope <= median(recSlope, 'omitnan');
    end

    % extract only spikes with/without theta oscillations
    if strcmp(param.thetaLvl, 'all')
        encThetaIdx     = ones(size(allEncPhase));
        recThetaIdx     = ones(size(allRecPhase));
    elseif strcmp(param.thetaLvl, 'theta')
        encThetaIdx     = encTheta == 1;
        recThetaIdx     = recTheta == 1;
    elseif strcmp(param.thetaLvl, 'notheta')
        encThetaIdx     = encTheta == 0;
        recThetaIdx     = recTheta == 0;
    end

    % inclusion index
    encIncIdx           = encRecallIdx & encPowerIdx & encSlopeIdx & encThetaIdx;
    recIncIdx           = recRecallIdx & recPowerIdx & recSlopeIdx & recThetaIdx;

    %% separate successful and unsuccessful segments

    % successful and unsuccessful encoding
    encPhase            = allEncPhase(encIncIdx);
    encPhaseSucc        = allEncPhase(encIncIdx & bGoodMemEnc);
    encPhaseFail        = allEncPhase(encIncIdx & ~bGoodMemEnc);

    % successful and unsuccessful recall
    recPhase            = allRecPhase(recIncIdx);
    recPhaseSucc        = allRecPhase(recIncIdx & bGoodMemRec);
    recPhaseFail        = allRecPhase(recIncIdx & ~bGoodMemRec);

    % concatenate successful and failed spikes
    encPhaseSuccFail    = cat(1, encPhaseSucc, encPhaseFail);
    recPhaseSuccFail    = cat(1, recPhaseSucc, recPhaseFail);

    % class variables
    encSuccFailClass    = cat(1 , ones(size(encPhaseSucc)), zeros(size(encPhaseFail)));
    recSuccFailClass    = cat(1 , ones(size(recPhaseSucc)), zeros(size(recPhaseFail)));

    try
        % get mean resultant vector length and mean angle
        [encMrv(iCell, 1), encAngle(iCell, 1)]         = circ_axialmean(encPhase);
        [recMrv(iCell, 1), recAngle(iCell, 1)]         = circ_axialmean(recPhase);
        [encSuccMrv(iCell, 1), encSuccAngle(iCell, 1)] = circ_axialmean(encPhaseSucc);
        [encFailMrv(iCell, 1), encFailAngle(iCell, 1)] = circ_axialmean(encPhaseFail);
        [recSuccMrv(iCell, 1), recSuccAngle(iCell, 1)] = circ_axialmean(recPhaseSucc);
        [recFailMrv(iCell, 1), recFailAngle(iCell, 1)] = circ_axialmean(recPhaseFail);

        % get Rayleigh z-value
        [~, encZval(iCell, 1)]       = circ_rtest(encPhase);
        [~, recZval(iCell, 1)]       = circ_rtest(recPhase);
        [~, encSuccZval(iCell, 1)]   = circ_rtest(encPhaseSucc);
        [~, encFailZval(iCell, 1)]   = circ_rtest(encPhaseFail);
        [~, recSuccZval(iCell, 1)]   = circ_rtest(recPhaseSucc);
        [~, recFailZval(iCell, 1)]   = circ_rtest(recPhaseFail);

        %% surrogate data

        % preallocate
        cellSurEncSuccZval          = nan(1, param.nSur);
        cellSurEncFailZval          = nan(1, param.nSur);
        cellSurRecSuccZval          = nan(1, param.nSur);
        cellSurRecFailZval          = nan(1, param.nSur);
        cellSmpEncSuccZval          = nan(1, param.nSur);
        cellSmpEncFailZval          = nan(1, param.nSur);
        cellSmpRecSuccZval          = nan(1, param.nSur);
        cellSmpRecFailZval          = nan(1, param.nSur);
        cellSurEncSuccAngle         = nan(1, param.nSur);
        cellSurEncFailAngle         = nan(1, param.nSur);
        cellSurRecSuccAngle         = nan(1, param.nSur);
        cellSurRecFailAngle         = nan(1, param.nSur);

        % loop through surrogates
        for iSur = 1:param.nSur

            % random shuffle of spike class assignment
            surEncSuccFailClass = datasample(encSuccFailClass, size(encSuccFailClass, 1), 'Replace', false);
            surRecSuccFailClass = datasample(recSuccFailClass, size(recSuccFailClass, 1), 'Replace', false);

            % create surrogates
            surEncPhaseSucc     = encPhaseSuccFail(surEncSuccFailClass == 1);
            surEncPhaseFail     = encPhaseSuccFail(surEncSuccFailClass == 0);
            surRecPhaseSucc     = recPhaseSuccFail(surRecSuccFailClass == 1);
            surRecPhaseFail     = recPhaseSuccFail(surRecSuccFailClass == 0);

            % get mean resultant vector length and mean angle
            [~, cellSurEncSuccAngle(1, iSur)] = circ_axialmean(surEncPhaseSucc);
            [~, cellSurEncFailAngle(1, iSur)] = circ_axialmean(surEncPhaseFail);
            [~, cellSurRecSuccAngle(1, iSur)] = circ_axialmean(surRecPhaseSucc);
            [~, cellSurRecFailAngle(1, iSur)] = circ_axialmean(surRecPhaseFail);

            % get Rayleigh z-value
            [~, cellSurEncSuccZval(1, iSur)]  = circ_rtest(surEncPhaseSucc);
            [~, cellSurEncFailZval(1, iSur)]  = circ_rtest(surEncPhaseFail);
            [~, cellSurRecSuccZval(1, iSur)]  = circ_rtest(surRecPhaseSucc);
            [~, cellSurRecFailZval(1, iSur)]  = circ_rtest(surRecPhaseFail);

            %% subsampling to correct for spike numbers - encoding

            % compare size of successful and unsuccessful encoding phases
            if size(encPhaseSucc, 1) > size(encPhaseFail, 1)

                % draw a subsample from the successful phases
                smpEncPhaseSucc = datasample(encPhaseSucc, size(encPhaseFail, 1), 'Replace', false);
                smpEncPhaseFail = encPhaseFail;

            elseif size(encPhaseSucc, 1) < size(encPhaseFail, 1)

                % draw a subsample from the unsuccessful phases
                smpEncPhaseSucc = encPhaseSucc;
                smpEncPhaseFail = datasample(encPhaseFail, size(encPhaseSucc, 1), 'Replace', false);

            elseif size(encPhaseSucc, 1) == size(encPhaseFail, 1)

                % no subsampling
                smpEncPhaseSucc = encPhaseSucc;
                smpEncPhaseFail = encPhaseFail;
            end

            % Rayleigh z-value
            [~, cellSmpEncSuccZval(1, iSur)] = circ_rtest(smpEncPhaseSucc);
            [~, cellSmpEncFailZval(1, iSur)] = circ_rtest(smpEncPhaseFail);

            %% subsampling to correct for spike numbers - recall

            % compare size of successful and unsuccessful recall phases
            if size(recPhaseSucc, 1) > size(recPhaseFail, 1)

                % draw a subsample from the successful phases
                smpRecPhaseSucc = datasample(recPhaseSucc, size(recPhaseFail, 1), 'Replace', false);
                smpRecPhaseFail = recPhaseFail;

            elseif size(recPhaseSucc, 1) < size(recPhaseFail, 1)

                % draw a subsample from the unsuccessful phases
                smpRecPhaseSucc = recPhaseSucc;
                smpRecPhaseFail = datasample(recPhaseFail, size(recPhaseSucc, 1), 'Replace', false);

            elseif size(recPhaseSucc, 1) == size(recPhaseFail, 1)

                % no subsampling
                smpRecPhaseSucc = recPhaseSucc;
                smpRecPhaseFail = recPhaseFail;
            end

            % Rayleigh z-value
            [~, cellSmpRecSuccZval(1, iSur)] = circ_rtest(smpRecPhaseSucc);
            [~, cellSmpRecFailZval(1, iSur)] = circ_rtest(smpRecPhaseFail);
        end
        
        %% test for significant shift for this cell

        % difference in z-values
        encSuccFailDiff             = encSuccZval(iCell, 1) - encFailZval(iCell, 1);
        recSuccFailDiff             = recSuccZval(iCell, 1) - recFailZval(iCell, 1);

        % surrogate differences in z-values
        surEncSuccFailDiff          = cellSurEncSuccZval - cellSurEncFailZval;
        surRecSuccFailDiff          = cellSurRecSuccZval - cellSurRecFailZval;

        % ranks of empirical z-value difference in surrogate differences
        encCellRank(iCell, :)       = sum(encSuccFailDiff > surEncSuccFailDiff) / size(surEncSuccFailDiff, 2);
        recCellRank(iCell, :)       = sum(recSuccFailDiff > surRecSuccFailDiff) / size(surRecSuccFailDiff, 2);

        %% collect surrogates across cells
        surEncSuccZval(iCell, :)    = cellSurEncSuccZval;
        surEncFailZval(iCell, :)    = cellSurEncFailZval;
        surRecSuccZval(iCell, :)    = cellSurRecSuccZval;
        surRecFailZval(iCell, :)    = cellSurRecFailZval;
        smpEncSuccZval(iCell, :)    = cellSmpEncSuccZval;
        smpEncFailZval(iCell, :)    = cellSmpEncFailZval;
        smpRecSuccZval(iCell, :)    = cellSmpRecSuccZval;
        smpRecFailZval(iCell, :)    = cellSmpRecFailZval;
        surEncSuccAngle(iCell, :)   = cellSurEncSuccAngle;
        surEncFailAngle(iCell, :)   = cellSurEncFailAngle;
        surRecSuccAngle(iCell, :)   = cellSurRecSuccAngle;
        surRecFailAngle(iCell, :)   = cellSurRecFailAngle;
    catch
        continue;
    end
end

%% t-tests for phase locking during encoding and recall across units

% empirical t-tests
[~, ~, ~, encStat]      = ttest(encSuccZval, encFailZval);
encT                    = encStat.tstat;
[~, ~, ~, recStat]      = ttest(recSuccZval, recFailZval);
recT                    = recStat.tstat;

% surrogate t-tests
[~, ~, ~, surEncStat]   = ttest(surEncSuccZval, surEncFailZval);
surEncT                 = surEncStat.tstat;
[~, ~, ~, surRecStat]   = ttest(surRecSuccZval, surRecFailZval);
surRecT                 = surRecStat.tstat;

% subsampling t-tests
[~, smpEncPval, ~, smpEncStat] = ttest(mean(smpEncSuccZval, 2, 'omitnan'), mean(smpEncFailZval, 2, 'omitnan'));
[~, smpRecPval, ~, smpRecStat] = ttest(mean(smpRecSuccZval, 2, 'omitnan'), mean(smpRecFailZval, 2, 'omitnan'));

% Bonferroni correction (successful and unsuccessful)
smpEncPval              = smpEncPval * 2;
smpRecPval              = smpRecPval * 2;

%% Watson-William test

% inclusion index
encIncIdx               = ~isnan(encSuccAngle) & ~isnan(encFailAngle);
recIncIdx               = ~isnan(recSuccAngle) & ~isnan(recFailAngle);

% mean difference
meanDiffEnc             = circ_mean(angdiff(encSuccAngle(encIncIdx), encFailAngle(encIncIdx)));
meanDiffRec             = circ_mean(angdiff(recSuccAngle(recIncIdx), recFailAngle(recIncIdx)));

% empirical Watson-William test
[~, wwEncTable]         = circ_wwtest(encSuccAngle(encIncIdx), encFailAngle(encIncIdx));
[~, wwRecTable]         = circ_wwtest(recSuccAngle(recIncIdx), recFailAngle(recIncIdx));
encWw                   = wwEncTable{2, 5};
recWw                   = wwRecTable{2, 5};

% surrogate Watson-William test
surEncWw = nan(param.nSur, 1);
surRecWw = nan(param.nSur, 1);
for iSur = 1:param.nSur

    % inclusion index
    surIncIdxEnc        = ~isnan(surEncSuccAngle(:, iSur)) & ~isnan(surEncFailAngle(:, iSur));
    surIncIdxRec        = ~isnan(surRecSuccAngle(:, iSur)) & ~isnan(surRecFailAngle(:, iSur));

    % surrogate Watson-William test
    [~, surEncWwTable]  = circ_wwtest(surEncSuccAngle(surIncIdxEnc, iSur), surEncFailAngle(surIncIdxEnc, iSur));
    [~, surRecWwTable]  = circ_wwtest(surRecSuccAngle(surIncIdxRec, iSur), surRecFailAngle(surIncIdxRec, iSur));
    surEncWw(iSur, 1)   = surEncWwTable{2, 5};
    surRecWw(iSur, 1)   = surRecWwTable{2, 5};
end

% surrogate test results - t-tests
encRank         = sum(encT > surEncT) / param.nSur;
encPval         = (1 - encRank) * 2; % Bonferroni correction (successful and unsuccessful)
recRank         = sum(recT > surRecT) / param.nSur;
recPval         = (1 - recRank) * 2;

% surrogate test results - Watson-Williams tests
encWwRank       = sum(encWw > surEncWw) / param.nSur;
encWwPval       = (1 - encWwRank) * 2;
recWwRank       = sum(recWw > surRecWw) / param.nSur;
recWwPval       = (1 - recWwRank) * 2;

% additional Bonferroni correction, if two groups are compared
if ~strcmp(param.recallType, 'all') || ...
        ~strcmp(param.phaseLocking, 'all') || ...
        ~strcmp(param.unitType, 'all') || ...
        ~strcmp(param.cellType, 'all') || ...
        ~strcmp(param.powLvl, 'all') || ...
        ~strcmp(param.slopeLvl, 'all') || ...
        ~strcmp(param.thetaLvl, 'all')
    encPval     = encPval * 2;
    recPval     = recPval * 2;
    encWwPval   = encWwPval * 2;
    recWwPval   = recWwPval * 2;
end

%% plots for all spikes
if strcmp(param.recallType, 'all') && ...
        strcmp(param.phaseLocking, 'all') && ...
        strcmp(param.unitType, 'all') && ...
        strcmp(param.cellType, 'all') && ...
        strcmp(param.powLvl, 'all') && ...
        strcmp(param.slopeLvl, 'all') && ...
        strcmp(param.thetaLvl, 'all')

    %% plot mean angle and zval distributions

    % plot all phases for successful encoding units
    encSuccPolar = figure;
    polarplot([encSuccAngle, encSuccAngle]', [zeros(size(encSuccAngle, 1), 1), encSuccMrv]', 'Color', [0.1, 0.6, 0.1, 0.6], 'LineWidth', 1);
    hold on;
    polarplot([circ_mean(encSuccAngle), circ_mean(encSuccAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    rticks([0, 0.5, 1]);
    title('Preffered theta phase during successful encoding');
    print(encSuccPolar, fullfile(paths.save, 'allCells_encSucc_PolarHistogram'), '-dsvg', '-r300');

    % plot all phases for unsuccessful encoding units
    encFailPolar = figure;
    polarplot([encFailAngle, encFailAngle]', [zeros(size(encFailAngle, 1), 1), encFailMrv]', 'Color', [0.6, 0.1, 0.1, 0.6], 'LineWidth', 1);
    hold on;
    polarplot([circ_mean(encFailAngle), circ_mean(encFailAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    rticks([0, 0.5, 1]);
    title('Preffered theta phase during unsuccessful encoding');
    print(encFailPolar, fullfile(paths.save, 'allCells_encFail_PolarHistogram'), '-dsvg', '-r300');

    % plot all phases for successful recall units
    recSuccPolar = figure;
    polarplot([recSuccAngle, recSuccAngle]', [zeros(size(recSuccAngle, 1), 1), recSuccMrv]', 'Color', [0.1, 0.6, 0.1, 0.6], 'LineWidth', 1);
    hold on;
    polarplot([circ_mean(recSuccAngle), circ_mean(recSuccAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    rticks([0, 0.5, 1]);
    title('Preffered theta phase during successful recall');
    print(recSuccPolar, fullfile(paths.save, 'allCells_recSucc_PolarHistogram'), '-dsvg', '-r300');

    % plot all phases for unsuccessful recall units
    recFailPolar = figure;
    polarplot([recFailAngle, recFailAngle]', [zeros(size(recFailAngle, 1), 1), recFailMrv]', 'Color', [0.6, 0.1, 0.1, 0.6], 'LineWidth', 1);
    hold on;
    polarplot([circ_mean(recFailAngle), circ_mean(recFailAngle)], [0, 1], 'Color', [0, 0, 0], 'LineWidth', 2);
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    rticks([0, 0.5, 1]);
    title('Preffered theta phase during unsuccessful recall');
    print(recFailPolar, fullfile(paths.save, 'allCells_recFail_PolarHistogram'), '-dsvg', '-r300');

    %% plot histogram of surrogate and empirical t-value - encoding
    encSurFigure       = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    encSurHist         = histogram(surEncT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(encT, 'color', 'k', 'LineWidth', 1);
    % xlim([-10, 10]);
    % ylim([0, 500]);
    set(gca,'TickDir','out');
    box off;
    saveas(encSurFigure, fullfile(paths.save, strcat(param.resultName, '_EncodingSurDistr.svg')));

    %% plot histogram of surrogate and empirical t-value - recall
    recSurFigure       = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    recSurHist         = histogram(surRecT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(recT, 'color', 'k', 'LineWidth', 1);
    % xlim([-10, 10]);
    % ylim([0, 500]);
    set(gca,'TickDir','out');
    box off;
    saveas(recSurFigure, fullfile(paths.save, strcat(param.resultName, '_RecallSurDistr.svg')));
end

%% plot mean zval of successful and unsuccessful encoding

% figure
encZvalFig      = figure('units', 'centimeters', 'position', [2, 2, 10, 15]);

% mean
encSuccZvalMean = mean(mean(smpEncSuccZval, 2, 'omitnan'), 'omitnan');
encFailZvalMean = mean(mean(smpEncFailZval, 2, 'omitnan'), 'omitnan');

% standard error of the mean
semEncSucc      = std(mean(smpEncSuccZval, 2, 'omitnan'), 'omitnan') / sqrt(size(encSuccZval(encIncIdx), 1));
semEncFail      = std(mean(smpEncFailZval, 2, 'omitnan'), 'omitnan') / sqrt(size(encFailZval(encIncIdx), 1));

% bar plot
b1              = bar([encSuccZvalMean; encFailZvalMean], 'BarWidth', 0.6);
b1.FaceColor    = 'flat';
b1.CData        = [0.1, 0.6, 0.1; 0.6, 0.1, 0.1];
hold on;
e1              = errorbar([1; 2], [encSuccZvalMean; encFailZvalMean], [semEncSucc; semEncFail]);
e1.Color        = [0 0 0];
e1.LineStyle    = 'none';
if encPval < 0.05
    hold on;
    l1          = line([1, 2], [0.4, 0.4], 'LineWidth', 2, 'Color', 'k');
    a1          = plot(1.5, 0.41, '*', 'Color', 'k');
end
xlim([0.3, 2.7]);
set(gca, 'xticklabel', {'Successful', 'Unsuccessful'});
ylabel('Mean zval');
title({'zval (mean across units)',  'during encoding', ...
    strcat('(p = ', num2str(encPval), ')')});

% set y-axis
ax          = gca;
yMax        = ceil(max(ax.YTick));
ylim([0, yMax]);
if ~(strcmp(param.recallType, 'all') && ...
        strcmp(param.phaseLocking, 'all') && ...
        strcmp(param.unitType, 'all') && ...
        strcmp(param.cellType, 'all') && ...
        strcmp(param.powLvl, 'all') && ...
        strcmp(param.slopeLvl, 'all') && ...
        strcmp(param.thetaLvl, 'all'))
    ax.YTick    = [0, yMax];
end
ax.TickDir  = 'out';
box off;

% print figure
print(encZvalFig, fullfile(paths.save, strcat(param.resultName, '_EncodingZval')), '-dsvg', '-r300');

%% plot mean zval of successful and unsuccessful recall

% figure
recZvalFig      = figure('units', 'centimeters', 'position', [2, 2, 10, 15]);

% mean
recSuccZvalMean = mean(mean(smpRecSuccZval, 2, 'omitnan'), 'omitnan');
recFailZvalMean = mean(mean(smpRecFailZval, 2, 'omitnan'), 'omitnan');

% standard error of the mean
semRecSucc      = std(mean(smpRecSuccZval, 2, 'omitnan'), 'omitnan') / sqrt(size(recSuccZval(recIncIdx), 1));
semRecFail      = std(mean(smpRecFailZval, 2, 'omitnan'), 'omitnan') / sqrt(size(recFailZval(recIncIdx), 1));

% bar plot
b1              = bar([recSuccZvalMean; recFailZvalMean], 'BarWidth', 0.6);
b1.FaceColor    = 'flat';
b1.CData        = [0.1, 0.6, 0.1; 0.6, 0.1, 0.1];
hold on;
e1              = errorbar([1; 2], [recSuccZvalMean; recFailZvalMean], [semRecSucc; semRecFail]);
e1.Color        = [0 0 0];
e1.LineStyle    = 'none';
if recPval < 0.05
    hold on;
    l1          = line([1, 2], [0.4, 0.4], 'LineWidth', 2, 'Color', 'k');
    a1          = plot(1.5, 0.41, '*', 'Color', 'k');
end
xlim([0.3, 2.7]);
set(gca, 'xticklabel', {'Successful', 'Unsuccessful'});
ylabel('Mean zval');
title({'zval (mean across units)',  'during recall', ...
    strcat('(p = ', num2str(recPval), ')')});

% set y-axis
ax          = gca;
yMax        = ceil(max(ax.YTick));
ylim([0, yMax]);
if ~(strcmp(param.recallType, 'all') && ...
        strcmp(param.phaseLocking, 'all') && ...
        strcmp(param.unitType, 'all') && ...
        strcmp(param.cellType, 'all') && ...
        strcmp(param.powLvl, 'all') && ...
        strcmp(param.slopeLvl, 'all') && ...
        strcmp(param.thetaLvl, 'all'))
    ax.YTick = [0, yMax];
end
ax.TickDir  = 'out';
box off;

% print figure
print(recZvalFig, fullfile(paths.save, strcat(param.resultName, '_RecallZval')), '-dsvg', '-r300');