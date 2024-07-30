%=======================================================================
% This script tests for different phase-locking angles between
% encoding and retrieval.
%
% Tim Guth, 2023
%=======================================================================

% start
clear; close all; clc;

% add paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisShifts_20231011';

% add functions
addpath(genpath('D:\External\Functions\'));

% load result
additionalResults       = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));

%% settings

% set random seed
randSeedNum             = 444;
rng(randSeedNum, 'twister');
randSeedNum             = randi(100000, 100000, 1); % for randomseed

% parameters
param                   = [];
param.nSur              = 10001;
param.polarHistEdges    = linspace(-pi, pi, 24);
param.polarHistLabels   = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
param.unitType          = 'all'; % 'all', 'single', or 'multi'
param.cellType          = 'all'; % 'all', 'object', or 'noobject'
param.powLvl            = 'all'; % 'all', 'high', or 'low'
param.slopeLvl          = 'all'; % 'all', 'high', or 'low'
param.thetaLvl          = 'all'; % 'all', 'theta' or 'notheta'
param.resultName        = strcat(param.unitType, 'Units', '_', ...
    param.cellType, 'Cells', '_', ...
    param.powLvl, 'Powers', '_', ...
    param.slopeLvl, 'Slopes', '_', ...
    param.thetaLvl, 'Osc', '_');

% include only specified units (all, single units or multi units)
if strcmp(param.unitType, 'all') % all cells
    bUnit = true(1, size([additionalResults.phaseRes], 1));
elseif strcmp(param.unitType, 'single')
    bUnit = [additionalResults.phaseRes.bSingleUnit];
elseif strcmp(param.unitType, 'multi')
    bUnit = ~[additionalResults.phaseRes.bSingleUnit];
end

% include only specified cells
if strcmp(param.cellType, 'all') % all cells
    bSel = true(1, size([additionalResults.phaseRes], 1));
elseif strcmp(param.cellType, 'object') % cells with high firing rate during encoding
    bSel = [additionalResults.phaseRes.bObjectCell];
elseif strcmp(param.cellType, 'noobject')
    bSel = ~[additionalResults.phaseRes.bObjectCell];
end

% figure names
figureDir               = dir(fullfile(paths.phaseRes, 'UnitFigures'));
allFigureNames          = {figureDir.name}';

% extract phase results
bInc                    = bUnit & bSel;
phaseRes                = additionalResults.phaseRes(bInc);

% extract structure fields of interest
allIdx                  = cat(1, phaseRes.idx);

% preallocate mean phases for succ. and unsucc. encoding and recall
exclusionIndex          = false(size(phaseRes, 1), 1);
encAngle                = nan(size(phaseRes, 1), 1);
encSuccAngle            = nan(size(phaseRes, 1), 1);
encFailAngle            = nan(size(phaseRes, 1), 1);
recAngle                = nan(size(phaseRes, 1), 1);
recSuccAngle            = nan(size(phaseRes, 1), 1);
recFailAngle            = nan(size(phaseRes, 1), 1);
wwEncRec                = nan(size(phaseRes, 1), 1);
wwEncRecSucc            = nan(size(phaseRes, 1), 1);
wwEncRecFail            = nan(size(phaseRes, 1), 1);
surWwEncRec             = nan(size(phaseRes, 1), param.nSur);
surWwEncRecSucc         = nan(size(phaseRes, 1), param.nSur);
surWwEncRecFail         = nan(size(phaseRes, 1), param.nSur);
encRecRank              = nan(size(phaseRes, 1), 1);
encRecRankSucc          = nan(size(phaseRes, 1), 1);
encRecRankFail          = nan(size(phaseRes, 1), 1);

%% loop through units and unfold the trial-wise information
parfor iCell = 1:size(phaseRes, 1)

    % some cells may not have any spikes for a certain condition
    try

        % initialize variables
        encRecallIdx    = [];
        recRecallIdx    = [];
        encCatPowerIdx  = [];
        recCatPowerIdx  = [];
        encCatSlopeIdx  = [];
        recCatSlopeIdx  = [];
        encCatThetaIdx  = [];
        recCatThetaIdx  = [];
        smpEncPhaseSucc = [];
        smpEncPhaseFail = [];
        smpRecPhaseSucc = [];
        smpRecPhaseFail = [];

        % set random seed
        rng(randSeedNum(iCell, 1), 'twister');

        % display progress
        disp(iCell);

        %% number of spikes

        % number of spikes per trial
        encNumSpikes        = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).encSpikePower, 'UniformOutput', 0));
        recNumSpikes        = cell2mat(cellfun(@(x) size(x, 1), phaseRes(iCell).recSpikePower, 'UniformOutput', 0));

        %% extract phases, power and slope

        % phases
        allEncPhase         = phaseRes(iCell).encPhase;
        allRecPhase         = phaseRes(iCell).recPhase;

        % concatenate phases
        encCatPhase         = cat(1, allEncPhase{:});
        recCatPhase         = cat(1, allRecPhase{:});

        % concatenate power
        encCatPower         = cat(1, phaseRes(iCell).encSpikePower{:});
        recCatPower         = cat(1, phaseRes(iCell).recSpikePower{:});

        % theta
        encCatTheta         = cat(1, phaseRes(iCell).encThetaIdx{:});
        recCatTheta         = cat(1, phaseRes(iCell).recThetaIdx{:});

        % slope
        encCatSlope         = cat(1, phaseRes(iCell).encSpikeSlope{:});
        recCatSlope         = cat(1, phaseRes(iCell).recSpikeSlope{:});

        %% inclusion index for spikes according to settings

        % extract only spikes with specific power
        if strcmp(param.powLvl, 'all')
            encCatPowerIdx     = ones(size(encCatPower));
            recCatPowerIdx     = ones(size(recCatPower));
        elseif strcmp(param.powLvl, 'high')
            encCatPowerIdx     = encCatPower > median(encCatPower);
            recCatPowerIdx     = recCatPower > median(recCatPower);
        elseif strcmp(param.powLvl, 'low')
            encCatPowerIdx     = encCatPower <= median(encCatPower);
            recCatPowerIdx     = recCatPower <= median(recCatPower);
        end

        % extract only spikes with specific slope
        if strcmp(param.slopeLvl, 'all')
            encCatSlopeIdx     = ones(size(encCatSlope));
            recCatSlopeIdx     = ones(size(recCatSlope));
        elseif strcmp(param.slopeLvl, 'high')
            encCatSlopeIdx     = encCatSlope > median(encCatSlope, 'omitnan');
            recCatSlopeIdx     = recCatSlope > median(recCatSlope, 'omitnan');
        elseif strcmp(param.slopeLvl, 'low')
            encCatSlopeIdx     = encCatSlope <= median(encCatSlope, 'omitnan');
            recCatSlopeIdx     = recCatSlope <= median(recCatSlope, 'omitnan');
        end

        % extract only spikes with/without theta oscillations
        if strcmp(param.thetaLvl, 'all')
            encCatThetaIdx     = ones(size(encCatPhase));
            recCatThetaIdx     = ones(size(recCatPhase));
        elseif strcmp(param.thetaLvl, 'theta')
            encCatThetaIdx     = encCatTheta == 1;
            recCatThetaIdx     = recCatTheta == 1;
        elseif strcmp(param.thetaLvl, 'notheta')
            encCatThetaIdx     = encCatTheta == 0;
            recCatThetaIdx     = recCatTheta == 0;
        end
        
        %% include only selected results

        % concatenated inclusion index
        encCatIncIdx        = encCatPowerIdx & encCatSlopeIdx & encCatThetaIdx;
        recCatIncIdx        = recCatPowerIdx & recCatSlopeIdx & recCatThetaIdx;

        % recreate segment-wise results
        encIncIdx           = mat2cell(encCatIncIdx, encNumSpikes, 1);
        recIncIdx           = mat2cell(recCatIncIdx, recNumSpikes, 1);
        
        % include only specific phases
        encPhase            = cellfun(@(x, y) x(y == 1), allEncPhase, encIncIdx, 'Uni', 0);
        recPhase            = cellfun(@(x, y) x(y == 1), allRecPhase, recIncIdx, 'Uni', 0);

        %% separate successful and unsuccessful segments

        % memory performance index per spike
        bGoodMemEnc         = phaseRes(iCell).bGoodMemEnc;
        bGoodMemRec         = phaseRes(iCell).bGoodMemRec;

        % successful and unsuccessful encoding
        encPhaseSucc        = encPhase(bGoodMemEnc);
        encPhaseFail        = encPhase(~bGoodMemEnc);

        % successful and unsuccessful recall
        recPhaseSucc        = recPhase(bGoodMemRec);
        recPhaseFail        = recPhase(~bGoodMemRec);

        % get mean resultant vector length
        [encMrv, encAngle(iCell, 1)]            = circ_axialmean(cat(1, encPhase{:}));
        [encMrvSucc, encSuccAngle(iCell, 1)]    = circ_axialmean(cat(1, encPhaseSucc{:}));
        [encMrvFail, encFailAngle(iCell, 1)]    = circ_axialmean(cat(1, encPhaseFail{:}));
        [recMrv, recAngle(iCell, 1)]            = circ_axialmean(cat(1, recPhase{:}));
        [recMrvSucc, recSuccAngle(iCell, 1)]    = circ_axialmean(cat(1, recPhaseSucc{:}));
        [recMrvFail, recFailAngle(iCell, 1)]    = circ_axialmean(cat(1, recPhaseFail{:}));

        %% Watson-Williams test for phase shifts

        % concatenate phases of two groups
        encRecPhases                = cat(1, encPhase, recPhase);
        encRecPhasesSucc            = cat(1, encPhaseSucc, recPhaseSucc);
        encRecPhasesFail            = cat(1, encPhaseFail, recPhaseFail);

        % class variable
        classEncRec                 = [ones(size(encPhase, 1), 1); ones(size(recPhase, 1), 1) * 2];
        classEncRecSucc             = [ones(size(encPhaseSucc, 1), 1); ones(size(recPhaseSucc, 1), 1) * 2];
        classEncRecFail             = [ones(size(encPhaseFail, 1), 1); ones(size(recPhaseFail, 1), 1) * 2];

        % empirical Watson-Williams test
        [~, wwEncRecTable]          = circ_wwtest(cat(1, encRecPhases{classEncRec == 1}), cat(1, encRecPhases{classEncRec == 2}));
        [~, wwEncRecTableSucc]      = circ_wwtest(cat(1, encRecPhasesSucc{classEncRecSucc == 1}), cat(1, encRecPhasesSucc{classEncRecSucc == 2}));
        [~, wwEncRecTableFail]      = circ_wwtest(cat(1, encRecPhasesFail{classEncRecFail == 1}), cat(1, encRecPhasesFail{classEncRecFail == 2}));

        % F-values
        wwEncRec(iCell, 1)          = wwEncRecTable{2, 5};
        wwEncRecSucc(iCell, 1)      = wwEncRecTableSucc{2, 5};
        wwEncRecFail(iCell, 1)      = wwEncRecTableFail{2, 5};

        % surrogate tests
        cellSurWwEncRec             = nan(param.nSur, 1);
        cellSurWwEncRecSucc         = nan(param.nSur, 1);
        cellSurWwEncRecFail         = nan(param.nSur, 1);
        for iSur = 1:param.nSur

            % shuffle class for surrogate dataset
            surClassEncRec                  = datasample(classEncRec, numel(classEncRec), 'replace', false);
            surClassEncRecSucc              = datasample(classEncRecSucc, numel(classEncRecSucc), 'replace', false);
            surClassEncRecFail              = datasample(classEncRecFail, numel(classEncRecFail), 'replace', false);

            % surrogate Watson-Williams test
            [~, surWwEncRecTable]           = circ_wwtest(cat(1, encRecPhases{surClassEncRec == 1}), cat(1, encRecPhases{surClassEncRec == 2}));
            [~, surWwEncRecTableSucc]       = circ_wwtest(cat(1, encRecPhasesSucc{surClassEncRecSucc == 1}), cat(1, encRecPhasesSucc{surClassEncRecSucc == 2}));
            [~, surWwEncRecTableFail]       = circ_wwtest(cat(1, encRecPhasesFail{surClassEncRecFail == 1}), cat(1, encRecPhasesFail{surClassEncRecFail == 2}));

            % surrogate F-values
            cellSurWwEncRec(iSur, 1)        = surWwEncRecTable{2, 5};
            cellSurWwEncRecSucc(iSur, 1)    = surWwEncRecTableSucc{2, 5};
            cellSurWwEncRecFail(iSur, 1)    = surWwEncRecTableFail{2, 5};
        end

        % collect surrogates across cells
        surWwEncRec(iCell, :)       = cellSurWwEncRec;
        surWwEncRecSucc(iCell, :)   = cellSurWwEncRecSucc;
        surWwEncRecFail(iCell, :)   = cellSurWwEncRecFail;

        % rank of empirical WW-test F-value in surrogate dataset
        encRecRank(iCell, 1)        = sum(wwEncRec(iCell, 1) > cellSurWwEncRec) / param.nSur;
        encRecRankSucc(iCell, 1)    = sum(wwEncRecSucc(iCell, 1) > cellSurWwEncRecSucc) / param.nSur;
        encRecRankFail(iCell, 1)    = sum(wwEncRecFail(iCell, 1) > cellSurWwEncRecFail) / param.nSur;

        %% plot distribution of phase shifts

        % phase shift examples for figure
        examples = [5, 0, 32, 2; ...
            5, 0, 46, 2; ...
            10, 1, 34, 2; ...
            10, 1, 40, 4];

        % plot distribution
        if ismember(phaseRes(iCell).idx, examples, 'rows')  && ...
                strcmp(param.powLvl, 'all') && ...
                strcmp(param.slopeLvl, 'all') && ...
                strcmp(param.thetaLvl, 'all') && ...
                strcmp(param.unitType, 'all') && ...
                strcmp(param.cellType, 'all')

            % p-value
            wwP         = 1 - (sum(wwEncRec(iCell, 1) > surWwEncRec(iCell, :)) / param.nSur);

            % correct low p-values of permutation tests (see Phipson and Smyth, 2010)
            wwP(wwP < (1 / param.nSur)) = 1 / param.nSur;

            % plot
            permFig     = figure;
            permDistr   = histogram(surWwEncRec(iCell, :), 'FaceColor', [0.5, 0.5, 0.5]);
            hold on;
            empVal      = xline(wwEncRec(iCell, 1), 'k');
            title({strcat('{Permutation test (P = }', num2str(wwP, 3), ')')});
            xlabel('f-value Watson-Williams test');
            ylabel('Number of surrogates');
            box off;
            set(gca, 'tickDir', 'out');

            % print figure
            unitName    = strjoin(string(phaseRes(iCell).idx), '_');
            figureName  = strcat(unitName, '_SurrogateDistribution');
            set(permFig, 'renderer', 'painters');
            saveas(permFig, fullfile(paths.save, strcat(figureName, '.svg')), 'svg');
        end

        % move figures of significant cells to separate folder
        if encRecRank(iCell, 1) > 0.95
            figIdx = strcmp(allFigureNames, strcat(regexprep(num2str(phaseRes(iCell).idx),'\s+','_'), '_both_PhaseAngleFigure.jpg'));
            copyfile(fullfile(figureDir(figIdx).folder, figureDir(figIdx).name), fullfile(paths.save, 'ShiftUnitFigures', figureDir(figIdx).name));
        elseif phaseRes(iCell).encBothRank > 0.95 && phaseRes(iCell).recBothRank > 0.95
            figIdx = strcmp(allFigureNames, strcat(regexprep(num2str(phaseRes(iCell).idx),'\s+','_'), '_both_PhaseAngleFigure.jpg'));
            copyfile(fullfile(figureDir(figIdx).folder, figureDir(figIdx).name), fullfile(paths.save, 'PhaseLockingUnitFigures', figureDir(figIdx).name));
        end
    catch
        exclusionIndex(iCell, 1) = true;
        continue;
    end
end

%% find cells with significant phase shift
if strcmp(param.powLvl, 'all') && ...
        strcmp(param.slopeLvl, 'all') && ...
        strcmp(param.thetaLvl, 'all') && ...
        strcmp(param.unitType, 'all') && ...
        strcmp(param.cellType, 'all')

    %% identify cells with significant phase locking and phase shifts

    % significant phase locking - encoding
    bEncPL              = cat(1, phaseRes.encBothRank) > 0.95;
    bEncSuccPL          = cat(1, phaseRes.encSuccRank) > 0.95;
    bEncFailPL          = cat(1, phaseRes.encFailRank) > 0.95;

    % significant phase locking - recall
    bRecPL              = cat(1, phaseRes.recBothRank) > 0.95;
    bRecSuccPL          = cat(1, phaseRes.recSuccRank) > 0.95;
    bRecFailPL          = cat(1, phaseRes.recFailRank) > 0.95;

    % phase locking during encoding and recall
    bBothPL             = bEncPL & bRecPL;
    bSuccPL             = bEncSuccPL & bRecSuccPL;
    bFailPL             = bEncFailPL & bRecFailPL;

    % phase difference
    bShift              = encRecRank > 0.95;
    bShiftSucc          = encRecRankSucc > 0.95;
    bShiftFail          = encRecRankFail > 0.95;

    % phase locking and phase difference
    bPlShift            = bBothPL & bShift;
    bPlShiftSucc        = bSuccPL & bShiftSucc;
    bPlShiftFail        = bFailPL & bShiftFail;

    % index only for phase locking neurons
    bOnlyPlShift        = bShift(bBothPL);
    bOnlyPlShiftSucc    = bShiftSucc(bSuccPL);
    bOnlyPlShiftFail    = bShiftFail(bFailPL);

    % correlation between WW-test ranks of succ and unsucc memory performance
    [succFailShiftRho, succFailShiftPval]       = corr(encRecRankSucc, encRecRankFail);

    % binomial tests
    allPval             = myBinomTest(sum(bShift), size(bShift, 1), 0.05, 'one');
    shiftPval           = myBinomTest(sum(bPlShift), sum(bBothPL), 0.05, 'one');
    succPval            = myBinomTest(sum(bShiftSucc), size(encRecRankSucc, 1), 0.05, 'one') * 2; % Bonferroni corrected
    succShiftPval       = myBinomTest(sum(bPlShiftSucc), sum(bSuccPL), 0.05, 'one') * 2;
    failPval            = myBinomTest(sum(bShiftFail), size(encRecRankFail, 1), 0.05, 'one') * 2;
    failShiftPval       = myBinomTest(sum(bPlShiftFail), sum(bFailPL), 0.05, 'one') * 2;
    
    %% permutation test for differences between successful and unsuccessful performance
    
    % difference in shifting units
    shiftDiff           = (sum(bShiftSucc) / size(bShiftSucc, 1)) - (sum(bShiftFail) / size(bShiftSucc, 1));
    shiftPlDiff         = (sum(bPlShiftSucc) / size(bOnlyPlShiftSucc, 1)) - (sum(bPlShiftFail) / size(bOnlyPlShiftFail, 1));
    
    % class variable
    succFailClass       = [ones(size(bShiftSucc)); ones(size(bShiftFail)) * 2];
    succFailPlClass     = [ones(size(bOnlyPlShiftSucc)); ones(size(bOnlyPlShiftFail)) * 2];
    
    % create surrogates
    surShiftDiff        = nan(param.nSur, 1);
    surShiftPlDiff      = nan(param.nSur, 1);
    for iSur = 1:param.nSur
        
        % concatenate successful and unsuccessful logicals
        catSuccFail             = cat(1, bShiftSucc, bShiftFail);
        catPlSuccFail           = cat(1, bOnlyPlShiftSucc, bOnlyPlShiftFail);

        % randomly permute class variable
        surSuccFailClass        = succFailClass(randperm(size(catSuccFail, 1)));
        surPlSuccFailClass      = succFailPlClass(randperm(size(catPlSuccFail, 1)));

        % assign shift indices to surrogate class variables
        surSucc                 = catSuccFail(surSuccFailClass == 1);
        surFail                 = catSuccFail(surSuccFailClass == 2);
        surPlSucc               = catPlSuccFail(surPlSuccFailClass == 1);
        surPlFail               = catPlSuccFail(surPlSuccFailClass == 2);

        % difference in shifting units
        surShiftDiff(iSur, 1)   = (sum(surSucc) / size(bShiftSucc, 1)) - (sum(surFail) / size(bShiftSucc, 1));
        surShiftPlDiff(iSur, 1) = (sum(surPlSucc) / size(bOnlyPlShiftSucc, 1)) - (sum(surPlFail) / size(bOnlyPlShiftFail, 1));
    end

    % rank
    permShiftRank       = sum(shiftDiff > surShiftDiff) / param.nSur;
    permShiftPlRank     = sum(shiftPlDiff > surShiftPlDiff) / param.nSur;

    % p-value
    permShiftPval       = 1 - permShiftRank;
    permShiftPlPval     = 1 - permShiftPlRank;

    %% figure of mean angle differences of cells with significant phase difference

    % all
    circDistSig         = figure;
    circDiff            = angdiff(encAngle(bShift), recAngle(bShift));
    circDiffHist        = polarhistogram(circDiff, param.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
    hold on;
    circDiffSig         = angdiff(encAngle(bPlShift), recAngle(bPlShift));
    circDiffSigHist     = polarhistogram(circDiffSig, param.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16]);
    rlim([0, 10]);
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    title(strcat('Mean phase difference', 32, '(', num2str(sum(bPlShift)), 32, 'out of', 32, num2str(sum(bBothPL)), ')'));

    % print figure
    print(circDistSig, fullfile(paths.save, strcat(param.resultName, '_BothPhaseDifference')), '-dsvg', '-r300');

    % successful memory performance
    circDistSigSucc     = figure;
    circDiffSucc        = angdiff(encSuccAngle(bShiftSucc), recSuccAngle(bShiftSucc));
    circDiffHistSucc    = polarhistogram(circDiffSucc, param.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1], 'FaceAlpha', 0.2);
    hold on;
    circDiffSigSucc     = angdiff(encSuccAngle(bPlShiftSucc), recSuccAngle(bPlShiftSucc));
    circDiffSigHistSucc = polarhistogram(circDiffSigSucc, param.polarHistEdges, 'FaceColor', [0.1, 0.6, 0.1]);
    rlim([0, 10]);
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    title(strcat('Mean phase difference', 32, '(', num2str(sum(bPlShiftSucc)), 32, 'out of', 32, num2str(sum(bSuccPL)), ')'));

    % print figure
    print(circDistSigSucc, fullfile(paths.save, strcat(param.resultName, '_SuccPhaseDifference')), '-dsvg', '-r300');

    % unsuccessful memory performance
    circDistSigFail     = figure;
    circDiffFail        = angdiff(encFailAngle(bShiftFail), recFailAngle(bShiftFail));
    circDiffHistFail    = polarhistogram(circDiffFail, param.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1], 'FaceAlpha', 0.2);
    hold on;
    circDiffSigFail     = angdiff(encFailAngle(bPlShiftFail), recFailAngle(bPlShiftFail));
    circDiffSigHistFail = polarhistogram(circDiffSigFail, param.polarHistEdges, 'FaceColor', [0.6, 0.1, 0.1]);
    rlim([0, 10]);
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    title(strcat('Mean phase difference', 32, '(', num2str(sum(bPlShiftFail)), 32, 'out of', 32, num2str(sum(bFailPL)), ')'));

    % print figure
    print(circDistSigFail, fullfile(paths.save, strcat(param.resultName, '_FailPhaseDifference')), '-dsvg', '-r300');

    %% circular mean differences

    % circular means
    meanCircDiff        = rad2deg(circ_mean(abs(circDiff)));
    meanCircDiffSig     = rad2deg(circ_mean(abs(circDiffSig)));
    meanCircDiffSucc    = rad2deg(circ_mean(abs(circDiffSucc)));
    meanCircDiffSigSucc = rad2deg(circ_mean(abs(circDiffSigSucc)));
    meanCircDiffFail    = rad2deg(circ_mean(abs(circDiffFail)));
    meanCircDiffSigFail = rad2deg(circ_mean(abs(circDiffSigFail)));
    
    % standard deviations
    stdCircDiff         = rad2deg(circ_std(abs(circDiff)));
    stdCircDiffSig      = rad2deg(circ_std(abs(circDiffSig)));
    stdCircDiffSucc     = rad2deg(circ_std(abs(circDiffSucc)));
    stdCircDiffSigSucc  = rad2deg(circ_std(abs(circDiffSigSucc)));
    stdCircDiffFail     = rad2deg(circ_std(abs(circDiffFail)));
    stdCircDiffSigFail  = rad2deg(circ_std(abs(circDiffSigFail)));

    %% bar plot

    % define figure for the bar plot
    shiftFig        = figure;

    % define categories for the x-axis
    shiftBounds     = categorical({'all shifts', 'PL shifts', 'succ shifts', 'succ PL shifts', 'fail shifts', 'fail PL shifts'});

    % define shiftBounds with explicit order
    shiftBounds     = categorical(shiftBounds, shiftBounds);

    % reorder the categories according to the desired sequence
    desiredOrder    = {'all shifts', 'PL shifts', 'succ shifts', 'succ PL shifts', 'fail shifts', 'fail PL shifts'};
    shiftBounds     = reordercats(shiftBounds, desiredOrder);

    % calculate ratio for each category
    shiftRatio      = [
        (sum(bShift) / size(bShift, 1)), ...
        (sum(bPlShift) / sum(bBothPL)), ...
        (sum(bShiftSucc) / size(bShiftSucc, 1)), ...
        (sum(bPlShiftSucc) / sum(bSuccPL)), ...
        (sum(bShiftFail) / size(bShiftFail, 1)), ...
        (sum(bPlShiftFail) / sum(bFailPL))
        ];

    % define colors for each category
    shiftColors     = [
        0.16, 0.16, 0.16; ...
        0.16, 0.16, 0.16; ...
        0.1, 0.6, 0.1; ...
        0.1, 0.6, 0.1; ...
        0.6, 0.1, 0.1; ...
        0.6, 0.1, 0.1
        ];

    % define transparency levels for each category
    shiftAlphas     = [0.2, 0.6, 0.2, 0.6, 0.2, 0.6];

    % initialize bar objects
    b = gobjects(size(shiftBounds));

    % clear current axes and begin plotting
    cla();
    hold on;

    % loop through each category to create bars with corresponding properties
    for i = 1:numel(shiftBounds)
        b(i) = bar(shiftBounds(i), shiftRatio(i));
        set(b(i), 'FaceColor', shiftColors(i, :), 'FaceAlpha', shiftAlphas(i));
    end

    % adjustments
    ylabel('Ratio');
    ylim([0, 0.1801]);
    yticks(0:0.05:1);
    set(gca, 'TickDir', 'out');

    % print figure
    print(shiftFig, fullfile(paths.save, strcat(param.resultName, '_ShiftCellPercentage')), '-dsvg', '-r300');
else
    
    % phase locking and phase difference
    bShift              = encRecRank > 0.95;

    % number of included cells
    numInc              = sum(exclusionIndex == 0);

    % binomial tests
    allPval             = myBinomTest(sum(bShift), numInc, 0.05, 'one') * 2; % Bonferroni correction

    %% figure of mean angle differences of cells with significant phase difference

    % all
    circDistSig         = figure;
    circDiff            = angdiff(encAngle(bShift), recAngle(bShift));
    circDiffHist        = polarhistogram(circDiff, param.polarHistEdges, 'FaceColor', [0.16, 0.16, 0.16], 'FaceAlpha', 0.2);
    rlim([0, 10]);
    set(gca, 'ThetaTickLabel', param.polarHistLabels);
    title(strcat('Mean phase difference', 32, '(', num2str(sum(bShift)), 32, 'out of', 32, num2str(numInc), ')'));
    
    % print figure
    print(circDistSig, fullfile(paths.save, strcat(param.resultName, '_BothPhaseDifference')), '-dsvg', '-r300');
end
