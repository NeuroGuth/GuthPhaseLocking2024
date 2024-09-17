%==========================================================================
% This script is for analyzing the relationship between fooof slope
% and neuronal phase locking
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisAperiodic_20230703'; % save folder

% parameters
param                   = [];
param.polarHistEdges    = linspace(-pi, pi, 21);
param.polarHistLabels   = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
param.cond              = {'rest'; 'enc'; 'rec'}; % {'rest'; 'enc'; 'rec'}; % or {'enc'; 'rec'};
param.nSur              = 10001; % 10001

% set random seed
randSeed    = 444;
rng(randSeed, 'twister');
randSeedNum = randi(100000, 100000, 1);

% own functions
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% load phase results data
results         = load(fullfile(paths.phaseRes, 'additionalResultsSmall.mat'));
if sum(strcmp(param.cond, 'rest')) > 0
    allResults  = load(fullfile(paths.phaseRes, 'allResults.mat'));
end

%% loop through encoding and recall
for iCond = 1:size(param.cond, 1)

    %% get phases and slope index for all spikes
    if strcmp(param.cond{iCond, 1}, 'rest')

        % get data
        allComplex      = {allResults.allRes.allComplex}';
        allEncRecIdx    = {allResults.allRes.allEncRecIdx}';
        allSlope        = {allResults.allRes.allSpikeSlope}';
        allOscIdx       = {allResults.allRes.allThetaIdx}';

        % exclude cells
        allComplex      = allComplex(~results.bExclude);
        allEncRecIdx    = allEncRecIdx(~results.bExclude);
        allSlope        = allSlope(~results.bExclude);
        allOscIdx       = allOscIdx(~results.bExclude);

        % exclude encoding and recall spikes
        allComplex      = cellfun(@(x, y) x(y == 0), allComplex, allEncRecIdx, 'Uni', 0);
        allSlope        = cellfun(@(x, y) x(y == 0), allSlope, allEncRecIdx, 'Uni', 0);
        allOscIdx       = cellfun(@(x, y) x(y == 0), allOscIdx, allEncRecIdx, 'Uni', 0);

        % get phase and power data
        allPhases       = cellfun(@(x) angle(x), allComplex, 'Uni', 0);
        allPower        = cellfun(@(x) abs(x) .^ 2, allComplex, 'Uni', 0);

        % define color of histograms
        colorData       = [0.6, 0.4, 0.2];

    elseif strcmp(param.cond{iCond, 1}, 'enc')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encPhase}', 'Uni', 0);
        allSlope        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikeSlope}', 'Uni', 0);
        allOscIdx       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encThetaIdx}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikePower}', 'Uni', 0);
        colorData       = [0, 0.447, 0.741];
    elseif strcmp(param.cond{iCond, 1}, 'rec')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recPhase}', 'Uni', 0);
        allSlope        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikeSlope}', 'Uni', 0);
        allOscIdx       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recThetaIdx}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikePower}', 'Uni', 0);
        colorData       = [0.85, 0.325, 0.098];
    end

    %% calculate slope index and power index
    allSlopeIdx     = cellfun(@(x) double(x > median(x, 'omitnan')), allSlope, 'Uni', 0);
    allPowerIdx     = cellfun(@(x) double(x > median(x, 'omitnan')), allPower, 'Uni', 0);

    %% set slope index to NaN were slope results had NaNs
    for iCell = 1:size(allPhases, 1)
        
        % get index of NaNs
        slopeNanIdx                         = isnan(allSlope{iCell, 1});

        % adjust slope index
        allSlopeIdx{iCell, 1}(slopeNanIdx)  = NaN;
    end

    %% correlation

    % spearman correlation between slope and power
    slopePowCorr    = cell2mat(cellfun(@(x, y) ...
        corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)), 'Type', 'Spearman'), ...
        allSlope, allPower, 'Uni', 0));
    slopePowRho     = mean(slopePowCorr, 'omitnan');

    % spearman correlation between slope and theta index
    slopeThetaCorr  = cell2mat(cellfun(@(x, y) ...
        corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)), 'Type', 'Spearman'), ...
        allSlope, allOscIdx, 'Uni', 0));
    slopeThetaRho   = mean(slopeThetaCorr, 'omitnan');

    % spearman correlation between theta index and power
    oscPowerCorr    = cell2mat(cellfun(@(x, y) ...
        corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)), 'Type', 'Spearman'), ...
        allOscIdx, allPower, 'Uni', 0));
    oscPowerRho     = mean(oscPowerCorr, 'omitnan');

    %% calculate zvals and differences between different slopes

    % preallocate
    topPhases       = cell(size(allPhases, 1), 1);
    lowPhases       = cell(size(allPhases, 1), 1);
    topZval         = nan(size(allPhases, 1), 1);
    lowZval         = nan(size(allPhases, 1), 1);
    surTopZval      = nan(size(allPhases, 1), param.nSur);
    surLowZval      = nan(size(allPhases, 1), param.nSur);
    topAngle        = nan(size(allPhases, 1), 1);
    lowAngle        = nan(size(allPhases, 1), 1);
    surTopAngle     = nan(size(allPhases, 1), param.nSur);
    surLowAngle     = nan(size(allPhases, 1), param.nSur);

    % loop through cells
    parfor iCell = 1:size(allPhases, 1)

        % set random seed
        rng(randSeedNum(iCell, 1), 'twister');

        % display progress
        disp(iCell);

        % get slope index
        thisCellSlopeIdx         = allSlopeIdx{iCell, 1};

        % select spikes by slope index
        topPhases{iCell, 1}      = allPhases{iCell, 1}(thisCellSlopeIdx == 1);
        lowPhases{iCell, 1}      = allPhases{iCell, 1}(thisCellSlopeIdx == 0);

        % get mean and mean resultant vector length for each unit
        [~, topAngle(iCell, 1)]  = circ_axialmean(topPhases{iCell, 1});
        [~, lowAngle(iCell, 1)]  = circ_axialmean(lowPhases{iCell, 1});

        % get Rayleigh z-value (n * (r ^ 2))
        try
            [~, topZval(iCell, 1)]   = circ_rtest(topPhases{iCell, 1});
            [~, lowZval(iCell, 1)]   = circ_rtest(lowPhases{iCell, 1});
        catch
            topZval(iCell, 1)        = NaN;
            lowZval(iCell, 1)        = NaN;
            cellSurTopZval           = nan(1, param.nSur);
            cellSurLowZval           = nan(1, param.nSur);
            cellSurTopAngle          = nan(1, param.nSur);
            cellSurLowAngle          = nan(1, param.nSur);
            continue;
        end

        %% create surrogates
        cellSurTopZval      = nan(1, param.nSur);
        cellSurLowZval      = nan(1, param.nSur);
        cellSurTopAngle     = nan(1, param.nSur);
        cellSurLowAngle     = nan(1, param.nSur);
        for iSur = 1:param.nSur
            
            % random assignment of top or low slope
            surAllSlopeIdx  = datasample(thisCellSlopeIdx, size(thisCellSlopeIdx, 1), 'Replace', false);
            
            % create surrogates
            surTopPhases    = allPhases{iCell, 1}(surAllSlopeIdx == 1);
            surLowPhases    = allPhases{iCell, 1}(surAllSlopeIdx == 0);
            
            % get mean resultant vector length and mean angle
            [~, cellSurTopAngle(1, iSur)]  = circ_axialmean(cat(1, surTopPhases));
            [~, cellSurLowAngle(1, iSur)]  = circ_axialmean(cat(1, surLowPhases));

            % get Rayleigh z-value
            [~, cellSurTopZval(1, iSur)]   = circ_rtest(cat(1, surTopPhases));
            [~, cellSurLowZval(1, iSur)]   = circ_rtest(cat(1, surLowPhases));
        end

        % collect surrogates across cells
        surTopZval(iCell, :)    = cellSurTopZval;
        surLowZval(iCell, :)    = cellSurLowZval;
        surTopAngle(iCell, :)   = cellSurTopAngle;
        surLowAngle(iCell, :)   = cellSurLowAngle;
    end

    %% t-tests for phase locking across units

    % empirical t-tests
    [~, ~, ~, topLowStat]       = ttest(topZval, lowZval);
    topLowT                     = topLowStat.tstat;

    % surrogate t-tests
    [~, ~, ~, surTopLowStat]    = ttest(surTopZval, surLowZval);
    surTopLowT                  = surTopLowStat.tstat;

    %% Watson-William test

    % inclusion index
    incIdx              = ~isnan(topAngle) & ~isnan(lowAngle);

    % mean difference
    meanDiff            = circ_mean(angdiff(topAngle(incIdx), lowAngle(incIdx)));

    % empirical Watson-William test
    [~, wwTable]        = circ_wwtest(topAngle(incIdx), lowAngle(incIdx));
    topLowWw            = wwTable{2, 5};

    % surrogate Watson-William test
    surTopLowWw = nan(param.nSur, 1);
    for iSur = 1:param.nSur

        % surrogate inclusion index
        surIncIdx               = ~isnan(surTopAngle(:, iSur)) & ~isnan(surLowAngle(:, iSur));
        
        % surrogate Watson-William test
        [~, surWwTable]         = circ_wwtest(surTopAngle(surIncIdx, iSur), surLowAngle(surIncIdx, iSur));
        surTopLowWw(iSur, 1)    = surWwTable{2, 5};
    end

    % surrogate test results
    topLowRank         = sum(topLowT > surTopLowT) / param.nSur;
    topLowWwRank       = sum(topLowWw > surTopLowWw) / param.nSur;
    topLowPval         = (1 - topLowRank) * 3; % Bonferroni correction
    topLowWwPval       = (1 - topLowWwRank) * 3;
    
    % plot histogram of surrogate and empirical t-value
    surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    surHist            = histogram(surTopLowT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xline(topLowT, 'color', colorData, 'LineWidth', 1);
    set(gca,'TickDir','out');
    box off;
    xlim([-10, 10]);
    ylim([0, 1000]);
    saveas(surFigure, fullfile(paths.save, strcat(param.cond{iCond}, '_surDistrZval .svg')));

    %% histogram to compare zvals of different slopes

    % calculate mean and standard error of the mean
    meanZval    = mean([topZval, lowZval], 'omitnan');
    sem         = std([topZval, lowZval], 'omitnan') / sqrt(length([topZval, lowZval]));

    % create figure
    zvalFig     = figure;
    zvalMean    = bar(meanZval, 'FaceColor', colorData);
    hold on;
    zvalErr     = errorbar(meanZval, sem, 'LineStyle', 'none', 'Color', 'k');
    set(gca, 'XTickLabel', {'high', 'low'});
    set(gca,'TickDir','out');
    box off;
    % ylim([0, 0.25]);
    title(strcat(param.cond{iCond}, 32, 'P =', 32, num2str(topLowPval)));
    saveas(zvalFig, fullfile(paths.save, strcat(param.cond{iCond}, '_meanZval.svg')));
end
