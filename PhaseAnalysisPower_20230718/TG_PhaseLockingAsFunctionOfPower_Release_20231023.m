%==========================================================================
% This script analyzes the relationship between power and
% neuronal phase locking
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisPower_20230718'; % save folder

% parameters
param                   = [];
param.polarHistEdges    = linspace(-pi, pi, 21);
param.polarHistLabels   = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
param.cond              = {'baseline'; 'enc'; 'rec'};
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
if sum(strcmp(param.cond, 'baseline')) > 0
    allResults  = load(fullfile(paths.phaseRes, 'allResults.mat'));
end

%% loop through encoding and recall
for iCond = 1:size(param.cond, 1)
    
    %% get phases and power index for all spikes
    if strcmp(param.cond{iCond, 1}, 'baseline')
        
        % get data
        allComplex      = {allResults.allRes.allComplex}';
        allEncRecIdx    = {allResults.allRes.allEncRecIdx}';
        
        % exclude cells
        allComplex      = allComplex(~results.bExclude);
        allEncRecIdx    = allEncRecIdx(~results.bExclude);

        % exclude encoding and recall spikes
        allComplex      = cellfun(@(x, y) x(y == 0), allComplex, allEncRecIdx, 'Uni', 0);
        
        % get phase and power data
        allPhases       = cellfun(@(x) angle(x), allComplex, 'Uni', 0);
        allPower        = cellfun(@(x) abs(x) .^ 2, allComplex, 'Uni', 0);

        % define color of histograms
        colorData       = [0.6, 0.4, 0.2];

    elseif strcmp(param.cond{iCond, 1}, 'enc')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encPhase}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikePower}', 'Uni', 0);
        colorData       = [0, 0.447, 0.741];
    elseif strcmp(param.cond{iCond, 1}, 'rec')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recPhase}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikePower}', 'Uni', 0);
        colorData       = [0.85, 0.325, 0.098];
    end

    %% calculate power index
    allPowerIdx     = cellfun(@(x) x > median(x), allPower, 'Uni', 0);

    %% calculate zvals and differences for oscillations present or not

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

        % select spikes by power index
        topPhases{iCell, 1}     = allPhases{iCell, 1}(allPowerIdx{iCell, 1} == 1);
        lowPhases{iCell, 1}     = allPhases{iCell, 1}(allPowerIdx{iCell, 1} == 0);
        
        % get mean and mean resultant vector length for each unit
        [~, topAngle(iCell, 1)] = circ_axialmean(topPhases{iCell, 1});
        [~, lowAngle(iCell, 1)] = circ_axialmean(lowPhases{iCell, 1});

        % get Rayleigh z-value (n * (r ^ 2))
        [~, topZval(iCell, 1)]  = circ_rtest(topPhases{iCell, 1});
        [~, lowZval(iCell, 1)]  = circ_rtest(lowPhases{iCell, 1});

        %% create surrogates
        cellSurTopZval      = nan(1, param.nSur);
        cellSurLowZval      = nan(1, param.nSur);
        cellSurTopAngle     = nan(1, param.nSur);
        cellSurLowAngle     = nan(1, param.nSur);
        for iSur = 1:param.nSur
            
            % random assignment of top or low slope
            surAllPowerIdx                  = datasample(allPowerIdx{iCell, 1}, size(allPowerIdx{iCell, 1}, 1), 'Replace', false);

            % create surrogates
            surTopPhases                    = allPhases{iCell, 1}(surAllPowerIdx == 1);
            surLowPhases                    = allPhases{iCell, 1}(surAllPowerIdx == 0);

            % get mean resultant vector length and mean angle
            [~, cellSurTopAngle(1, iSur)]   = circ_axialmean(cat(1, surTopPhases));
            [~, cellSurLowAngle(1, iSur)]   = circ_axialmean(cat(1, surLowPhases));

            % get Rayleigh z-value
            [~, cellSurTopZval(1, iSur)]    = circ_rtest(cat(1, surTopPhases));
            [~, cellSurLowZval(1, iSur)]    = circ_rtest(cat(1, surLowPhases));
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

    % mean angular difference
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
    % limsx              = get(gca, 'XLim');
    % set(gca,'Xlim', [param.surHistEdges(1), limsx(2)]);
    xlim([-15, 15]);
    ylim([0, 1000]);
    set(gca,'TickDir','out');
    box off;
    saveas(surFigure, fullfile(paths.save, strcat(param.cond{iCond}, '_surDistr.svg')));

    %% histogram to compare z-values of different powers

    % calculate mean and standard error of the mean
    meanZval    = mean([topZval, lowZval], 'omitnan');
    sem         = std([topZval, lowZval], 'omitnan') / sqrt(length([topZval, lowZval]));

    % create figure
    zvalFig     = figure;
    zvalMean    = bar(meanZval, 'FaceColor', colorData);
    hold on;
    zvalErr     = errorbar(meanZval, sem, 'LineStyle', 'none', 'Color', 'k');
    set(gca, 'XTickLabel', {'high', 'low'});
    if strcmp(param.cond{iCond, 1}, 'baseline')
        ylim([0, 200]);
    elseif strcmp(param.cond{iCond, 1}, 'enc')
        ylim([0, 10]);
    elseif strcmp(param.cond{iCond, 1}, 'rec')
        ylim([0, 40]);
    end
    set(gca,'TickDir','out');
    box off;
    title(strcat(param.cond{iCond}, 32, 'P =', 32, num2str(topLowPval)));
    saveas(zvalFig, fullfile(paths.save, strcat(param.cond{iCond}, '_meanMRVL.svg')));
end
