%==========================================================================
% This script is for analyzing the relationship between the occurrence of
% theta oscillations and neuronal phase locking
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.phaseRes          = 'D:\TreasureHunt\PhaseAnalysis_20230921'; % results folder
paths.save              = 'D:\TreasureHunt\PhaseAnalysisOscillations_20230718'; % save folder

% parameters
param                   = [];
param.polarHistEdges    = linspace(-pi, pi, 21);
param.polarHistLabels   = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
param.cond              = {'rest'; 'enc'; 'rec'}; % or {'enc'; 'rec'}
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

    %% get phases and oscillation index for all spikes
    if strcmp(param.cond{iCond, 1}, 'rest')

        % get data
        allComplex      = {allResults.allRes.allComplex}';
        allEncRecIdx    = {allResults.allRes.allEncRecIdx}';
        allOscIdx       = {allResults.allRes.allThetaIdx}';

        % exclude cells
        allComplex      = allComplex(~results.bExclude);
        allEncRecIdx    = allEncRecIdx(~results.bExclude);
        allOscIdx       = allOscIdx(~results.bExclude);

        % exclude encoding and recall spikes
        allComplex      = cellfun(@(x, y) x(y == 0), allComplex, allEncRecIdx, 'Uni', 0);
        allOscIdx       = cellfun(@(x, y) x(y == 0), allOscIdx, allEncRecIdx, 'Uni', 0);

        % get phase and power data
        allPhases       = cellfun(@(x) angle(x), allComplex, 'Uni', 0);
        allPower        = cellfun(@(x) abs(x) .^ 2, allComplex, 'Uni', 0);

        % define color of histograms
        colorData       = [0.6, 0.4, 0.2];
        
    elseif strcmp(param.cond{iCond, 1}, 'enc')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encPhase}', 'Uni', 0);
        allOscIdx       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encThetaIdx}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.encSpikePower}', 'Uni', 0);
        colorData       = [0, 0.447, 0.741];
    elseif strcmp(param.cond{iCond, 1 }, 'rec')
        allPhases       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recPhase}', 'Uni', 0);
        allOscIdx       = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recThetaIdx}', 'Uni', 0);
        allPower        = cellfun(@(x) cat(1, x{:}), {results.phaseRes.recSpikePower}', 'Uni', 0);
        colorData       = [0.85, 0.325, 0.098];
    end

    %% correlation between oscillation and power

    % Spearman correlation between power and theta index
    oscPowerCorr  = cell2mat(cellfun(@(x, y) ...
        corr(x(~isnan(x) & ~isnan(y)), y(~isnan(x) & ~isnan(y)), 'Type', 'Spearman'), ...
        allOscIdx, allPower, 'Uni', 0));
    oscPowerRho   = mean(oscPowerCorr, 'omitnan');

    %% calculate zvals and differences for oscillations present or not

    % preallocate
    oscPhases       = cell(size(allPhases, 1), 1);
    nonPhases       = cell(size(allPhases, 1), 1);
    oscZval         = nan(size(allPhases, 1), 1);
    nonZval         = nan(size(allPhases, 1), 1);
    surOscZval      = nan(size(allPhases, 1), param.nSur);
    surNonZval      = nan(size(allPhases, 1), param.nSur);
    smpOscZval      = nan(size(allPhases, 1), param.nSur);
    smpNonZval      = nan(size(allPhases, 1), param.nSur);
    oscAngle        = nan(size(allPhases, 1), 1);
    nonAngle        = nan(size(allPhases, 1), 1);
    surOscAngle     = nan(size(allPhases, 1), param.nSur);
    surNonAngle     = nan(size(allPhases, 1), param.nSur);

    % loop through cells
    parfor iCell = 1:size(allPhases, 1)
        
        % set random seed
        rng(randSeedNum(iCell, 1), 'twister');
        
        % display progress
        disp(iCell);

        % select spikes by oscillation index
        oscPhases{iCell, 1}         = allPhases{iCell, 1}(allOscIdx{iCell, 1} == 1);
        nonPhases{iCell, 1}         = allPhases{iCell, 1}(allOscIdx{iCell, 1} == 0);

        % get mean and mean resultant vector length for each unit
        [~, oscAngle(iCell, 1)]     = circ_axialmean(oscPhases{iCell, 1});
        [~, nonAngle(iCell, 1)]     = circ_axialmean(nonPhases{iCell, 1});

        % get Rayleigh z-value
        try
            [~, oscZval(iCell, 1)]      = circ_rtest(oscPhases{iCell, 1});
            [~, nonZval(iCell, 1)]      = circ_rtest(nonPhases{iCell, 1});
        catch
            oscZval(iCell, 1)           = NaN;
            nonZval(iCell, 1)           = NaN;
        end

        %% create surrogates
        smpOscPhases        = [];
        smpNonPhases        = [];
        cellSurOscZval      = nan(1, param.nSur);
        cellSurNonZval      = nan(1, param.nSur);
        cellSmpOscZval      = nan(1, param.nSur);
        cellSmpNonZval      = nan(1, param.nSur);
        cellSurOscAngle     = nan(1, param.nSur);
        cellSurNonAngle     = nan(1, param.nSur);
        for iSur = 1:param.nSur
            
            %% surrogate test

            % random assignment of theta index
            surAllOscIdx                     = datasample(allOscIdx{iCell, 1}, size(allOscIdx{iCell, 1}, 1), 'Replace', false);
            
            % create surrogates
            surOscPhases                     = allPhases{iCell, 1}(surAllOscIdx == 1);
            surNonPhases                     = allPhases{iCell, 1}(surAllOscIdx == 0);

            % get mean resultant vector length and mean angle
            [~, cellSurOscAngle(1, iSur)]    = circ_axialmean(surOscPhases);
            [~, cellSurNonAngle(1, iSur)]    = circ_axialmean(surNonPhases);
            
            % Rayleigh z-value
            try
                [~, cellSurOscZval(1, iSur)] = circ_rtest(surOscPhases);
                [~, cellSurNonZval(1, iSur)] = circ_rtest(surNonPhases);
            catch
                cellSurOscZval(1, iSur)      = NaN;
                cellSurNonZval(1, iSur)      = NaN;
            end

            %% subsampling to correct for spike numbers

            % compare size of oscillation and non-oscillation phases
            if size(oscPhases{iCell, 1}, 1) > size(nonPhases{iCell, 1}, 1)
                
                % draw a subsample from the oscillation phases
                smpOscPhases = datasample(oscPhases{iCell, 1}, size(nonPhases{iCell, 1}, 1), 'Replace', false);
                smpNonPhases = nonPhases{iCell, 1};

            elseif size(oscPhases{iCell, 1}, 1) < size(nonPhases{iCell, 1}, 1)
                
                % draw a subsample from the non-oscillation phases
                smpOscPhases = oscPhases{iCell, 1};
                smpNonPhases = datasample(nonPhases{iCell, 1}, size(oscPhases{iCell, 1}, 1), 'Replace', false);
            end
            
            % Rayleigh z-value
            try
                [~, cellSmpOscZval(1, iSur)] = circ_rtest(smpOscPhases);
                [~, cellSmpNonZval(1, iSur)] = circ_rtest(smpNonPhases);
            catch
                cellSmpOscZval(1, iSur)      = NaN;
                cellSmpNonZval(1, iSur)      = NaN;
            end
        end

        % collect surrogates across cells
        surOscAngle(iCell, :)   = cellSurOscAngle;
        surNonAngle(iCell, :)   = cellSurNonAngle;
        surOscZval(iCell, :)    = cellSurOscZval;
        surNonZval(iCell, :)    = cellSurNonZval;
        smpOscZval(iCell, :)    = cellSmpOscZval;
        smpNonZval(iCell, :)    = cellSmpNonZval;
    end
    
    %% t-tests for phase locking across units

    % empirical t-tests
    [~, ~, ~, oscNonStat]       = ttest(oscZval, nonZval);
    oscNonT                     = oscNonStat.tstat;
    
    % surrogate t-tests
    [~, ~, ~, surOscNonStat]    = ttest(surOscZval, surNonZval);
    surOscNonT                  = surOscNonStat.tstat;

    % subsampling t-tests
    [~, smpOscNonPval, ~, ~]    = ttest(mean(smpOscZval, 2, 'omitnan'), mean(smpNonZval, 2, 'omitnan'));

    %% Watson-William test

    % inclusion index
    incIdx              = ~isnan(oscAngle) & ~isnan(nonAngle);

    % mean difference
    meanDiff            = circ_mean(angdiff(oscAngle(incIdx), nonAngle(incIdx)));

    % empirical Watson-William test
    [~, wwTable]        = circ_wwtest(oscAngle(incIdx), nonAngle(incIdx));
    oscNonWw            = wwTable{2, 5};

    % surrogate Watson-William test
    surOscNonWw = nan(param.nSur, 1);
    for iSur = 1:param.nSur
        
        % surrogate inclusion index
        surIncIdx               = ~isnan(surOscAngle(:, iSur)) & ~isnan(surNonAngle(:, iSur));

        % surrogate Watson-William test
        [~, surWwTable]         = circ_wwtest(surOscAngle(surIncIdx, iSur), surNonAngle(surIncIdx, iSur));
        surOscNonWw(iSur, 1)    = surWwTable{2, 5};
    end
    
    % surrogate test results
    oscNonRank         = sum(oscNonT > surOscNonT) / param.nSur;
    oscNonWwRank       = sum(oscNonWw > surOscNonWw) / param.nSur;
    oscNonPval         = (1 - oscNonRank) * 3; % Bonferroni correction
    smpOscNonPval      = smpOscNonPval * 3;
    oscNonWwPval       = (1 - oscNonWwRank) * 3;
    
    % plot histogram of surrogate and empirical t-value
    surFigure          = figure('units', 'centimeters', 'Position', [10, 10, 6, 4]);
    surHist            = histogram(surOscNonT, 'EdgeColor', 'none', 'FaceColor', [0.4, 0.4, 0.4]);
    xlim([3, 5]);
    ylim([0, 800]);
    xline(oscNonT, 'color', colorData, 'LineWidth', 1);
    set(gca,'TickDir','out');
    box off;
    saveas(surFigure, fullfile(paths.save, strcat(param.cond{iCond}, '_surDistrZval.svg')));
    
    %% histogram to compare zvals of theta and no theta
    
    % calculate mean and standard error of the mean
    meanZval    = mean([mean(smpOscZval, 2, 'omitnan'), mean(smpNonZval, 2, 'omitnan')], 'omitnan');
    sem         = std([mean(smpOscZval, 2, 'omitnan'), mean(smpNonZval, 2, 'omitnan')], 'omitnan') / sqrt(length([oscZval(incIdx), nonZval(incIdx)]));

    % create figure
    zvalFig     = figure;
    zvalMean    = bar(meanZval, 'FaceColor', colorData);
    hold on;
    zvalErr     = errorbar(meanZval, sem, 'LineStyle', 'none', 'Color', 'k');
    set(gca, 'XTickLabel', {'theta', 'no theta'});
    set(gca,'TickDir', 'out');
    box off;
    if contains('rest', param.cond{iCond})
            ylim([0, 40]);
    elseif contains('enc', param.cond{iCond})
            ylim([0, 3]);
    elseif contains('rec', param.cond{iCond})
            ylim([0, 8]);
    end
    title(strcat(param.cond{iCond}, 32, 'P =', 32, num2str(oscNonPval)));
    saveas(zvalFig, fullfile(paths.save, strcat(param.cond{iCond}, '_meanZval.svg')));
end

