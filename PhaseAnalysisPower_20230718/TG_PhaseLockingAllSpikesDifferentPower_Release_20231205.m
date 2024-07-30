%==========================================================================
% This script creates figures of phase locking of spikes at different power
% levels
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% paths
paths.load              = 'D:\TreasureHunt\PhaseAnalysis_20230921';
paths.save              = 'D:\TreasureHunt\PhaseAnalysisPower_20230718'; % save folder
if ~exist(paths.save, 'dir')
    mkdir(paths.save);
end

% own functions
addpath(genpath('D:\External\Functions'));

% add fieltrip toolbox
addpath('D:\External\Toolboxes\Fieldtrip');
ft_defaults; % default fieldtrip settings

% set random seed
randSeed                = 444;
rng(randSeed, 'twister');
randSeedNum             = randi(100000, 100000, 1); % for randomseed

% parameters
param                   = [];
param.polarHistEdges    = linspace(-pi, pi, 21);
param.polarHistLabels   = {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'};
param.nSur              = 10001;

%% load previously saved phase analysis results
results                 = load(fullfile(paths.load, 'additionalResultsSmall.mat'));
allResults              = load(fullfile(paths.load, 'allResults.mat'), 'allRes');
allResults              = allResults.allRes;
fprintf('Total number of cells: %d.\n', size(allResults, 1));

% get spike-specific data for all cells
allComplex              = {allResults.allComplex}';

% exclude cells
allComplex              = allComplex(~results.bExclude);

% get phase data
allPhases               = cellfun(@(x) angle(x), allComplex, 'Uni', 0);
catAllPhases            = cat(1, allPhases{:});

% get power data
allPower                = cellfun(@(x) abs(x) .^ 2, allComplex, 'Uni', 0);

% calculate power index
allPowerIdx             = cellfun(@(x) x > median(x), allPower, 'Uni', 0);
catAllPowerIdx          = cat(1, allPowerIdx{:});

% select spikes by power index - top phases
topPhases               = cellfun(@(x, y) x(y == 1), allPhases, allPowerIdx, 'Uni', 0);
catTopPhases            = cat(1, topPhases{:});

% select spikes by power index - low phases
lowPhases               = cellfun(@(x, y) x(y == 0), allPhases, allPowerIdx, 'Uni', 0);
catLowPhases            = cat(1, lowPhases{:});

%% polar histogram for all power
allPowFig   = figure;
pAll        = polarhistogram(catAllPhases, 60, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', param.polarHistLabels);
rlim([0, 0.025]);
legend([num2str(size(catAllPhases, 1)), newline, strcat('all-power spikes - ', 'all')], 'Location', 'South');

% print figure
print(allPowFig, fullfile(paths.save, 'FigAllPower_20231017'), '-dsvg', '-r300');

%% polar histogram for high power
topPowFig   = figure;
pTop        = polarhistogram(catTopPhases, 60, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', param.polarHistLabels);
rlim([0, 0.025]);
legend([num2str(size(catTopPhases, 1)), newline, strcat('high-power spikes - ', 'all')], 'Location', 'South');

% print figure
print(topPowFig, fullfile(paths.save, 'FigTopPower_20231017'), '-dsvg', '-r300');

%% polar histogram for low power
lowPowFig   = figure;
pLow        = polarhistogram(catLowPhases, 60, 'FaceColor', [0.4, 0.4, 0.4], 'FaceAlpha', 1, 'Normalization', 'probability');
set(gca, 'ThetaTickLabel', param.polarHistLabels);
rlim([0, 0.025]);
legend([num2str(size(catLowPhases, 1)), newline, strcat('low-power spikes - ', 'all')], 'Location', 'South');

% print figure
print(lowPowFig, fullfile(paths.save, 'FigLowPower_20231017'), '-dsvg', '-r300');

%% Rayleigh z-values

% all phases
[~, allZ]  = circ_rtest(catAllPhases);

% top phases
[~, topZ]  = circ_rtest(catTopPhases);

% low phases
[~, lowZ]  = circ_rtest(catLowPhases);

%% difference in z-values between high and low power

% empirical z-value difference between top and low power
zDiff      = topZ - lowZ;

% create surrogate z-value differences
surZDiff = nan(1, param.nSur);
parfor iSur = 1:param.nSur
    
    % set random seed
    rng(randSeedNum(iSur, 1), 'twister');

    % display progress
    disp(iSur);

    % create surrogate power index
    surPowerIdx         = datasample(catAllPowerIdx, size(catAllPowerIdx, 1), 'Replace', false);
    
    % create surrogate top and low power phases
    surTopPhases        = catAllPhases(surPowerIdx == 1);
    surLowPhases        = catAllPhases(surPowerIdx == 0);

    % calculate surrogate z-values
    [~, surTopZ]        = circ_rtest(surTopPhases);
    [~, surLowZ]        = circ_rtest(surLowPhases);

    % surrogate z-value difference between top and low power
    surZDiff(1, iSur)   = surTopZ - surLowZ;
end

