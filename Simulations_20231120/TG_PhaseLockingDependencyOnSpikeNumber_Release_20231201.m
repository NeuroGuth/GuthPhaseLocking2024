%==========================================================================
% This script tests whether the MRVL, its rank in a distribution
% of surrogate MRVLs and a calculated z-value depend on the number of
% spikes
%
% Tim Guth, 2023
%==========================================================================

% start
clear; close all; clc;

% add paths
paths.save      = 'D:\TreasureHunt\Simulations_20231120';

% add functions
addpath(genpath('D:\External\Functions\'));

% settings
spikeNum        = [10; 100; 1000];
cellNum         = 1000;
nSur            = 1000;
randSeedNum     = 444; % set random seed
rng(randSeedNum, 'twister');

% preallocate
zval            = nan(size(spikeNum, 1), cellNum);
subZ            = nan(1, cellNum);
mrvl            = nan(size(spikeNum, 1), cellNum);
rank            = nan(size(spikeNum, 1), cellNum);
zscored         = nan(size(spikeNum, 1), cellNum);

% loop through example cells
for iC = 1:cellNum
    
    % display progress
    disp(iC);

    % random phase-locking neuron
    randPhasesPL = cat(1, vmrand(0, 2, [250, 1]), 2 * pi * rand(750, 1) - pi);

    % subsampling 1000 spikes to 100 spikes
    cellSubZ = nan(nSur, 1);
    for iSp = 1:nSur

        % draw a 100 spikes subsample from the 1000 spikes
        subPhases               = datasample(randPhasesPL, spikeNum(2), 'Replace', false);
        [~, cellSubZ(iSp, 1)]   = circ_rtest(subPhases);
    end
    
    % calculate mean of the subsamples
    subZ(1, iC) = mean(cellSubZ);

    % loop through sample sizes
    for iSn = 1:size(spikeNum, 1)

        % draw n random spikes of neuron
        smpPhases            = datasample(randPhasesPL, spikeNum(iSn), 'Replace', false);

        % Rayleigh z-value
        [~, zval(iSn, iC)]   = circ_rtest(smpPhases);

        % mean resultant vector length
        mrvl(iSn, iC)        = round(circ_axialmean(smpPhases), 5);

        % surrogates
        surDistr             = 2 * pi * rand(spikeNum(iSn), nSur) - pi;
        surMrvl              = circ_axialmean(surDistr);

        % rank of mrvl
        rank(iSn, iC)        = sum(mrvl(iSn, iC) > surMrvl) / nSur;

        % z-score mrvl
        zscored(iSn, iC)     = (mrvl(iSn, iC) - mean(surMrvl)) / std(surMrvl);
    end
end

%% random phase-locking neurons

% 10 spikes
f1 = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
p1 = polarhistogram(randPhasesPL(1:spikeNum(1)), 24, 'FaceColor', [0.2, 0.2, 1]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
rlim([0, 2]);
print(f1, fullfile(paths.save, '10Spikes'), '-dsvg', '-r300');

% 100 spikes
f2 = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
p2 = polarhistogram(randPhasesPL(1:spikeNum(2)), 24, 'FaceColor', [0.2, 1, 0.2]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
print(f2, fullfile(paths.save, '100Spikes'), '-dsvg', '-r300');

% 1000 spikes
f3 = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
p3 = polarhistogram(randPhasesPL(1:spikeNum(3)), 24, 'FaceColor', [1, 0.2, 0.2]);
set(gca, 'ThetaTickLabel', {'0'; '30'; '60'; '90'; '120'; '150'; '\pm180'; '-150'; '-120'; '-90'; '-60'; '-30'});
print(f3, fullfile(paths.save, '1000Spikes'), '-dsvg', '-r300');

%% histograms for mrvl
mF = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% mrvl distributions
histogram(mrvl(1, :), 'FaceColor',[0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(mrvl(2, :), 'FaceColor', [0.2, 1, 0.2], 'EdgeColor', 'none');
histogram(mrvl(3, :), 'FaceColor', [1, 0.2, 0.2], 'EdgeColor', 'none');

% means
xline(mean(mrvl(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 3);
xline(mean(mrvl(2, :)), 'Color', [0.2, 1, 0.2], 'LineWidth', 3);
xline(mean(mrvl(3, :)), 'Color', [1, 0.2, 0.2], 'LineWidth', 3);

% settings
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of surrogates');
xlabel('MRVL');
print(mF, fullfile(paths.save, 'MRVL'), '-dsvg', '-r300');

%% histograms for Rayleigh's z
mRz = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% distributions
histogram(zval(1, :), 'FaceColor', [0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(zval(2, :), 'FaceColor', [0.2, 1, 0.2], 'EdgeColor', 'none');
histogram(zval(3, :), 'FaceColor', [1, 0.2, 0.2], 'EdgeColor', 'none');

% means
xline(mean(zval(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 3);
xline(mean(zval(2, :)), 'Color', [0.2, 1, 0.2], 'LineWidth', 3);
xline(mean(zval(3, :)), 'Color', [1, 0.2, 0.2], 'LineWidth', 3);

% settings
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of surrogates');
xlabel('Rayleigh z-value');
print(mRz, fullfile(paths.save, 'Rayleigh'), '-dsvg', '-r300');

%% histograms for rank
mRank = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% distributions
histogram(rank(1, :), 'FaceColor', [0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(rank(2, :), 'FaceColor', [0.2, 1, 0.2], 'EdgeColor', 'none');
% histogram(rank(3, :), 'FaceColor', [1, 0.2, 0.2], 'EdgeColor', 'none');

% means
xline(mean(rank(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 3);
xline(mean(rank(2, :)), 'Color', [0.2, 1, 0.2], 'LineWidth', 3);
xline(mean(rank(3, :)), 'Color', [1, 0.2, 0.2], 'LineWidth', 3);

% settings
ylim([0, nSur]);
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of surrogates');
xlabel('Rank in 1000 surrogates');
print(mRank, fullfile(paths.save, 'Rank'), '-dsvg', '-r300');

%% histograms for zval
mZval = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);

% distributions
histogram(zscored(1, :), 'FaceColor', [0.2, 0.2, 1], 'EdgeColor', 'none');
hold on;
histogram(zscored(2, :), 'FaceColor', [0.2, 1, 0.2], 'EdgeColor', 'none');
histogram(zscored(3, :), 'FaceColor', [1, 0.2, 0.2], 'EdgeColor', 'none');

% means
xline(mean(zscored(1, :)), 'Color', [0.2, 0.2, 1], 'LineWidth', 3);
xline(mean(zscored(2, :)), 'Color', [0.2, 1, 0.2], 'LineWidth', 3);
xline(mean(zscored(3, :)), 'Color', [1, 0.2, 0.2], 'LineWidth', 3);

% settings
box off;
set(gca, 'tickdir', 'out');
ylabel('Number of surrogates');
xlabel('Z-score (against 1001 surrogates)');
print(mZval, fullfile(paths.save, 'Zscore'), '-dsvg', '-r300');

%% subsampling demonstration of Rayleigh z-value

% mean and SEM of Rayleigh z-value
meanZ           = mean(zval, 2);
semZ            = std(zval, [], 2) / sqrt(size(zval, 2));

% plot Rayleigh's z-values for 100 and 1000 spikes
rayFig          = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
b1               = bar([meanZ(2), meanZ(3)], 'FaceColor', 'flat', 'FaceAlpha', 0.6);
set(gca, 'XTickLabel', {'100 spikes'; '1000 spikes'});
b1.CData(1, :)   = [0.2, 1, 0.2];
b1.CData(2, :)   = [1, 0.2, 0.2];
hold on;
e1               = errorbar([1; 2], [meanZ(2), meanZ(3)], [semZ(2), semZ(3)], '');
e1.Color         = [0, 0, 0];
e1.LineStyle     = 'none';
ylabel('Mean Rayleigh''s z');
ylim([0, 40]);
box off;
set(gca, 'tickdir', 'out');
print(rayFig, fullfile(paths.save, 'Rayleigh100vs1000'), '-dsvg', '-r300');

% mean and SEM of subsampled Rayleigh z-value
meanSubZ        = mean(subZ);
semSubZ         = std(subZ) / sqrt(size(subZ, 2));

% plot Rayleigh's z-values for 100 and 1000 spikes
raySubFig        = figure('units', 'centimeters', 'position', [15, 15, 6, 6]);
b2               = bar([meanZ(2), meanSubZ], 'FaceColor', 'flat', 'FaceAlpha', 0.6);
set(gca, 'XTickLabel', {'100 spikes'; '1000 spikes'});
b2.CData(1, :)   = [0.2, 1, 0.2];
b2.CData(2, :)   = [1, 0.2, 0.2];
hold on;
e2               = errorbar([1; 2], [meanZ(2), meanSubZ], [semZ(2), semSubZ], '');
e2.Color         = [0, 0, 0];
e2.LineStyle     = 'none';
ylabel('Mean Rayleigh''s z');
ylim([0, 5]);
box off;
set(gca, 'tickdir', 'out');
print(raySubFig, fullfile(paths.save, 'Rayleigh100vs1000subsampled'), '-dsvg', '-r300');
