%==========================================================================
% This script examines the memory performance in Treasure Hunt
%
% Tim Guth, 2023
%==========================================================================

%% settings
clc; close all; clear;

% set random seed
rng(444);

% add functions
addpath('C:\Sciebo\GitCode\NeuroGuth\TreasureHunt\Functions');

% paths
paths.behlog         = 'D:\TreasureHunt\Beh_20210111'; % behavioral logfile
paths.SUData         = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303'; % single-unit data from wave_clus
paths.save           = 'D:\TreasureHunt\BehAnalysis_20200528';

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

% configurations for random locations in the Treasure Hunt arena to normalize drop error
RandLoccfg            = [];
RandLoccfg.maxR       = 50;
RandLoccfg.minR       = 0;
RandLoccfg.N          = 1000000;
RandLoccfg.centerX    = 370;
RandLoccfg.centerY    = 360;
randLocs              = TG_RandomPointsInCircle(RandLoccfg);

% preallocate
objRecallPerf         = cell(size(subjects, 1), 1); % all
objRecallPerfThree    = cell(size(subjects, 1), 1); % trials with three chests
objRecallPerfFour     = cell(size(subjects, 1), 1); % trials with four chests
locRecallPerf         = cell(size(subjects, 1), 1); % all
locRecallPerfTwo      = cell(size(subjects, 1), 1); % trials with two chests
locRecallPerfThree    = cell(size(subjects, 1), 1); % trials with three chest
succLocRecDuration    = cell(size(subjects, 1), 1);
failLocRecDuration    = cell(size(subjects, 1), 1);
distractorDuration    = cell(size(subjects, 1), 1);
sessionDuration       = cell(size(subjects, 1), 1);

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % get sessions
    sessions                = dir(fullfile(paths.SUData, subjects{iSub}, 'session*'));
    
    % preallocate performance for this subject
    subjObjRecallPerf       = cell(size(sessions, 1), 1);
    subjObjRecallPerfThree  = cell(size(sessions, 1), 1);
    subjObjRecallPerfFour   = cell(size(sessions, 1), 1);
    subjLocRecallPerf       = cell(size(sessions, 1), 1);
    subjLocRecallPerfTwo    = cell(size(sessions, 1), 1);
    subjLocRecallPerfThree  = cell(size(sessions, 1), 1);
    subjSuccLocRecDur       = cell(size(sessions, 1), 1);
    subjFailLocRecDur       = cell(size(sessions, 1), 1);
    subjDistrDur            = cell(size(sessions, 1), 1);
    subjSessDur             = cell(size(sessions, 1), 1);
    bGoodMemLoc             = cell(size(sessions, 1), 1);
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % original session index
        sessIdx     = split(sessions(iSess).name, '_');
        sessIdx     = str2double(sessIdx{2});
        
        % display session information
        fprintf('\n==================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        % load behavioral logfile
        trialInfo                           = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo                           = trialInfo.trialInfo;
        
        % memory performance during object recall
        subjObjRecallPerf{iSess, 1}         = strcmp(trialInfo.audioResponse, 'CORRECT');

        % differentiate between 3 and 4 chest object recalls
        numObjIdx                           = repelem(trialInfo.numObjRecallPerTrial, trialInfo.numObjRecallPerTrial);
        subjObjRecallPerfThree{iSess, 1}    = subjObjRecallPerf{iSess, 1}(numObjIdx == 3);
        subjObjRecallPerfFour{iSess, 1}     = subjObjRecallPerf{iSess, 1}(numObjIdx == 4);
        
        % memory performance during location recall
        thisLocRecallPerf                   = nan(size(trialInfo.CORRECT_TEST_POSITION, 1), 1);
        for iTrial = 1:size(trialInfo.CORRECT_TEST_POSITION, 1)
            dropError                       = pdist2(trialInfo.CHOSEN_TEST_POSITION(iTrial, [1, 3]), trialInfo.CORRECT_TEST_POSITION(iTrial, [1, 3]));
            potDropErrors                   = pdist2(trialInfo.CORRECT_TEST_POSITION(iTrial, [1, 3]), randLocs);
            thisLocRecallPerf(iTrial, 1)    = sum(dropError < potDropErrors) / numel(potDropErrors);
        end
        subjLocRecallPerf{iSess, 1}         = thisLocRecallPerf;
        
        % differentiate between 2 and 3 chest location recalls
        numLocIdx                           = repelem(trialInfo.numLocRecallPerTrial, trialInfo.numLocRecallPerTrial);
        subjLocRecallPerfTwo{iSess, 1}      = subjLocRecallPerf{iSess, 1}(numLocIdx == 2);
        subjLocRecallPerfThree{iSess, 1}    = subjLocRecallPerf{iSess, 1}(numLocIdx == 3);

        % location recall durations
        bGoodMemLoc                         = thisLocRecallPerf > median(thisLocRecallPerf, 'omitnan'); % logical index for good location recall performance
        subjSuccLocRecDur{iSess, 1}         = trialInfo.timeLocRecall(bGoodMemLoc, 1) - trialInfo.timeCue4LocRecall(bGoodMemLoc, 1);
        subjFailLocRecDur{iSess, 1}         = trialInfo.timeLocRecall(~bGoodMemLoc, 1) - trialInfo.timeCue4LocRecall(~bGoodMemLoc, 1);
        
        % distractor game durations
        subjDistrDur{iSess, 1}              = trialInfo.DISTRACTOR_GAME_ENDED - trialInfo.DISTRACTOR_GAME_STARTED;

        % session duration
        subjSessDur{iSess, 1}               = max(trialInfo.FEEDBACK_ENDED) - min(trialInfo.HOMEBASE_TRANSPORT_STARTED);
    end
    
    % collapse across subjects
    objRecallPerf{iSub, 1}      = subjObjRecallPerf;
    objRecallPerfThree{iSub, 1} = subjObjRecallPerfThree;
    objRecallPerfFour{iSub, 1}  = subjObjRecallPerfFour;
    locRecallPerf{iSub, 1}      = subjLocRecallPerf;
    locRecallPerfTwo{iSub, 1}   = subjLocRecallPerfTwo;
    locRecallPerfThree{iSub, 1} = subjLocRecallPerfThree;
    succLocRecDuration{iSub, 1} = subjSuccLocRecDur;
    failLocRecDuration{iSub, 1} = subjFailLocRecDur;
    distractorDuration{iSub, 1} = subjDistrDur;
    sessionDuration{iSub, 1}    = subjSessDur;
end

%% collect results - object recall

% session-wise object recall performance
allSessObjRecallPerf        = cat(1, objRecallPerf{:});
allSessObjRecallPerfThree   = cat(1, objRecallPerfThree{:});
allSessObjRecallPerfFour    = cat(1, objRecallPerfFour{:});

% exclusion index
exclSessIdx                 = cell2mat(cellfun(@(x) isempty(x), allSessObjRecallPerf, 'UniformOutput', false));

% excluded sessions
allSessObjRecallPerf        = allSessObjRecallPerf(~exclSessIdx);
allSessObjRecallPerfThree   = allSessObjRecallPerfThree(~exclSessIdx);
allSessObjRecallPerfFour    = allSessObjRecallPerfFour(~exclSessIdx);

% average per session
avgObjRecallPerf            = cell2mat(cellfun(@(x) mean(x), allSessObjRecallPerf, 'UniformOutput', false));
avgObjRecallPerfThree       = cell2mat(cellfun(@(x) mean(x), allSessObjRecallPerfThree, 'UniformOutput', false));
avgObjRecallPerfFour        = cell2mat(cellfun(@(x) mean(x), allSessObjRecallPerfFour, 'UniformOutput', false));

% t-test to compare three versus four chest object recall
[~, objLoadP, ~, objLoadT]  = ttest(avgObjRecallPerfThree, avgObjRecallPerfFour);

%% collect results - location recall

% session-wise location recalls
allSessLocRecallPerf        = cat(1, locRecallPerf{:});
allSessLocRecallPerfTwo     = cat(1, locRecallPerfTwo{:});
allSessLocRecallPerfThree   = cat(1, locRecallPerfThree{:});

% excluded sessions
allSessLocRecallPerf        = allSessLocRecallPerf(~exclSessIdx);
allSessLocRecallPerfTwo     = allSessLocRecallPerfTwo(~exclSessIdx);
allSessLocRecallPerfThree   = allSessLocRecallPerfThree(~exclSessIdx);

% average per session
avgLocRecallPerf            = cell2mat(cellfun(@(x) mean(x), allSessLocRecallPerf, 'UniformOutput', false));
avgLocRecallPerfTwo         = cell2mat(cellfun(@(x) mean(x), allSessLocRecallPerfTwo, 'UniformOutput', false));
avgLocRecallPerfThree       = cell2mat(cellfun(@(x) mean(x), allSessLocRecallPerfThree, 'UniformOutput', false));

% t-test to compare two versus three chest location recall
[~, locLoadP, ~, locLoadT]  = ttest(avgLocRecallPerfTwo, avgLocRecallPerfThree);

%% figure for performance during (location-cued) object recall

% histogram of subject-wise memory accuracies
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.7, 3.8, 4]);
hold on;
objH = histogram(avgObjRecallPerf, 0:0.05:1, ...
    'FaceColor', [0.7, 0.7, 0.7]);
hold off;
tmpAx = get(gca);
set(gca, ...
    'xlim', [-0.05, 1.05], ...
    'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Memory performance');
yl = ylabel('# sessions', ...
    'units', 'normalized', 'position', [-0.05, 0.5, 0]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(gcf, 'PaperPositionMode', 'auto');
print(f, fullfile(paths.save, 'AllSubj_ObjectRecall_MemoryPerformance'), '-dsvg', '-r450');

%% figure for performance during (object-cued) location recall

% trial-wise location recall performance
trialLocRecallPerf  = cell2mat(tg_unwrapNestedCell(locRecallPerf));

% histogram of subject-wise memory accuracies
f = figure('units', 'centimeters', 'position', [2, 2, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.7, 3.8, 4]);
hold on;
locH = histogram(trialLocRecallPerf, 0:0.05:1, ...
    'FaceColor', [0.7, 0.7, 0.7]);
plot([0.5, 0.5], [0, max(locH.Values)], ':', ...
    'Color', [1 0 0], 'LineWidth', 2);
hold off;
tmpAx = get(gca);
set(gca, ...
    'xlim', [-0.05, 1.05], ...
    'ylim', [0, max(tmpAx.YLim)], 'ytick', [0, max(tmpAx.YLim)], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'box', 'off');
xl = xlabel('Memory performance');
yl = ylabel('# location recalls', ...
    'units', 'normalized', 'position', [-0.05, 0.5, 0]);
set(gca, ...
    'ylim', [0, max(locH.Values)], 'ytick', [0, max(locH.Values)]);
set([gca, xl, yl], ...
    'fontunits', 'centimeters', 'fontsize', 0.5);
% save figure
set(gcf, 'PaperPositionMode', 'auto');
print(f, fullfile(paths.save, 'AllSubj_LocationRecall_MemoryPerformance'), '-dsvg');

%% memory performance over time - object recall

% memory performance first versus second half
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.55, 3.8, 4]);
hold on;
objFirstSecond = zeros(size(allSessObjRecallPerf, 1), 2);
for iSub = 1:size(allSessObjRecallPerf, 1)
    objNumTrials            = size(allSessObjRecallPerf{iSub, 1}, 1);
    objHalf                 = ceil(objNumTrials / 2);
    objFirst                = allSessObjRecallPerf{iSub, 1}(1:objHalf);
    objSecond               = allSessObjRecallPerf{iSub, 1}((objHalf + 1):end);
    objFirstSecond(iSub, :) = [mean(objFirst), mean(objSecond)];
    plot([1, 2], objFirstSecond(iSub, :), '-', ...
        'Color', [0.5, 0.5, 0.5]);
end
plot([1, 2], mean(objFirstSecond), '-', ...
    'Color', 'blue', 'LineWidth', 3); % average across sessions
hold off;
set(gca, ...
    'xlim', [0.9, 2.1], 'xtick', [1, 2], 'xticklabel', {'First half', 'Second half'}, ...
    'ylim', [0, 1], 'ytick', [0, 1], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
xlabel('Trial', ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
ylabel('Memory performance', ...
    'units', 'normalized', 'position', [-0.15, 0.5, 0], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);

% t-test first vs second half
[~, objP, ~, objT] = ttest(objFirstSecond(:, 1), objFirstSecond(:, 2));
title(strcat('{P = }', num2str(objP, 3)));

% save figure
set(gcf, 'PaperPositionMode', 'auto');
print(f, fullfile(paths.save, 'AllSubj_ObjectRecall_MemoryPerformanceOverTime'), '-dsvg');

%% memory performance over time - location recall

% memory performance first versus second half
f = figure('units', 'centimeters', 'position', [5, 5, 6, 6]);
axes('units', 'centimeters', 'position', [1.3, 1.55, 3.8, 4]);
hold on;
locFirstSecond = zeros(size(allSessLocRecallPerf, 1), 2);
for iSub = 1:size(allSessLocRecallPerf, 1)
    locNumTrials            = size(allSessLocRecallPerf{iSub, 1}, 1);
    locHalf                 = ceil(locNumTrials / 2);
    locFirst                = allSessLocRecallPerf{iSub, 1}(1:locHalf);
    locSecond               = allSessLocRecallPerf{iSub, 1}((locHalf + 1):end);
    locFirstSecond(iSub, :) = [mean(locFirst), mean(locSecond)];
    plot([1, 2], locFirstSecond(iSub, :), '-', ...
        'Color', [0.5, 0.5, 0.5]);
end
plot([1, 2], mean(locFirstSecond), '-', ...
    'Color', 'blue', 'LineWidth', 3); % average across sessions
hold off;
set(gca, ...
    'xlim', [0.9, 2.1], 'xtick', [1, 2], 'xticklabel', {'First half', 'Second half'}, ...
    'ylim', [0, 1], 'ytick', [0, 1], ...
    'tickdir', 'out', 'ticklength', [0.02, 0.02], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
xlabel('Trial', ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);
ylabel('Memory performance', ...
    'units', 'normalized', 'position', [-0.15, 0.5, 0], ...
    'FontUnits', 'centimeters', 'FontSize', 0.5);

% t-test first vs second half
[~, locP, ~, locT] = ttest(locFirstSecond(:, 1), locFirstSecond(:, 2));
title(strcat('{P = }', num2str(locP, 3)));

% save figure
set(gcf, 'PaperPositionMode', 'auto');
print(f, fullfile(paths.save, 'AllSubj_LocationRecall_MemoryPerformanceOverTime'), '-dsvg');

%% trial duration successful vs. unsuccessful location recalls

% successful recalls
allSuccLocRecDuration      = cat(1, succLocRecDuration{:});
allSuccLocRecDuration      = allSuccLocRecDuration(~exclSessIdx);
avgSuccLocRecDuration      = cell2mat(cellfun(@(x) mean(x), allSuccLocRecDuration, 'UniformOutput', false));
grandAvgSuccLocRecDuration = mean(avgSuccLocRecDuration);

% unsuccessful recalls
allFailLocRecDuration      = cat(1, failLocRecDuration{:});
allFailLocRecDuration      = allFailLocRecDuration(~exclSessIdx);
avgFailLocRecDuration      = cell2mat(cellfun(@(x) mean(x), allFailLocRecDuration, 'UniformOutput', false));
grandAvgFailLocRecDuration = mean(avgFailLocRecDuration);

% t-test comparing the durations of successful and unsuccessful location recall
[~, durP, ~, durT]         = ttest(avgSuccLocRecDuration, avgFailLocRecDuration);

%% duration of distractor game
allDistractorDuration      = cat(1, distractorDuration{:});
meanDistractorDuration     = mean(cell2mat(allDistractorDuration));

%% session duration
alllSessDuration           = cat(1, sessionDuration{:});

% mean
meanSessionDuration        = mean(cell2mat(alllSessDuration));
meanSessionDurationMinutes = meanSessionDuration / 60;

% sem
semSessionDuration         = std(cell2mat(alllSessDuration)) / sqrt(length(cell2mat(alllSessDuration)));
semSessionMinutes          = semSessionDuration / 60;

%% correlation between session-wise object-recall performance and location-recall performance

% average recall performance per session
sessLocRecallPerf          = cellfun(@mean, tg_unwrapNestedCell(locRecallPerf));
sessLocRecallPerf          = sessLocRecallPerf(~isnan(sessLocRecallPerf)); % excluded sessions

% grand average object recall performance
grandAvgObjRecallPerf      = mean(avgObjRecallPerf);

% report
fprintf('\nNumber of subjects: %d.\n', size(objRecallPerf, 1));
fprintf('Number of sessions: %d.\n', size(avgObjRecallPerf, 1));
fprintf('Average object-recall performance: %.3f.\n', mean(avgObjRecallPerf));
fprintf('Percentage of location-recall trials above chance level: %.3f%%.\n', 100 * sum(trialLocRecallPerf > 0.5) / size(trialLocRecallPerf, 1));

% correlation
[rho, corrPval] = corr(avgObjRecallPerf, sessLocRecallPerf, 'type', 'spearman');
fprintf('Correlation between object-recall performance and location-recall performance: rho = %.3f, p = %.3f.\n', rho, corrPval);
