%==========================================================================
% This script gets the behavioral data of TH into an adequate format.
% It also creates a new timeline with resampled behavioral data.
%--------------------------------------------------------------------------
% Output:
% - allBehData: behavioral logfile as cell
% - trialInfo:  time points and additional information for salient events
% - allBehInfo: for each logfile line, important information
% - behInfo:    for each time point, important information
% - behInfo10Hz: for each time point at 10 Hz, important information
%
% Tim Guth, 2021
%==========================================================================

% start
clc; close all; clear;
addpath(genpath('D:\External\Functions\'));

% paths
behdata_path           = 'D:\TreasureHunt\Beh_20210111\';
diagnostics_path       = strcat(behdata_path, 'Diagnostics_20200303\'); % for sanity checks
if ~exist(diagnostics_path, 'dir')
    mkdir(diagnostics_path);
end

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

%==========================================================================
%--------------------------------------------------------------------------
%% convert logfile into cell (allBehInfo) and get salient info (trialInfo)

%% loop through subjects
for iSub = 1:size(subjects, 1)
    
    % available sessions
    sessions = dir(strcat(behdata_path, subjects{iSub}, filesep, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        %% skip, if output data already exist
        
        % check whether "trialInfo" exists already
        if exist(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'), 'file') > 0
            fprintf('Subject: %s, session: %s. "allBehData" and "trialInfo" already exist, thus skipping ...\n', subjects{iSub}, sessions(iSess).name);
            continue;
        end
        
        %% read in data
        
        % report progress
        fprintf('\nSubject: %s, session: %s. Converting logfile into matlab cell ...\n', subjects{iSub}, sessions(iSess).name);
        tic;
        
        % open file
        file        = dir(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'data.txt'));
        numLines    = countlines(fullfile(file.folder, file.name)); % get number of lines of this file for preallocation
        tmpBehData  = cell(numLines, 1);
        
        % read first line as a start
        fid         = fopen(fullfile(file.folder, file.name));
        lineText    = fgetl(fid);
        idxCount 	= 1;
        
        % loop through remaining lines
        while ischar(lineText)
            
            % concatenate data
            tmpBehData{idxCount, 1} = split(lineText, char(9))';
            
            % read next line and count read lines
            lineText    = fgetl(fid);
            idxCount    = idxCount + 1;
        end
        
        % close file
        fclose(fid);
        
        % transform into better format (stable number of columns)
        numColumns      = cellfun(@(x) size(x, 2), tmpBehData);
        maxNumColumns   = max(numColumns);
        allBehData      = cell(numLines, maxNumColumns);
        for iT = 1:size(allBehData, 1)
            allBehData(iT, 1:numColumns(iT, 1)) = tmpBehData{iT, 1};
        end
        
        % save data as cell
        mkdir(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name));
        save(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'allBehData'), 'allBehData', '-v7.3');
        toc;
        
        %% extract important trial information
        
        % remove incomplete last trial for TH18, session_0
        if strcmp(subjects(iSub), 'TH18') && strcmp(sessions(iSess).name, 'session_0')
            trialEnds                             = find(contains(allBehData(:,4), 'FEEDBACK_ENDED'));
            allBehData(trialEnds(end) + 1:end, :) = [];
        end
        
        % structure to store trial information
        tI  = [];
        
        % trial index
        idx             = strcmp(allBehData(:, 3), 'Trial Info');
        numTrials       = sum(idx);
        tI.trialIdx     = transpose(1:numTrials);
        
        % HOMEBASE_TRANSPORT_STARTED
        idx                             = strcmp(allBehData(:, 4), 'HOMEBASE_TRANSPORT_STARTED') & cellfun(@isempty, allBehData(:, 5));
        tI.HOMEBASE_TRANSPORT_STARTED   = str2double(allBehData(idx, 1)) ./ 1000;
        
        % HOMEBASE_TRANSPORT_ENDED
        idx                         = strcmp(allBehData(:, 4), 'HOMEBASE_TRANSPORT_ENDED');
        tI.HOMEBASE_TRANSPORT_ENDED = str2double(allBehData(idx, 1)) ./ 1000;
        
        % TRIAL_NAVIGATION_STARTED
        idx                         = strcmp(allBehData(:, 4), 'TRIAL_NAVIGATION_STARTED');
        tI.TRIAL_NAVIGATION_STARTED = str2double(allBehData(idx, 1)) ./ 1000;
        
        % PLAYER_CHEST_ROTATION_STARTED
        idx                                 = strcmp(allBehData(:, 4), 'PLAYER_CHEST_ROTATION_STARTED');
        tmpTime                             = str2double(allBehData(idx, 1)) ./ 1000;
        tI.numChestsPerTrial                = transpose(histcounts(tmpTime, [tI.TRIAL_NAVIGATION_STARTED; Inf]));
        tI.chestIdx                         = transpose(1:size(tmpTime, 1));
        tI.PLAYER_CHEST_ROTATION_STARTED    = tmpTime;
        
        % PLAYER_CHEST_ROTATION_ENDED
        idx                             = strcmp(allBehData(:, 4), 'PLAYER_CHEST_ROTATION_ENDED');
        tI.PLAYER_CHEST_ROTATION_ENDED  = str2double(allBehData(idx, 1)) ./ 1000;
        
        % chest destroyed
        idx                     = contains(allBehData(:, 3), 'TreasureChest') & strcmp(allBehData(:, 4), 'DESTROYED');
        tI.chestDestroyedTime   = str2double(allBehData(idx, 1)) ./ 1000;
        
        % chest location
        idx             = ~cellfun(@isempty, regexp(allBehData(:, 3), 'TreasureChest')) & strcmp(allBehData(:, 4), 'POSITION');
        tI.chestLoc     = unique(str2double(allBehData(idx, 5:7)), 'rows', 'stable');
        
        % TREASURE_LABEL
        idx                     = strcmp(allBehData(:, 4), 'TREASURE_LABEL');
        tI.TREASURE_LABEL_GER   = allBehData(idx, 5);
        
        % English TREASURE_LABEL
        idx                     = strcmp(allBehData(:, 5), 'SpecialObject') & cellfun(@isempty, regexp(allBehData(:, 3), 'UICopy'));
        tI.TREASURE_LABEL_EN    = allBehData(idx, 3);
        
        % TRIAL_NAVIGATION_ENDED
        idx                         = strcmp(allBehData(:, 4), 'TRIAL_NAVIGATION_ENDED');
        tI.TRIAL_NAVIGATION_ENDED   = str2double(allBehData(idx, 1)) ./ 1000;
        
        % TOWER_TRANSPORT_STARTED
        idx                         = strcmp(allBehData(:, 4), 'TOWER_TRANSPORT_STARTED');
        tI.TOWER_TRANSPORT_STARTED  = str2double(allBehData(idx, 1)) ./ 1000;
        
        % TOWER_TRANSPORT_ENDED
        idx                         = strcmp(allBehData(:, 4), 'TOWER_TRANSPORT_ENDED');
        tI.TOWER_TRANSPORT_ENDED    = str2double(allBehData(idx, 1)) ./ 1000;
        
        % DISTRACTOR_GAME_STARTED
        idx                         = strcmp(allBehData(:, 4), 'DISTRACTOR_GAME_STARTED');
        tI.DISTRACTOR_GAME_STARTED  = str2double(allBehData(idx, 1)) ./ 1000;
        
        % distractor game: correctness of the response
        idx                     = ~cellfun(@isempty, regexp(allBehData(:, 3), 'BoxAudio')) & strcmp(allBehData(:, 4), 'AUDIO_PLAYING');
        tI.correctnessDisGame   = allBehData(idx, 3);
        
        % DISTRACTOR_GAME_ENDED
        idx                         = strcmp(allBehData(:, 4), 'DISTRACTOR_GAME_ENDED');
        tI.DISTRACTOR_GAME_ENDED    = str2double(allBehData(idx, 1)) ./ 1000;
        
        % recall type
        idx             = strcmp(allBehData(:, 3), 'TRIAL TYPE');
        tI.recallType   = allBehData(idx, 4); % 'LOCATION' ==> location has to be recalled; 'OBJECT' ==> object has to be recalled
        
        % RECALL_PHASE_STARTED
        idx                         = strcmp(allBehData(:, 4), 'RECALL_PHASE_STARTED');
        tI.RECALL_PHASE_STARTED     = str2double(allBehData(idx, 1)) ./ 1000;
        
        % recall of location: time of cueing object
        idx                         = ~cellfun(@isempty, regexp(allBehData(:, 5), 'Wo war'));
        tmpTime                     = str2double(allBehData(idx, 1)) ./ 1000;
        tI.numLocRecallPerTrial     = transpose(histcounts(tmpTime, [tI.TRIAL_NAVIGATION_STARTED; inf]));
        tI.timeCue4LocRecall        = tmpTime;
        
        % recall of location: time
        idx                 = strcmp(allBehData(:, 4), 'CHOSEN_TEST_POSITION');
        tI.timeLocRecall    = str2double(allBehData(idx, 1)) ./ 1000;
        
        % recall of location: CHOSEN_TEST_POSITION
        tI.CHOSEN_TEST_POSITION     = str2double(allBehData(idx, 5:7));
        
        % recall of location: CORRECT_TEST_POSITION
        idx                         = ~cellfun(@isempty, regexp(allBehData(:, 4), 'CORRECT_TEST_POSITION', 'once'));
        tI.CORRECT_TEST_POSITION    = str2double(allBehData(idx, 5:7));
        
        % recall of location: object that is used for cueing
        tI.cueingObject     = allBehData(idx, 8);
        
        % recall of object: OBJECT_RECALL_CHOICE_STARTED
        idx                                 = strcmp(allBehData(:, 4), 'OBJECT_RECALL_CHOICE_STARTED');
        tmpTime                             = str2double(allBehData(idx, 1)) ./ 1000;
        tI.numObjRecallPerTrial             = transpose(histcounts(tmpTime, [tI.TRIAL_NAVIGATION_STARTED; Inf]));
        tI.OBJECT_RECALL_CHOICE_STARTED     = tmpTime;
        
        % RECORDING_STARTED and index of recording
        idx                     = strcmp(allBehData(:, 3), 'RECORDING_STARTED');
        tI.RECORDING_STARTED    = str2double(allBehData(idx, 1)) ./ 1000;
        tI.idxRecording         = allBehData(idx, 4);
        
        % RECORDING_ENDED
        idx                 = strcmp(allBehData(:, 3), 'RECORDING_ENDED');
        tI.RECORDING_ENDED  = str2double(allBehData(idx, 1)) ./ 1000;
        
        % recall of object: OBJECT_RECALL_CHOICE_ENDED
        idx                             = strcmp(allBehData(:, 4), 'OBJECT_RECALL_CHOICE_ENDED');
        tI.OBJECT_RECALL_CHOICE_ENDED   = str2double(allBehData(idx, 1)) ./ 1000;
        
        % recall of object: location that is used for cueing
        tmpLines    = str2double(allBehData(strcmp(allBehData(:, 4), 'OBJECT_RECALL_CHOICE_STARTED'), 2));
        idx       	= ismember(str2double(allBehData(:, 2)), tmpLines) & ...
            strcmp(allBehData(:, 3), 'ObjectSelectorVisuals') & strcmp(allBehData(:, 4), 'POSITION');
        tI.cueingLoc    = unique(str2double(allBehData(idx, 5:7)), 'rows', 'stable');
        
        % recall of object: correct or incorrect as defined by Cortana
        idx                 = strcmp(allBehData(:, 3), 'CORTANA_RESPONSE');
        tI.CORTANA_RESPONSE = allBehData(idx, 4);
        
        % TEMPORAL_RETRIEVAL_STARTED
        idx                             = strcmp(allBehData(:, 4), 'TEMPORAL_RETRIEVAL_STARTED');
        tI.TEMPORAL_RETRIEVAL_STARTED   = str2double(allBehData(idx, 1)) ./ 1000;
        
        % OPTION_A_vs_B
        idx                     = strcmp(allBehData(:, 4), 'OPTION_A');
        tmpTime                 = str2double(allBehData(idx, 1)) ./ 1000;
        tI.numTempRetPerTrial   = transpose(histcounts(tmpTime, [tI.TRIAL_NAVIGATION_STARTED; Inf]));
        tI.timeTempRet          = tmpTime;
        tI.OPTION_A_vs_B        = allBehData(idx, [5, 7]);
        
        % temporal retrieval: correctness of the response
        idx                     = strcmp(allBehData(:, 3), 'TEMPORAL_RETRIEVAL') & ~cellfun(@isempty, regexp(allBehData(:, 4), '_RESPONSE'));
        tI.correctnessTempRet   = allBehData(idx, 4);
        
        % TEMPORAL_RETRIEVAL_ENDED
        idx                         = strcmp(allBehData(:, 4), 'TEMPORAL_RETRIEVAL_ENDED');
        tI.TEMPORAL_RETRIEVAL_ENDED = str2double(allBehData(idx, 1)) ./ 1000;
        
        % RECALL_PHASE_ENDED
        idx                     = strcmp(allBehData(:, 4), 'RECALL_PHASE_ENDED');
        tI.RECALL_PHASE_ENDED   = str2double(allBehData(idx, 1)) ./ 1000;
        
        % FEEDBACK_STARTED
        idx                 = strcmp(allBehData(:, 4), 'FEEDBACK_STARTED');
        tI.FEEDBACK_STARTED = str2double(allBehData(idx, 1)) ./ 1000;
        
        % FEEDBACK_ENDED
        idx                 = strcmp(allBehData(:, 4), 'FEEDBACK_ENDED');
        tI.FEEDBACK_ENDED   = str2double(allBehData(idx, 1)) ./ 1000;
        
        % rename main variable
        trialInfo   = tI;
        
        %% sanity-check figures
        
        % produces several sanity-check-figues
        try
            cfg                     = [];
            cfg.diagnostics_path    = diagnostics_path;
            cfg.subject             = subjects{iSub, 1};
            cfg.session             = sessions(iSess).name;
            tg_trialInfoDiagnostics_20200511(cfg, trialInfo);
        catch
            fprintf('WARNING: Could not create all diagnostic images!\n');
        end
        
        %% save behavioral output
        save(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'trialInfo'), 'trialInfo');
        
    end
end

%==========================================================================
%--------------------------------------------------------------------------
%% extract timepoint-specific information from original logfile (original timeline)

%% loop through subjects
 for iSub = 1:size(subjects, 1)

    % available sessions
    sessions 	= dir(strcat(behdata_path, subjects{iSub}, filesep, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        %% skip, if output data already exist
        
        % check whether "behInfo" exists already
        if exist(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'behInfo.mat'), 'file') > 0
            fprintf('Subject: %s, session: %s. "allBehInfo" and "behInfo" already exist, thus skipping ...\n', subjects{iSub}, sessions(iSess).name);
            continue;
        end
        
        %% read in data
        
        % report progress
        fprintf('Subject: %s, session: %s. Creating timepoint specific behavioral information file ...\n', subjects{iSub}, sessions(iSess).name);
        
        % open original logfile
        file        = dir(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'data.txt'));
        numLines    = countlines(fullfile(file.folder, file.name)); % get number of lines of this file for preallocation
        fid         = fopen(fullfile(file.folder, file.name));
        
        % go through behavioral file line-by-line
        thisPlayerPosition  = NaN(1, 3); % at the beginning, you don't know the current player position
        thisPlayerRotation  = NaN(1, 3); % at the beginning, you don't know the current player rotation
        idxTrial            = 0; % trial index
        idxPart             = 0; % trial part index (i.e., current chest number)
        thisPhase           = 'none'; % at the beginning, you don't know the current phase
        thisTreas           = 'none'; % at the beginning, you don't know the current treasure
        recallType          = 'none';
        
        % preallocate "allBehInfo" variable, which contains the relevant
        % behavioral information for all lines of the logfile
        allBehInfo          = cell(numLines, 12);
        idxCount            = 1;
        lineText            = fgetl(fid); % read first line as a start
        while ischar(lineText)
            
            % preallocate behavioral information extracted from this line
            thisBehInfo     = cell(1, 12); % time, location (xyz), rotation (Euler xyz), trialidx, trialpart, trialphase, treasurename, recalltype
            
            % get information stored in this line
            lineSplits      = strsplit(lineText, '\t');
            
            % get time during this line
            thisTime        = str2double(lineSplits{1}) ./ 1000; % convert into seconds
            thisBehInfo{1}  = thisTime;
            
            % get player position and rotation during this line
            if strcmp(lineSplits{3}, 'Player') && strcmp(lineSplits{4}, 'POSITION') % update position if available
                thisPlayerPosition  = str2double(lineSplits(5:7)); % xyz
            elseif strcmp(lineSplits{3}, 'Player') && strcmp(lineSplits{4}, 'ROTATION') % update rotation if available
                thisPlayerRotation  = str2double(lineSplits(5:7)); % Euler xyz (= roll, yaw, pitch)
            end
            thisBehInfo(2:4)    = num2cell(thisPlayerPosition);
            thisBehInfo(5:7)    = num2cell(thisPlayerRotation);
            
            % get trial index
            if strcmp(lineSplits{3}, 'TRIAL TYPE')
                idxTrial  	= idxTrial + 1; % increase trial-index
                fprintf('\tTrial %d starts at line %d.\n', idxTrial, idxCount);
            end
            
            % get trial-phase during this line
            if strcmp(lineSplits{4}, 'INSTRUCTION_VIDEO_STARTED')
                thisPhase   = 'video'; % video before experiment
            elseif strcmp(lineSplits{4}, 'INSTRUCTION_VIDEO_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'HOMEBASE_TRANSPORT_STARTED')
                thisPhase   = 'homeTrans'; % allo2ego
            elseif strcmp(lineSplits{4}, 'HOMEBASE_TRANSPORT_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'TRIAL_NAVIGATION_STARTED')
                thisPhase   = 'navigation'; % navigation
                idxPart     = 1; % trial part index (first chest)
            elseif strcmp(lineSplits{4}, 'PLAYER_CHEST_ROTATION_STARTED')
                thisPhase   = 'chestRot'; % chest rotation
            elseif strcmp(lineSplits{4}, 'TREASURE_OPEN')
                thisPhase   = 'treasOpen';
            elseif contains(lineSplits{3}, 'TreasureChest') && strcmp(lineSplits{4}, 'DESTROYED')
                thisPhase   = 'navigation'; % navigation
                idxPart     = idxPart + 1; % towards next chest
                thisTreas   = 'none';
            elseif strcmp(lineSplits{4}, 'TRIAL_NAVIGATION_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'TOWER_TRANSPORT_STARTED')
                thisPhase   = 'towerTrans'; % ego2allo
                idxPart     = 0; % reset trial part index
            elseif strcmp(lineSplits{4}, 'TOWER_TRANSPORT_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'DISTRACTOR_GAME_STARTED')
                thisPhase   = 'distractor'; % distractor
            elseif strcmp(lineSplits{4}, 'DISTRACTOR_GAME_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'RANDOM_JITTER_STARTED')
                thisPhase   = 'jitter';
            elseif strcmp(lineSplits{4}, 'RANDOM_JITTER_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'LOCATION_RECALL_CHOICE_STARTED')
                thisPhase   = 'locRecall'; % retrieval of location
            elseif strcmp(lineSplits{4}, 'LOCATION_RECALL_CHOICE_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'OBJECT_RECALL_CHOICE_STARTED')
                thisPhase   = 'objRecall'; % retrieval of object
            elseif strcmp(lineSplits{4}, 'OBJECT_RECALL_CHOICE_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'TEMPORAL_RETRIEVAL_STARTED')
                thisPhase   = 'tempRet'; % temporal retrieval
            elseif strcmp(lineSplits{4}, 'TEMPORAL_RETRIEVAL_ENDED')
                thisPhase   = 'wait4Next';
            elseif strcmp(lineSplits{4}, 'FEEDBACK_STARTED')
                thisPhase   = 'feedback'; % feedback
            elseif strcmp(lineSplits{4}, 'FEEDBACK_ENDED')
                thisPhase   = 'wait4Next';
            end
            
            % get treasure label while chest is open
            if strcmp(lineSplits{4}, 'TREASURE_LABEL')
                thisTreas    = lineSplits{5};
            end
            
            % information about recall type
            if strcmp(lineSplits{3}, 'TRIAL TYPE')
                if strcmp(lineSplits{4}, 'OBJECT')
                    recallType  = 'objRecall'; % location-cued object recall
                elseif strcmp(lineSplits{4}, 'LOCATION')
                    recallType  = 'locRecall'; % object-cued location recall
                end
            end
            
            % fill in information
            thisBehInfo{8}      = idxTrial;
            thisBehInfo{9}      = idxPart;
            thisBehInfo{10}     = thisPhase;
            thisBehInfo{11}     = thisTreas;
            thisBehInfo{12}     = recallType;
            
            % collapse information
            allBehInfo(idxCount, :)    = thisBehInfo;
            
            % read next line and count read lines
            lineText    = fgetl(fid);
            idxCount   	= idxCount + 1;
        end
        
        % close logfile
        fclose(fid);
        
        % save entire behavioral output
        save(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'allBehInfo'), 'allBehInfo', '-v7.3');
        
        %% now condense information from identical timepoints
        
        % get start indices of time chunks
        startIdx    = find(diff([0; cell2mat(allBehInfo(:, 1))]) > 0);
        endIdx      = find(diff([cell2mat(allBehInfo(:, 1)); inf]) > 0);
        fprintf('There are %d start-indices and %d end-indices.\n', numel(startIdx), numel(endIdx));
        
        % preallocate
        behInfo     = cell(size(startIdx, 1), size(allBehInfo, 2));
        parfor (iTime = 1:size(startIdx, 1), 3) % M indicates maximum number of workers
            
            % report progress
            if mod(iTime, 100000) == 0
                fprintf('Converted %d time points (out of %d).\n', iTime, numel(startIdx))
            end
            
            % get most prevalent information during the chunk
            thisChunk   = allBehInfo(startIdx(iTime):endIdx(iTime), :);
            tmpBehInfo  = cell(1, size(allBehInfo, 2));
            for iColumn = 1:size(thisChunk, 2)
                
                % extract main data from this chunk
                if isa(thisChunk{1, iColumn}, 'char')
                    [s, ~, j]               = unique(thisChunk(:, iColumn));
                    tmpBehInfo{1, iColumn}  = s{mode(j)};
                else
                    tmpBehInfo{1, iColumn}  = mode(cell2mat(thisChunk(:, iColumn)));
                end
            end
            
            % collect information
            behInfo(iTime, :)   = tmpBehInfo;
        end
        
        % save condensed behavioral output
        save(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'behInfo'), 'behInfo');
        
        %% information about real-world start and end time
        
        % start time and end time
        starttime   = min(cell2mat(behInfo(:, 1)));
        realworld_starttime = datetime(starttime, ...
            'ConvertFrom', 'posixTime', ...
            'TimeZone', 'Europe/Zurich', ...
            'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
        endtime     = max(cell2mat(behInfo(:, 1)));
        realworld_endtime = datetime(endtime, ...
            'ConvertFrom', 'posixTime', ...
            'TimeZone', 'Europe/Zurich', ...
            'Format', 'dd-MMM-yyyy HH:mm:ss.SSS');
        fprintf('\nReal world start time = %s, real world end time = %s ...\n', realworld_starttime, realworld_endtime);
        
        % get experiment duration
        fprintf('Experiment duration = %.3f minutes ...\n\n', range(cell2mat(behInfo(:, 1))) / 60); % convert to minutes
        fprintf('Number of trials = %d ...\n\n', max(cell2mat(behInfo(:, 8)))); % number of trials
        
        %% plot player positions during navigation-phase over time
        % the xz-coordinates represent the path
        
        % plot path of participant
        f = figure('units', 'centimeters', 'position', [1, 1, 10, 10]);
        logIdx  = strcmp(behInfo(:, 10), 'navigation') | strcmp(behInfo(:, 10), 'chestRotation');
        plot(cell2mat(behInfo(logIdx, 2)), cell2mat(behInfo(logIdx, 4)), ...
            '.', 'Color', [0.5, 0.5, 0.5]);
        box on;
        axis equal;
        set(gca, 'xdir', 'reverse', 'ydir', 'reverse', ...
            'xlim', [310, 430], 'ylim', [300, 420]);
        tl = title({'Navigation path', strrep([subjects{iSub}, ', ', sessions(iSess).name], '_', '\_')});
        xl = xlabel('x (vu)');
        yl = ylabel('z (vu)');
        set([gca, xl, yl, tl], ...
            'fontunits', 'centimeters', 'fontsize', 0.5, 'fontweight', 'normal');
        % print figure
        set(gcf, 'PaperPositionMode', 'auto');
        print(f, '-dtiff', ...
            fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'navigationPath.tiff'), ...
            '-r450');
        
        % close figure
        close(f);
    end
end


%==========================================================================
%--------------------------------------------------------------------------
%% process behavioral information into new timeline

%% loop through subjects
for iSub = 1:length(subjects)
    
    % sessions performed by this subject
    sessions 	= dir(strcat(behdata_path, subjects{iSub}, filesep, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        % report progress
        fprintf('\nBehavioral processing of %s, %s ...\n', subjects{iSub}, sessions(iSess).name);
        
        % load "behInfo" that contains, for each timepoint:
        % time, xyz, Euler xyz, trialidx, trialpartidx, trialphase, objectname, retrievaltype
        tmp     = load(fullfile(sessions(iSess).folder, sessions(iSess).name, 'behInfo.mat'));
        behInfo = tmp.behInfo;
        fprintf('Original size of "behInfo": \t%d x %d.\n', size(behInfo));
        fprintf('The data segment lasts %.2f min.\n', range(cell2mat(behInfo(:, 1))) / 60);
        
        %% original timeline
        
        % store original behavioral information
        orig                = [];
        orig.behInfo        = behInfo;
        
        % time, durations, and distances
        orig.time           = cell2mat(orig.behInfo(:, 1));
        orig.durations      = [diff(orig.time); median(diff(orig.time))]; % add an average duration at the end
        orig.xyz            = cell2mat(orig.behInfo(:, 2:4)); % note: the subject moves in the xz-plane
        orig.distances      = [sqrt(diff(orig.xyz(:, 1)) .^ 2 + diff(orig.xyz(:, 3)) .^ 2); 0]; % add a zero-distance at the end (no moving)
        orig.speed          = orig.distances ./ orig.durations;
        
        % get yaws, to this end convert Euler y values into yaw values that
        % fit with the xz plane using this conversion:
        % yaw = deg2rad(wrapTo180(-1 .* EulerY + 90))
        tmpEulerY           = cell2mat(orig.behInfo(:, 6));
        tmpYaws             = wrapTo180(-1 .* tmpEulerY + 90);
        orig.yaws           = deg2rad(tmpYaws); % map onto [-pi, pi]
        
        % calculate navigation direction based on subsequent xy locations
        % cave: navigation direction is similar, but not identical, to yaws
        % and can only be determined for navigation periods
        orig.navDir         = [atan2(diff(orig.xyz(:, 3)), diff(orig.xyz(:, 1))); 0]; % atan2(y, x); add a 0° navigation direction at the end
        
        % estimate angular speed
        orig.angDistances   = abs([angdiff(orig.yaws); 0]); % add a zero-distance at the end (no turning)
        orig.angSpeed       = orig.angDistances ./ orig.durations;
        
        %% new timeline
        
        % define new timeline with 0.1 sec time resolution
        new             = [];
        new.timeRes     = 0.1; % e.g., 100 msec for 10 Hz resolution
        newFirstTime    = ceil(min(orig.time)); % [sec]
        newLastTime     = floor(max(orig.time)); % [sec]
        newTimeLine     = transpose(newFirstTime:new.timeRes:newLastTime);
        
        % give information
        fprintf('Creating a new timeline with a temporal resolution of %.1f Hz.\n', 1 / new.timeRes);
        fprintf('Original first time: %.3f, new first time: %.3f.\n', min(orig.time), newFirstTime);
        fprintf('Original last time:  %.3f, new last time:  %.3f.\n', max(orig.time), newLastTime);
        fprintf('Overall duration of the new timeline: %.1f sec.\n', newLastTime - newFirstTime);
        
        % preallocate output and add parfor progressbar
        allNewBehInfo   = cell(size(newTimeLine, 1), size(orig.behInfo, 2));
        ppm             = ParforProgressbar(size(allNewBehInfo, 1) - 1, 'progressBarUpdatePeriod', 2);
        parfor (iT = 1:size(allNewBehInfo, 1) - 1, 7) % M indicates maximum number of workers
            
            % preallocate new data line
            thisNewBehInfo  = cell(1, size(orig.behInfo, 2));
            
            % get section of original behInfo-file with data from new time
            % window
            bTWOI           = cell2mat(orig.behInfo(:, 1)) >= newTimeLine(iT, 1) & cell2mat(orig.behInfo(:, 1)) < newTimeLine(iT + 1, 1);
            thisBehinfo     = orig.behInfo(bTWOI, :);
            if isempty(thisBehinfo)
                allNewBehInfo(iT, :)    = num2cell(nan(1, size(thisNewBehInfo, 2)));
                continue;
            end
            
            %% behInfo
            
            % x, y, z
            thisNewBehInfo(1, 2:4)   	= num2cell(mode(cell2mat(thisBehinfo(:, 2:4)), 1)); % mode along first dimension
            
            % Euler x, y, z (angular variables in the range [0, 360] !)
            thisEulerY      = deg2rad(cell2mat(thisBehinfo(:, 6))); % convert into radians
            thisEulerY    	= thisEulerY(~isnan(thisEulerY)); % remove nans
            if ~isempty(thisEulerY)
                thisNewBehInfo(1, 5:7)     = num2cell([nan, wrapTo360(rad2deg(circ_mean(thisEulerY))), nan]); % convert back into degrees; don't fill in values for Euler x and z
            else
                thisNewBehInfo(1, 5:7)     = num2cell([nan, nan, nan]); % necessary, because circ_mean gives 0 when input is empty
            end
            
            % trialidx, trialpart
            thisNewBehInfo(1, 8:9)  = num2cell(mode(cell2mat(thisBehinfo(:, 8:9)), 1)); % along first dimension
            
            % trialphase
            [s, ~, j]            	= unique(thisBehinfo(:, 10));
            thisNewBehInfo{1, 10} 	= s{mode(j, 1)};
            
            % objectname
            [s, ~, j]           	= unique(thisBehinfo(:, 11));
            thisNewBehInfo{1, 11} 	= s{mode(j, 1)};
            
            % retrievaltype
            [s, ~, j]               = unique(thisBehinfo(:, 12));
            thisNewBehInfo{1, 12}   = s{mode(j, 1)};
            
            % add to overall variable
            allNewBehInfo(iT, :) 	= thisNewBehInfo;
            
            % update parfor progressbar
            ppm.increment();
        end
        
        % delete parfor progressbar
        delete(ppm);
        
        % fill last line with nans
        allNewBehInfo(end, :)   = num2cell(nan(1, size(allNewBehInfo, 2)));
        
        % substitute lines with nans with content of previous lines
        nanIdx                  = find(isnan(cell2mat(allNewBehInfo(:, 6))));
        for iNan = 1:size(nanIdx, 1)
            if (nanIdx(iNan, 1) - 1) > 0
                allNewBehInfo(nanIdx(iNan, 1), :)   = allNewBehInfo(nanIdx(iNan, 1) - 1, :);
            end
        end
        
        % insert into new "behInfo"-structure
        new.behInfo             = allNewBehInfo;
        
        %% compute further relevant information in new timeline
        
        % time, durations, and distances
        new.time            = cell2mat(new.behInfo(:, 1));
        new.durations       = [diff(new.time); median(diff(new.time))]; % add an average duration at the end
        new.xyz             = cell2mat(new.behInfo(:, 2:4)); % note: the subject moves in the xz-plane
        new.distances       = [sqrt(diff(new.xyz(:, 1)) .^ 2 + diff(new.xyz(:, 3)) .^ 2); 0]; % add a zero-distance at the end (no moving)
        new.speed           = new.distances ./ new.durations;
        
        % get yaws, to this end convert Euler y values into yaw values that
        % fit with the xz plane using this conversion:
        % yaw = deg2rad(wrapTo180(-1 .* EulerY + 90))
        tmpEulerY           = cell2mat(new.behInfo(:, 6));
        tmpYaws             = wrapTo180(-1 .* tmpEulerY + 90); % range [-180, 180]
        new.yaws            = deg2rad(tmpYaws); % map onto [-pi, pi]
        
        % calculate navigation direction based on subsequent xy locations
        % cave: navigation direction is similar, but not identical, to yaws
        % and can only be determined for navigation periods
        new.navDir          = [atan2(diff(new.xyz(:, 3)), diff(new.xyz(:, 1))); 0]; % atan2(y, x); add a 0° navigation direction at the end
        
        % estimate angular speed
        new.angDistances    = abs([angdiff(new.yaws); 0]); % add a zero-distance at the end (no turning)
        new.angSpeed        = new.angDistances ./ new.durations;
        
        %% sanity-check figure
        
        f = figure('units', 'centimeters', 'position', [1, 1, 40, 16]);
        
        % original navDirs and yaws, no masking
        axes('units', 'centimeters', 'position', [1.5, 9.5, 6, 6]);
        plot(orig.navDir, orig.yaws, 'k.');
        xl = xlabel('orig.navDir');
        yl = ylabel('orig.yaws');
        set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.35);
        
        % original navDirs and yaws, only navigation periods with speed
        % threshold
        axes('units', 'centimeters', 'position', [1.5, 1.5, 6, 6]);
        bMask   = strcmp(orig.behInfo(:, 10), 'navigation') & orig.speed > 0.01;
        plot(orig.navDir(bMask), orig.yaws(bMask), 'k.');
        xl = xlabel('orig.navDir(navi & speed>0.01)');
        yl = ylabel('orig.yaws(navi & speed>0.01)');
        set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.35);
        
        % new navDirs and yaws, no masking
        axes('units', 'centimeters', 'position', [9.5, 9.5, 6, 6]);
        plot(new.navDir, new.yaws, 'b.');
        xl = xlabel('new.navDir');
        yl = ylabel('new.yaws');
        set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.35);
        
        % new navDirs and yaws, only navigation periods with speed
        % threshold
        axes('units', 'centimeters', 'position', [9.5, 1.5, 6, 6]);
        bMask   = strcmp(new.behInfo(:, 10), 'navigation') & new.speed > 0.01;
        plot(new.navDir(bMask), new.yaws(bMask), 'b.');
        xl = xlabel('new.navDir(navi & speed>0.01)');
        yl = ylabel('new.yaws(navi & speed>0.01)');
        set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.35);
        
        % direct comparison between navDirs
        axes('units', 'centimeters', 'position', [17.5, 9.5, 22, 6]);
        hold on;
        plot(orig.time - orig.time(1), orig.navDir);
        plot(new.time - orig.time(1), new.navDir, '.');
        axis tight;
        xlabel('Time (sec)');
        ylabel('navDir');
        set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.35);
        
        % direct comparison between yaws
        axes('units', 'centimeters', 'position', [17.5, 1.5, 22, 6]);
        hold on;
        plot(orig.time - orig.time(1), orig.yaws);
        plot(new.time - orig.time(1), new.yaws, '.');
        axis tight;
        xlabel('Time (sec)');
        ylabel('yaws');
        set([gca, xl, yl], 'fontunits', 'centimeters', 'fontsize', 0.35);
        
        % save figure
        set(gcf, 'PaperPositionMode', 'auto');
        print(f, '-dtiff', ...
            fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, 'newTimeLine_controlNavDirsAndYaws.tiff'), ...
            '-r450');
        
        %% save output
        
        % combine original and new behavioral data
        behInfo         = [];
        behInfo.orig    = orig;
        behInfo.new     = new;
        save(fullfile(behdata_path, subjects{iSub}, sessions(iSess).name, ['behInfo', num2str(1 / new.timeRes), 'Hz']), 'behInfo');
        
        % close all open figures
        close all;
    end
end


%==========================================================================
%--------------------------------------------------------------------------
% GENERAL INFORMATION ABOUT THE TH LOGFILE
%--------------------------------------------------------------------------
%==========================================================================

%% Trial Event
%==========================================================================
% INSTRUCTION_VIDEO_STARTED (==> before experiment)
% INSTRUCTION_VIDEO_ENDED
% SHOWING_INSTRUCTIONS
% ----------------------------------------------------------------- n-times
% HOMEBASE_TRANSPORT_STARTED (==> homebase)
% HOMEBASE_TRANSPORT_ENDED
% SHOWING_INSTRUCTIONS
%
% TRIAL_NAVIGATION_STARTED (==> navigation)
%--- n-times
% PLAYER_CHEST_ROTATION_STARTED (==> encoding)
% PLAYER_CHEST_ROTATION_ENDED (==> navigation)
%---
% TRIAL_NAVIGATION_ENDED
%
% TOWER_TRANSPORT_STARTED (==> translation egocentric --> allocentric)
% TOWER_TRANSPORT_ENDED
%
% DISTRACTOR_GAME_STARTED (==> distractor)
% DISTRACTOR_GAME_ENDED
%
% RECALL_PHASE_STARTED (==> location/object retrieval)
%--- n-times
% LOCATION_RECALL_CHOICE_STARTED    | OBJECT_RECALL_CHOICE_STARTED
% LOCATION_RECALL_CHOICE_ENDED      | OBJECT_RECALL_CHOICE_ENDED
% SHOWING_INSTRUCTIONS
%---
% TEMPORAL_RETRIEVAL_STARTED (==> temporal retrieval)
% TEMPORAL_RETRIEVAL_ENDED
% RECALL_PHASE_ENDED
%
% FEEDBACK_STARTED (==> feedback)
% (SHOWING_INSTRUCTIONS)
% TEMPORAL_FEEDBACK_STARTED
% TEMPORAL_FEEDBACK_ENDED
% (SHOWING_INSTRUCTIONS)
% SCORESCREEN_STARTED
% SCORESCREEN_ENDED
% FEEDBACK_ENDED
%--------------------------------------------------------------- in-between
% TASK_PAUSED True
% TASK_PAUSED False
%==========================================================================
