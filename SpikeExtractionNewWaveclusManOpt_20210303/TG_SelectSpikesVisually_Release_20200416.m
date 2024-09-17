%==========================================================================
% This script creates a file with the decision whether a given cluster
% shall be used for the analyses based on a visual inspection of its
% wave-clus output.
%
% Lukas Kunz, 2023 and Tim Guth, 2023
%==========================================================================

clear; clc; close all;
opengl hardware;

% paths
spike_path  = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303\';

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

% settings
myDate                      = '20210315_020120'; % current date and time
saveName                    = ['TG_', myDate];
% myDate                      = '20200416_113912';
% saveName                    = ['LK_', myDate]; % save name

% information about clusters
clusterInfo                 = [];

%% loop through subjects
for iSub = 1:length(subjects)
    
    % report
    fprintf('\nWorking on subject "%s" ...\n', subjects{iSub});
    cd(spike_path);
    
    % available sessions
    sess    = dir(strcat(spike_path, subjects{iSub}, filesep, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sess, 1)
        
        % get wires of this patient
        wires   = dir(fullfile(sess(iSess).folder, sess(iSess).name, 'chan*'));
        tmp     = split(transpose({wires.name}), 'n');
        [B, I]  = sort(cellfun(@str2num, tmp(:, 2)));
        wires   = wires(I);
        
        %% loop through wires
        for iWire = 1:size(wires, 1)
            
            % prellocate
            Cluster4Analysis    = [];
            
            % load spike data
            try
                % load wave-clus output
                t   = load(fullfile(wires(iWire).folder, wires(iWire).name, 'times_datacut.mat'));
            catch
                % save empty structure if no waveclus output available
                save(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', saveName), 'Cluster4Analysis'); % save empty structure
                warning('Wave-clus did not extract any spikes.');
                continue;
            end
            
            % load previous decision
            if exist(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', saveName, '.mat'), 'file') > 0
                prevDec = load(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', saveName, '.mat'));
            else
                clear prevDec; % no previous decision
            end
            
            %% loop through clusters
            for iClus = min(t.cluster_class(:, 1)):max(t.cluster_class(:, 1))
                
                % report
                fprintf('\n\n\nSubject: %s. Session: %s. Wire: %s. Cluster: %d.\n', ...
                    subjects{iSub}, sess(iSess).name, wires(iWire).name, iClus);
                
                % general information about this cluster
                tmpC4A          = [];
                tmpC4A.Cluster  = iClus; % cluster "name"
                tmpC4A.Date     = char(datetime("today")); % current date for bookkeeping
                
                % if you include all analyzable units
                if strcmp(saveName, 'All')
                    
                    % store decision
                    if iClus == 0
                        tmpC4A.Decision = 'no'; % do not include "rest" cluster
                    else
                        tmpC4A.Decision = 'yes';
                    end
                else
                    
                                        % perform visual control
                    if iClus == 0
                        tmpC4A.Decision = 'no'; % do not include "rest" cluster
                    else
                        
                        % decide whether to use this cluster or not
                        if exist('prevDec', 'var')
                            
                            % show previous decision
                            tmpPrevDec  = strcmp(prevDec.Cluster4Analysis(cell2mat({prevDec.Cluster4Analysis.Cluster}) == iClus).Decision, 'yes');
                            dec         = input(['Make a decision: shall this cluster be used for analyses (previous decision was ---', ...
                                num2str(tmpPrevDec), '---)? (1 = yes): ']);
                        else
                            
                            % show wave-clus image
                            f = figure('units', 'normalized', 'position', [0, 0, 1, 1]);
                            if iClus <= 3
                                Im  = imread(fullfile(wires(iWire).folder, wires(iWire).name, 'fig2print_datacut.png'));
                            elseif iClus <= 8
                                Im  = imread(fullfile(wires(iWire).folder, wires(iWire).name, 'fig2print_datacuta.png'));
                            else
                                Im  = imread(fullfile(wires(iWire).folder, wires(iWire).name, 'fig2print_datacutb.png'));
                            end
                            AxesH = axes('Units', 'pixels', 'position', [1, 1, 1536, 864], 'Visible', 'off'); % adjust according to screen size
                            image(Im, 'Parent', AxesH);
                            
                            % make a visual decision
                            dec = input('Make a decision: shall this cluster be used for analyses (previous decision was --- nothing ---)? (1 = yes): ');
                        end
                        if dec == 1
                            tmpC4A.Decision = 'yes';
                        else
                            tmpC4A.Decision = 'no';
                        end
                        
                        % close figures
                        close all;
                    end
                end
                
                % collect across all clusters of this wire
                Cluster4Analysis    = cat(1, Cluster4Analysis, tmpC4A);
                
                % collect info about all clusters
                if strcmp(tmpC4A.Decision, 'yes')
                    dec = 1; % recode decision information
                else
                    dec = 0;
                end
                clusterInfo = cat(1, clusterInfo, [iSub, iSess, iWire, iClus, dec]);
            end
            
            % save decision
            save(strcat(wires(iWire).folder, filesep, wires(iWire).name, filesep, 'Cluster4Analysis_', saveName), ...
                'Cluster4Analysis');
        end
    end
end

%% save summary output

% save
save(strcat(spike_path, saveName, '_clusterInfo'), 'clusterInfo');

% report
fprintf('Number of units offered by wave-clus: %d.\n', size(clusterInfo, 1));
fprintf('Number of units accepted for analysis: %d.\n', sum(clusterInfo(:, 5)));
