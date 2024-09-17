%==========================================================================
% This script performs spike-detection and -clustering using Wave-clus.
%
% Tim Guth, 2023
%==========================================================================

% start
clear; clc; close all;
addpath(genpath('D:\External\Toolboxes\Wave_clus')); % use new wave-clus (Chaure et al., 2018)

% paths
data_path   = 'D:\TreasureHunt\Micro_20210111\';
spike_path  = 'D:\TreasureHunt\SpikeExtractionNewWaveclusManOpt_20210303\';
pic_path    = strcat(spike_path, 'AllPics\');
if ~exist(pic_path, 'dir')
    mkdir(pic_path);
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

%% loop through subjects
parfor iSub = 1:size(subjects, 1)
    
    rng(1); % for reproducibility
    fprintf('\nWorking on subject %s ...\n', subjects{iSub});
    cd(spike_path);
    
    % available sessions
    sess    = dir(strcat(data_path, subjects{iSub}, filesep, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sess, 1)
        
        fprintf('\nWorking on subject: "%s", session: "%s" ...\n', subjects{iSub}, sess(iSess).name);
        
        %% list wires
        trigs   = dir(strcat(sess(iSess).folder, filesep, sess(iSess).name, filesep, 'ainp*'));
        wires   = dir(strcat(sess(iSess).folder, filesep, sess(iSess).name, filesep, 'chan*'));
        
        % sort files and combine
        tmp     = split(transpose({wires.name}), 'n');
        [~, I]  = sort(cellfun(@str2num, tmp(:, 2)));
        files   = [trigs; wires(I)];
        
        % report
        fprintf('\n------------ Performing wave_clus\n');
        
        %% loop through files and extract action potentials
        for iFile = 1:size(files, 1)
            
            % report
            fprintf('\tWorking on file "%s" ...\n', fullfile(files(iFile).folder, files(iFile).name));
            
            % storage destination
            spike_path_spec = fullfile(spike_path, subjects{iSub}, filesep, sess(iSess).name, filesep, files(iFile).name, filesep);
            
            % copy original data to storage destination          
            mkdir(spike_path_spec);
            status = copyfile(fullfile(files(iFile).folder, files(iFile).name, filesep, 'datacut.mat'), ...
                spike_path_spec);
            
            %% wave clus
            % needs input file containing "data" and "sr"
            
            if ~isempty(regexp(files(iFile).name, 'chan', 'once')) % if it is a channel-file
                
                % get file for wave-clus
                cd(spike_path_spec);
                file2use4wc     = {strcat(spike_path_spec, 'datacut.mat')};
                
                % run wave_clus to extract spikes and do clustering
                try
                    % spike extraction
                    myPar                   = [];
                    myPar.detection         = 'neg'; % determine detection type
                    myPar.randomseed        = 1; % define randomseed
                    Get_spikes(file2use4wc, 'par', myPar);
                    
                    % spike clustering
                    file2use4wc             = {strcat(spike_path_spec, 'datacut_spikes.mat')};
                    myPar                   = [];
                    myPar.randomseed        = 1; % define randomseed
                    myPar.template_sdnum    = 1.5; % for each unsorted spike, look up the closest template - but only if it is within its SD times 'par.template_sdnum'; (default 3)
                    myPar.min_clus          = 60; % minimum size of a cluster; (default 20)
                    myPar.max_clus        	= 10; % maximum number of clusters allowed (default 13)
                    myPar.mintemp           = 0.05; % minimum temperature
                    Do_clustering(file2use4wc, 'par', myPar);
                    
                    % copy output picture(s) to overview folder
                    pics = dir(fullfile(spike_path_spec, '*.png'));
                    for iPic = 1:length(pics)
                        copyfile(fullfile(pics(iPic).folder, pics(iPic).name), ...
                            strcat(pic_path, subjects{iSub}, '_', sess(iSess).name, '_', files(iFile).name, '_', pics(iPic).name));
                    end
                catch
                    warning('Wave-clus could not be performed.');
                end
                
                %% delete original input-file to save space
                delete(strcat(spike_path_spec, 'datacut.mat'));
            end
        end
    end
end
