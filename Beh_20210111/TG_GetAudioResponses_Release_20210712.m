%==========================================================================
% This script helps to classify the recorded audio responses during object
% recall as correct or incorrect and creates trialInfo.audioResponse.
%
% Tim Guth, 2021
%==========================================================================

%% settings
clc; close all; clear;

%% paths
paths.audio                = 'D:\TreasureHunt\Audio_20210111';
paths.behlog               = 'D:\TreasureHunt\Beh_20210111';

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

%% loop through subjects
for iSub = 1:length(subjects)
    
    % get sessions
    sessions = dir(fullfile(paths.audio, subjects{iSub}, filesep, 'session*'));
    
    %% loop through sessions
    for iSess = 1:size(sessions, 1)
        
        fprintf('\n===================== Subject: %s. Session: %s.\n', subjects{iSub}, sessions(iSess).name);
        
        % load behavioral logfile
        trialInfo               = load(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'));
        trialInfo               = trialInfo.trialInfo;
        
        if isfield(trialInfo, 'audioResponse')
            continue;
        end
        
        % preallocate audioResponse
        audioResponse           = cell(size(trialInfo.idxRecording, 1), 1);
        
        %% loop through object recalls
        for iRec = 1:size(trialInfo.idxRecording, 1)
            
            % read in audio file
            [audioFile, Fs] = audioread(fullfile(paths.audio, subjects{iSub}, sessions(iSess).name, strcat(trialInfo.idxRecording{iRec, 1}, '.wav')));
            sound(audioFile, Fs);
            
            % read text file with correct response and display it
            fid             = fopen(fullfile(paths.audio, subjects{iSub}, sessions(iSess).name, strcat(trialInfo.idxRecording{iRec, 1}, '.txt')));
            correctResponse = fread(fid, inf, 'uint8=>char')';
            fprintf(strcat('\n', trialInfo.idxRecording{iRec, 1},' Correct item: \n\n%s \n\n'), correctResponse);
            
            % ask user if correct or not
            bInput = input('Press "1" if correct and "0" if incorrect: ');
            if bInput == 1
                audioResponse{iRec, 1} = 'CORRECT';
            elseif bInput == 0
                audioResponse{iRec, 1} = 'INCORRECT';
            end
        end
        
        % add audioResponse to trialInfo
        trialInfo.audioResponse = audioResponse;
        
        % define position of audioResponse within structure
        oldFieldnames           = fieldnames(trialInfo);
        idxCortana              = find(contains(oldFieldnames, 'CORTANA_RESPONSE'));
        newFieldnames           = [oldFieldnames(1:idxCortana); oldFieldnames(end); oldFieldnames(idxCortana + 1:end - 1)];
        trialInfo               = orderfields(trialInfo, newFieldnames);
        
        % save new trailInfo-file
        mkdir(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name));
        save(fullfile(paths.behlog, subjects{iSub}, sessions(iSess).name, 'trialInfo.mat'), 'trialInfo');
        
    end
end
