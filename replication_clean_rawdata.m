% Replication of clean raw data
% Requisites: EEGlab with plugins bva-io, bids-matlab-tools and firfilt
% clean_rawdata
% Tested on Matlab 2020a, Windows 10
% Cristina Gil Avila, 14.03.2022, TUM

clear; close all;
% EEGLAB start and clean_rawdata
run(fullfile('..','..','eeglab','eeglab.m'));
addpath(fullfile('..','..','clean_rawdata'));

% Datapaths
inputfolder = fullfile('..','..','eeg_datasets','rawBIDS');
outputfolder = fullfile('..','..','eeg_datasets','preprocessed_crd');
studyname = 'CGX_MS';
task = 'closed';

% Hard-coded parameters for simplicity
nReps = 10; % Number of repetitions
nChans = 29; % Number of channels
nSubj = 6; % Number of recordings
% Initialize bad channels
badchans = nan(nChans,nReps,nSubj);

%% Parameters to try:
% Try a higher transition band
highpass = [0.25 0.75]; % default
% highpass = [1 1.5];

% Number of samples for RANSAC
numsamples = 50; % default
% numsamples = 100;
% numsamples = 500;
% numsamples = 1000;

% Correlation parameter (channelCriterion)
channelcriterion = 0.7;
% channelcriterion = 0.8; % default
% channelcriterion = 0.9;

%% Loop over repetitions
for iRep =1:nReps
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab; % If I comment this line results are replicable
    % call BIDS importer
    [STUDY, ALLEEG] = pop_importbids(inputfolder,'bidsevent','on','bidschanloc','on',...
        'studyName',studyname,'bidstask',task,'outputdir',outputfolder);
    CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = 1:length(EEG);
    
    % Look up channel positions
    if isempty(find(~cellfun('isempty', {ALLEEG(1).chanlocs.X }), 1))
        EEG = pop_chanedit(EEG, 'lookup','elecpos_standard_1005.elc');
    end
    
    % Remove bad channels with clean_rawdata (Note: the field NumSamples is custom-added)
    EEG = pop_clean_rawdata(EEG,'FlatlineCriterion',5,'ChannelCriterion',channelcriterion,...
        'LineNoiseCriterion',4,'NumSamples',numsamples,'Highpass',highpass,...
        'BurstCriterion','off','WindowCriterion','off','BurstRejection','off',...
        'Distance','Euclidian','WindowCriterionTolerances','off');
    
    % Reshape bad channels for plotting
    etc = cellfun(@(x) x,{EEG.etc},'UniformOutput',0);
    ids = cellfun(@(x) isfield(x,'clean_channel_mask'),etc);
    if ~all(ids)
        for ix = find(~ids)
            etc{ix}.clean_channel_mask = true(nChans,1); % In case CRD did not detect any bad channels
        end
    end
    badchans(:,iRep,:) = ~cell2mat(cellfun(@(x) x.clean_channel_mask,etc,'UniformOutput',0));
    
end

% plot
avgchan = squeeze(mean(badchans,2));
stdchan = squeeze(std(badchans,0,2));
f = figure('Position',[1000, 62, 1019, 1276]);
for iSubj=1:nSubj
    subplot(6,1,iSubj)
    er = errorbar(avgchan(:,iSubj),stdchan(:,iSubj),'Marker','square','MarkerSize',8,...
        'MarkerFaceColor','#0072BD','LineStyle','none');
    title(['Subject ' num2str(iSubj)]);
    ylim([-0.5 1.5])
    yticks([0 1]);
    yticklabels({'Good','Bad'});
    xticks(1:nChans);
    xticklabels({ALLEEG(1).chanlocs.labels});
end
sgtitle(['ChannelCriterion = ' num2str(channelcriterion)]);
saveas(f,['replication_CRD_' num2str(channelcriterion) '_nSamples' num2str(numsamples) '.png']);
%     saveas(f,['replication_CRD_' num2str(channelcriterion(iC)) '_hp1Hz.png']);
close(f);
