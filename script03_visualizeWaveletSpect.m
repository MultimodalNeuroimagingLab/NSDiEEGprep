%% This is the preprocessing script for the clean NSD analysis.
% Requires each subject to have annotated channels, events, and electrodes
% with SOZ annotated

%% Set paths, get preprocessed data for one subject

localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('functions');

%% Load BB output to check

subjects = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17'};

ss = 6;
sub_label = subjects{ss};
ses_label = 'ieeg01';

outdir = fullfile(localDataPath.output,'derivatives','preproc_car',['sub-' sub_label]);

load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_ieeg.mat']), 'tt', 'srate', 'Mdata', 'eventsST', 'all_channels');
load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCARBB_ieeg.mat']), 'tt', 'srate', 'Mbb', 'eventsST', 'all_channels');
readtable(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_events.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');

cm = ieeg_kjmloccolormap;


%% Normalize bb power and signal per run

% Indicate the interval for baseline, used in normalization
norm_int = find(tt>-.2 & tt<0);

% Baseline subtract signal
Mdata(strcmp(all_channels.type, 'SEEG'),:,:) = minus(Mdata(strcmp(all_channels.type, 'SEEG'),:,:), mean(Mdata(strcmp(all_channels.type, 'SEEG'),norm_int,:),2));

% Initialize normalized log power of BB
Mbb_norm = Mbb;
Mbb_norm(strcmp(all_channels.type, 'SEEG'),:) = log10(Mbb(strcmp(all_channels.type, 'SEEG'),:)); 

% Normalize per run
for run_idx = 1:max(eventsST.tasknumber)
    this_run = find(eventsST.tasknumber==run_idx); % out of 1500
    
    % find pre-stim events with 'good' status
    trials_norm = find(ismember(eventsST.pre_status,'good') & eventsST.tasknumber==run_idx);

    Mbb_norm(strcmp(all_channels.type, 'SEEG'),:,this_run) = minus(Mbb_norm(strcmp(all_channels.type, 'SEEG'),:,this_run),mean(Mbb_norm(strcmp(all_channels.type, 'SEEG'),norm_int,trials_norm),[2 3],'omitnan'));
end

clear Mbb

%% Calculate spectogram for one electrode 
 
el_name = 'ROC8';
el_nr = find(ismember(all_channels.name,el_name));

[S, f] = ieeg_getWaveletSpectrogram(squeeze(Mdata(el_nr,:,:)), srate, [1 175]);
S = log10(S);

% mean baseline power from -500:0 ms
S_base = mean(S(:,tt>-0.5 & tt<0,:),[3 2]);

% normalization wrt baseline is necessary (Power/baselinePower)
S_norm = S - repmat(S_base,1,size(S,2),size(S,3));

% downsample
% If X is a matrix, then resample treats each column as an independent channel.
% get size of downsampled data and time
[Y,tt_downsample] = resample(S_norm(:,:,1)',tt,120);
S_downsample = zeros(size(S_norm,1),length(tt_downsample),size(S_norm,3));
for kk = 1:size(S_norm,3) % trials
    [Y,Ty] = resample(S_norm(:,:,kk)',tt,120);
    S_downsample(:,:,kk) = Y';
end

% clip
S_downsample = S_downsample(:,tt_downsample>-0.5 & tt_downsample<1.5,:);
tt_downsample = tt_downsample(tt_downsample>-0.5 & tt_downsample<1.5);

S_downsample = single(S_downsample);
save(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_wavelet_el' el_name '.mat']),...
    'tt_downsample', 'S_downsample', 'f', 'el_name','-v7.3');

%% Plot broadband and spectogram

figure('Position',[0 0 400 600]),
subplot(3,1,1), hold on
plot(tt,mean(Mbb_norm(el_nr,:,:),3),'b')
xlabel('Time (s)'),ylabel('Broadband log power')
xlim([-.5 1.5])
yline(0)

subplot(3,1,2:3)
uimagesc(tt,f,mean(S_norm,3),[-.5 .5])
axis xy
xlim([-.5 1.5])
xlabel('Time (s)'),ylabel('Frequency (Hz)')
title([all_channels.name{el_nr}])
colormap(cm); a = colorbar;
ylabel(a,'log power change wrt baseline','FontSize',12,'Rotation',90)

%% Plot broadband and downsampled spectogram

figure('Position',[0 0 400 600]),
subplot(3,1,1), hold on
plot(tt,mean(Mbb_norm(el_nr,:,:),3),'b')
xlabel('Time (s)'),ylabel('Broadband log power')
xlim([-.5 1.5])
yline(0)

subplot(3,1,2:3)
uimagesc(tt_downsample,f,mean(S_downsample,3),[-.5 .5])
axis xy
xlim([-.5 1.5])
xlabel('Time (s)'),ylabel('Frequency (Hz)')
title([all_channels.name{el_nr}])
colormap(cm); a = colorbar;
ylabel(a,'log power change wrt baseline','FontSize',12,'Rotation',90)



