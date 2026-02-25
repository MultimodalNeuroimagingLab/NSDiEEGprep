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

%% load all NSD images 

im_dims = [425,425,3];

stimuliOrig = zeros(im_dims(1),im_dims(2),1000,im_dims(3));
imNames = cell(1000,1);
imPaths = dir(fullfile(localDataPath.output,'stimuli','shared1000','*.png'));
for kk = 1:length(imPaths)
    if mod(kk, 100) == 0, disp(['image ' int2str(kk) ' of ' int2str(length(imPaths))]), end
    this_im = imread(fullfile(imPaths(kk).folder,imPaths(kk).name));
    % im_gray = im2gray(this_im);
    stimuliOrig(:,:,kk,:) = this_im;
    imNames{kk,1} = imPaths(kk).name;
end

%% Load fmri data
% fMRI shared NSD (average across all subjects)
nsd_fmri = load(fullfile(localDataPath.fmri,'betaAvg_shared_allsub.mat'));
designFile = load(fullfile(localDataPath.fmri,'nsd_expdesign.mat'));
% make sure order is the same as iEEG
aa = zeros(1000,1);
for kk = 1:1000
    a_temp = extractBetween(imNames{kk},'_nsd','.png');
    aa(kk) = str2num(a_temp{1});
end
figure,plot(aa-designFile.sharedix')

%% reorder data according to images

eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);

bb_1000 = NaN(size(Mbb_norm,1),size(Mbb_norm,2),1000);
sig_1000 = NaN(size(Mdata,1),size(Mdata,2),1000);
bb_100 = NaN(size(Mbb_norm,1),size(Mbb_norm,2),100,6);
special100_idx = zeros(1000,1);
count_100 = 0;
for kk = 1:length(imNames) % 1000
    % get nsd_idx 
    this_shared_idx = extractBetween(imNames{kk},'shared','_nsd');
    this_shared_idx = str2num(this_shared_idx{1});
    
    % trials with this shared idx & good 
    these_trials = find(shared_idx==this_shared_idx & events_status(:,1)==0);
    
    % get data, average if multiple
    if ~isempty(these_trials)
        bb_1000(:,:,kk) = mean(Mbb_norm(:,:,these_trials),3);
        sig_1000(:,:,kk) = mean(Mdata(:,:,these_trials),3);
    end
    if length(these_trials)>1
        count_100 = count_100+1;
        special100_idx(kk) = count_100;
        bb_100(:,:,count_100,1:length(these_trials)) = Mbb_norm(:,:,these_trials);
    end
end


%%
%%
%% correlation over time per ROI
%%
%%

el_nr = find(ismember(all_channels.name,'ROC2'));
thisBB = squeeze(bb_1000(el_nr,:,:));

% fmri beta from one roi
hemi = 'r';
[v_Rosenke, v_Rosenke_label, cmap] = loadAtlas( localDataPath.fs, hemi, 'visf');
[v_BensonE, ~, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonE');
[v_BensonV, v_BensonV_label, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonV');
[v_Wang, v_Wang_label, ~] = loadAtlas( localDataPath.fs, hemi, 'Wang');
[v_HCP, v_HCP_label, ~] = loadAtlas( localDataPath.fs, hemi, 'HCP');

% fmri_beta = nsd_fmri.beta_avgR(v_Rosenke>=1 & v_Rosenke<=9,:);
% fmri_beta = nsd_fmri.beta_avgR(v_HCP==18 | v_HCP==22,:); % FFA | PIT 
% fmri_beta = nsd_fmri.beta_avgR(v_HCP==10 | v_HCP==11,:); % FEF
% fmri_beta = nsd_fmri.beta_avgR(v_BensonV==4,:); % hV4
% fmri_beta = nsd_fmri.beta_avgR(v_Wang==1,:); 

vertmask(1).idx = v_Wang==1;
vertmask(1).nr = length(find(vertmask(1).idx==1));
vertmask(2).idx = v_Wang==3;
vertmask(2).nr = length(find(vertmask(2).idx==1));
vertmask(3).idx = v_Wang==5;
vertmask(3).nr = length(find(vertmask(3).idx==1));
vertmask(4).idx = v_BensonV==4;
vertmask(4).nr = length(find(vertmask(4).idx==1));
vertmask(5).idx = v_Rosenke==8 | v_Rosenke==9; % places
vertmask(5).nr = length(find(vertmask(5).idx==1));

t_ints = [[-.020:0.020:0.500]' [0:0.020:0.520]'];
sum_stats = zeros(4,size(t_ints,1),length(vertmask));

for jj = 1:length(vertmask)
    fmri_beta = nsd_fmri.beta_avgR(vertmask(jj).idx==1,:); 

    corr_out = zeros(size(t_ints,1),size(fmri_beta,1));
    for ii = 1:size(t_ints,1)
        disp(['starting ' int2str(ii) ' of ' int2str(size(t_ints,1))])
        bb = mean(thisBB(tt>t_ints(ii,1) & tt<t_ints(ii,2),:),1)';
    
        for kk = 1:size(fmri_beta,1)
            if mod(kk,10000)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta,1))]), end
            R = faster_corr_mtrx([bb fmri_beta(kk,:)']);
            corr_out(ii,kk) = R(1,2);
        end
    end
    
    th = .3;
    sum_stats(:,:,jj) = ...
        [sum(corr_out>th,2)';...
        median(corr_out,2)';...
        mean(corr_out,2)';...
        max(corr_out')];
end       

figure('Position',[0 0 400 900]),
subplot(5,1,1)
plot(tt,mean(thisBB,2))
xlim([t_ints(1,1) t_ints(end,1)])
xlabel('time(s)'), ylabel('bb log power')

stat_name = {'# vertices with r>.3','median','mean','max',};
for ii_stat = 1:3    
    subplot(5,1,1+ii_stat)
    plot(t_ints(:,1),squeeze(sum_stats(ii_stat,:,:)))
    xlabel('time(s)'), ylabel(stat_name{ii_stat})
    xlim([t_ints(1,1) t_ints(end,1)])
end
subplot(5,1,2)
legend({'V1v','V2v','V3v','hV4','places'})

perc_vert = zeros(4,size(t_ints,1));
for v_area = 1:4
    perc_vert(v_area,:) = sum_stats(1,:,v_area) ./ vertmask(v_area).nr;   
end

subplot(5,1,5)
plot(t_ints(:,1),perc_vert)
xlabel('time(s)'), ylabel('% vertices with r>.3')
xlim([t_ints(1,1) t_ints(end,1)])

