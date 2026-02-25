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


%% Calculate spectograms
% 
el_name = 'ROC2';
el_nr = find(ismember(all_channels.name,el_name));

[S, f] = ieeg_getWaveletSpectrogram(squeeze(Mdata(el_nr,:,:)), srate, [1 175]);
S = log10(S);

% mean baseline power from -500:0 ms
S_base = mean(S(:,tt>-0.5 & tt<0,:),[3 2]);

% normalization wrt baseline is necessary (Power/baselinePower)
S_norm = S - repmat(S_base,1,size(S,2),size(S,3));

figure('Position',[0 0 300 300]),
uimagesc(tt,f,mean(S_norm,3),[-.3 .3])
axis xy
xlim([-.5 1.5])
xlabel('Time (s)'),ylabel('Frequency (Hz)')
title([all_channels.name{el_nr}])
colormap(cm); a = colorbar;
ylabel(a,'log power change wrt baseline','FontSize',12,'Rotation',90)

% reorder spectogram according to 1000 images
spect_1000 = NaN(size(S_norm,1),size(S_norm,2),length(imNames));
for kk = 1:length(imNames) % 1000
    % get nsd_idx 
    this_shared_idx = extractBetween(imNames{kk},'shared','_nsd');
    this_shared_idx = str2num(this_shared_idx{1});
    
    % trials with this shared idx & good 
    these_trials = find(shared_idx==this_shared_idx & events_status(:,1)==0);
    
    % get data, average if multiple
    if ~isempty(these_trials)
        spect_1000(:,:,kk) = mean(S_norm(:,:,these_trials),3);
    end
end

figure('Position',[0 0 300 300]),
uimagesc(tt,f,mean(spect_1000,3),[-.3 .3])
axis xy
xlim([-.5 1.5])
xlabel('Time (s)'),ylabel('Frequency (Hz)')
title([all_channels.name{el_nr}])
colormap(cm); a = colorbar;
ylabel(a,'log power change wrt baseline','FontSize',12,'Rotation',90)

%% correlate spectogram features with fmri

% right hemisphere
fmri_beta = nsd_fmri.beta_avgR;
 
ieeg_alpha = squeeze(mean(spect_1000(f>14 & f<20,tt>.3 & tt<.8,:),[1 2]));
ieeg_bb = squeeze(mean(spect_1000(f>80 & f<160,tt>.100 & tt<.120,:),[1 2]));

corr_alpha = NaN(size(fmri_beta,1),1);
corr_bb = NaN(size(fmri_beta,1),1);
for kk = 1:size(fmri_beta,1)
    if mod(kk,1000)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta,1))]), end
    R = faster_corr_mtrx([ieeg_alpha fmri_beta(kk,:)']);
    corr_alpha(kk) = R(1,2);
    R = faster_corr_mtrx([ieeg_bb fmri_beta(kk,:)']);
    corr_bb(kk) = R(1,2);
end


%%
% render and plot correlation BB and Alpha

% Load fs average surfaces
gI_L = gifti(fullfile(localDataPath.fs,'fsaverage','inflated.L.surf.gii'));
gI_R = gifti(fullfile(localDataPath.fs,'fsaverage','inflated.R.surf.gii'));
sulcal_labels_L = read_curv(fullfile(localDataPath.fs,'fsaverage','surf','lh.sulc'));
sulcal_labels_R = read_curv(fullfile(localDataPath.fs,'fsaverage','surf','rh.sulc'));

for hh = 2

    if hh==1
        hemi = 'l';
        g = gI_L;
        sulcal_labels = sulcal_labels_L;
        views_plot = {[40,-30],[-45,-10],[-90,20]};
    elseif hh==2
        hemi = 'r';
        g = gI_R;
        sulcal_labels = sulcal_labels_R;
        views_plot = {[-40,-30],[45,-10],[90,20]};
    end

    % load the color map
    cmap = cmapsign4(200);

    % make a plot with electrode dots
    for vv = 3%:length(views_plot)
        v_d = [views_plot{vv}(1),views_plot{vv}(2)];
        
        figure
        vert_amp = corr_bb;
        vert_amp = vert_amp/max(abs(vert_amp)); % set max to 1
        vert_amp(abs(vert_amp)<.1) = NaN;
        vert_amp = 100+round(vert_amp*100); % multiply by 100 and add 100

        tH = ieeg_RenderGiftiLabels(g,vert_amp,cmap,[],sulcal_labels);
        ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')

        els = g.vertices(loc_mni305.vertex_fsaverage(loc_ind),:);
        a_offset=.1*max(abs(els(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els = els + repmat(a_offset,size(els,1),1);
        plot3(els(:,1),els(:,2),els(:,3),'o','MarkerFaceColor', [.5 .5 1],'MarkerEdgeColor',[0 0 .3],'LineWidth',2,'MarkerSize',10)

    end
%     close all
end

%% stack images for alpha and plot over time
thisAA = squeeze(mean(spect_1000(f>8 & f<12,:,:),[1]));

% get MNI1305 position of electrode and check position in atlas
loc_mni305 = readtable(fullfile(localDataPath.input, ['sub-' sub_label],'ses-ieeg01','ieeg',['sub-' sub_label '_ses-ieeg01_space-MNI305_electrodes.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');
loc_ind = find(ismember(loc_mni305.name,el_name));

hemi = lower(el_name(1));

% get fmri for correct hemisphere
if isequal(hemi,'r')
    fmri_beta_all = nsd_fmri.beta_avgR;
elseif isequal(hemi,'l')
    fmri_beta_all = nsd_fmri.beta_avgL;
end

% check which visual area this electrode is positioned in
[v_Rosenke, v_Rosenke_label, cmap] = loadAtlas( localDataPath.fs, hemi, 'visf');
[v_BensonV, v_BensonV_label, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonV');
[v_Wang, v_Wang_label, ~] = loadAtlas( localDataPath.fs, hemi, 'Wang');
[v_HCP, v_HCP_label, ~] = loadAtlas( localDataPath.fs, hemi, 'HCP');
[v_BensonE, ~, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonE');

atlas_ind = 4;
if atlas_ind==1
    el_atlasLabels = v_Rosenke(loc_mni305.vertex_fsaverage(loc_ind));
    v_Rosenke_label(el_atlasLabels)
    roi_inds = find(v_Rosenke==el_atlasLabels);
elseif atlas_ind==2
    el_atlasLabels = v_BensonV(loc_mni305.vertex_fsaverage(loc_ind));
    roi_inds = find(v_BensonV==el_atlasLabels);
    v_BensonV_label(el_atlasLabels)
elseif atlas_ind==3
    el_atlasLabels = v_Wang(loc_mni305.vertex_fsaverage(loc_ind));
    roi_inds = find(v_Wang==el_atlasLabels);
    v_Wang_label(el_atlasLabels)
elseif atlas_ind==4    
    el_atlasLabels = v_HCP(loc_mni305.vertex_fsaverage(loc_ind));
    roi_inds = find(v_HCP==el_atlasLabels);    
    v_HCP_label(el_atlasLabels)
end

if el_atlasLabels==0
    disp('no atlas label, return')
    return
end

% fMRI beta within this roi
fmri_beta = fmri_beta_all(roi_inds,:);

% correlate fMRI/iEEG for one time-bin around maximum iEEG
t_int = find(tt>.2 & tt<.8); %+/- 25 ms
av_bb = mean(thisAA(t_int,:),1)';

% use the first 100 images to establish vertices that correlate
corr_out = zeros(1,size(fmri_beta,1));
for kk = 1:size(fmri_beta,1)
    if mod(kk,10000)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta,1))]), end
    R = faster_corr_mtrx([av_bb(special100_idx==0) fmri_beta(kk,special100_idx==0)']);
    corr_out(kk) = R(1,2);
end

corr_th = -.10;
disp([int2str(length(find(corr_out<corr_th))) ' voxels with r<' num2str(corr_th) ])

% vertices that correlate somewhat in this roi
corr_vert_indices = roi_inds(corr_out<corr_th);

% for plotting the fMRI images (testing set)
stim_900 = stimuliOrig(:,:,special100_idx==0,:);
[fmri_sort,ii_sort] = sort(mean(fmri_beta(corr_out<corr_th,special100_idx==0),1),'descend');
ff = makeimagestack_wrapper(stim_900(:,:,ii_sort(1:50),:),0,1,[10 5]);

% now we use the other 900 images to get image specificy
fmri_avbeta = fmri_beta(corr_out<corr_th,special100_idx>0); % average across correlated vertices
bb_tempcorr = thisAA(:,special100_idx>0);

figure('Position',[0 0 300 270])

subplot(2,2,1)
% plot(tt,mean(thisBB,2),'k')
ieeg_plotCurvConf(tt,thisAA',[0 0 1],.5)
xlim([-0.1 1.2]),yline(0)
set(gca,'YTick',[0:.1:1],'XTick',[0:.4:1],'fontname','arial')
title(el_name)
box off

subplot(2,2,[2 4])
imshow(uint8(ff))
title('predicted preferred images')
set(gca,'fontname','arial')

% now we show at which timepoint these two vector point in the same direction
subplot(2,2,3),hold on
% cross_proj = (bb_tempcorr*fmri_avbeta');
% plot(tt,cross_proj)
% rr = corr(bb_tempcorr',mean(fmri_avbeta,1)');
% plot(tt,rr,'k')
rr = corr(bb_tempcorr',fmri_avbeta');
fill([tt(1) tt(end) tt(end) tt(1)],[quantile(mean(rr(tt<0,:),2),.975) quantile(mean(rr(tt<0,:),2),.975) quantile(mean(rr(tt<0,:),2),.025) quantile(mean(rr(tt<0,:),2),.025)],...
    [1 .9 .8],'EdgeColor',[1 .9 .8],'FaceAlpha',.5)
plot(tt,mean(rr,2),'k','LineWidth',2)
% ieeg_plotCurvConf(tt,rr',[0 0 1],.5)
xlim([-0.1 1.2]),ylim([-.6 .6]),yline(0,'k')
set(gca,'YTick',[-0.4:.2:1],'XTick',[0:.4:1],'fontname','arial')
box off

%% plot alpha VS bb

t_int1 = find(tt>.2 & tt<.8); 
thisAA = squeeze(mean(spect_1000(f>15 & f<24,t_int1,:),[1 2]));
t_int2 = find(tt>.1 & tt<.5); %+/- 25 ms
thisBB = squeeze(mean(spect_1000(f>60 & f<180,t_int2,:),[1 2]));

% use the first 100 images to establish vertices that correlate
corr_out = zeros(1,size(fmri_beta,1));
for kk = 1:size(fmri_beta,1)
    if mod(kk,10000)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta,1))]), end
    R = faster_corr_mtrx([thisBB(special100_idx==0) fmri_beta(kk,special100_idx==0)']);
    corr_out(kk) = R(1,2);
end

corr_th = 0.10;
disp([int2str(length(find(corr_out>corr_th))) ' voxels with r>' num2str(corr_th) ])

% vertices that correlate somewhat in this roi
corr_vert_indices = roi_inds(corr_out>corr_th);

% now we use the other 900 images to get image specificy
% fmri_avbeta = fmri_beta(corr_out>corr_th,:); % average across correlated vertices
fmri_avbeta = mean(fmri_beta,1); % average across correlated vertices


regstats(fmri_avbeta,[thisBB thisAA],'linear')


