%% This is the preprocessing script for the clean NSD analysis.
% Requires each subject to have annotated channels, events, and electrodes
% with SOZ annotated

%% Set paths, get preprocessed data for one subject

localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('functions');

%% Load BB output to check

subjects = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17'};

ss = 17;
sub_label = subjects{ss};
ses_label = 'ieeg01';

outdir = fullfile(localDataPath.output,'derivatives','preproc_car',['sub-' sub_label]);

% load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_ieeg.mat']), 'tt', 'srate', 'Mdata', 'eventsST', 'all_channels');
load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCARBB_ieeg.mat']), 'tt', 'srate', 'Mbb', 'eventsST', 'all_channels');
readtable(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_events.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');

cm = ieeg_kjmloccolormap;


%% Normalize bb power and signal per run

% Indicate the interval for baseline, used in normalization
norm_int = find(tt>-.2 & tt<0);

% % Baseline subtract signal
% Mdata(strcmp(all_channels.type, 'SEEG'),:,:) = minus(Mdata(strcmp(all_channels.type, 'SEEG'),:,:), mean(Mdata(strcmp(all_channels.type, 'SEEG'),norm_int,:),2));

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

%% load all NSD images and fMRI data

%%%% Load NSD images
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

%%%% Load fmri data
% fMRI shared NSD (average across all subjects)
nsd_fmri = load(fullfile(localDataPath.fmri,'betaAvg_shared_allsub.mat'));
designFile = load(fullfile(localDataPath.fmri,'nsd_expdesign.mat'));
% make sure order is the same as iEEG
aa = zeros(1000,1);
for kk = 1:1000
    a_temp = extractBetween(imNames{kk},'_nsd','.png');
    aa(kk) = str2num(a_temp{1});
end
figure,plot(aa-designFile.sharedix'), title('should be all zeros')

%% reorder data according to images

eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);

bb_1000 = NaN(size(Mbb_norm,1),size(Mbb_norm,2),1000);
% sig_1000 = NaN(size(Mdata,1),size(Mdata,2),1000);
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
%         sig_1000(:,:,kk) = mean(Mdata(:,:,these_trials),3);
    end
    if length(these_trials)>1
        count_100 = count_100+1;
        special100_idx(kk) = count_100;
        bb_100(:,:,count_100,1:length(these_trials)) = Mbb_norm(:,:,these_trials);
    end
end


%% rank broadband over time bins

% sub-05: RPS 2 - body
% sub-05: RPS 3 - food

el_name = 'RT5';

el_nr = find(ismember(all_channels.name,el_name));
thisBB = squeeze(bb_1000(el_nr,:,:));
figure
subplot(3,3,1:3)
plot(tt,mean(thisBB,2))

t_int = [[-.1:.05:1]' [-.05:.05:1.05]'];
xlim([min(t_int(:)) max(t_int(:))])

f_out = [];
for kk = 1:size(t_int,1)
    [bb_sort,ii_sort] = sort(mean(thisBB(tt>t_int(kk,1) & tt<t_int(kk,2),:),1),'descend');
    f = makeimagestack_wrapper(stimuliOrig(:,:,ii_sort(1:20),:),0,1,[20 1]);
    f_out = cat(2,f_out,256*ones(size(f,1),20,size(f,3)),f);
end
subplot(3,3,4:9)
imshow(uint8(f_out))



%% find correlated fMRI vertices within the area where the electrode is located
% then we plot over time when this selectivity arises

el_name = 'LOC4';
el_nr = find(ismember(all_channels.name,el_name));
thisBB = squeeze(bb_1000(el_nr,:,:));

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
aa = mean(thisBB(:,special100_idx==0),2);
[~,max_bb_tint] = max(aa(tt>0 & tt<.8));
max_bb_tint = find(tt>0,1)+max_bb_tint;
t_int = max_bb_tint-srate/40:max_bb_tint+srate/40; %+/- 25 ms
av_bb = mean(thisBB(t_int,:),1)';

% use the first 100 images to establish vertices that correlate
corr_out = zeros(1,size(fmri_beta,1));
for kk = 1:size(fmri_beta,1)
    if mod(kk,10000)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta,1))]), end
    R = faster_corr_mtrx([av_bb(special100_idx==0) fmri_beta(kk,special100_idx==0)']);
    corr_out(kk) = R(1,2);
end

corr_th = .10;
disp([int2str(length(find(corr_out>corr_th))) ' voxels with r>' num2str(corr_th) ])

% vertices that correlate somewhat in this roi
corr_vert_indices = roi_inds(corr_out>corr_th);

% for plotting the fMRI images (testing set)
stim_900 = stimuliOrig(:,:,special100_idx==0,:);
[fmri_sort,ii_sort] = sort(mean(fmri_beta(corr_out>corr_th,special100_idx==0),1),'descend');
f = makeimagestack_wrapper(stim_900(:,:,ii_sort(1:50),:),0,1,[10 5]);

% now we use the other 900 images to get image specificy
fmri_avbeta = fmri_beta(corr_out>corr_th,special100_idx>0); % average across correlated vertices
bb_tempcorr = thisBB(:,special100_idx>0);

figure('Position',[0 0 300 270])

subplot(2,2,1)
% plot(tt,mean(thisBB,2),'k')
ieeg_plotCurvConf(tt,thisBB',[0 0 1],.5)
xlim([-0.1 1.2]),yline(0)
set(gca,'YTick',[0:.1:1],'XTick',[0:.4:1],'fontname','arial')
title(el_name)
box off

subplot(2,2,[2 4])
imshow(uint8(f))
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
xlim([-0.1 1.2]),ylim([-.2 .6]),yline(0,'k')
set(gca,'YTick',[-0.4:.2:1],'XTick',[0:.4:1],'fontname','arial')
box off

 
print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','VSS',...
    ['sub-' sub_label '_' el_name '_bb_corr_roi']))
print('-depsc','-painters','-r300',fullfile(localDataPath.input,'derivatives','VSS',...
    ['sub-' sub_label '_' el_name '_bb_corr_roi']))


%% render and plot atlas, correlated voxels and electrode

% Load fs average surfaces
gI_L = gifti(fullfile(localDataPath.fs,'fsaverage','inflated.L.surf.gii'));
gI_R = gifti(fullfile(localDataPath.fs,'fsaverage','inflated.R.surf.gii'));
sulcal_labels_L = read_curv(fullfile(localDataPath.fs,'fsaverage','surf','lh.sulc'));
sulcal_labels_R = read_curv(fullfile(localDataPath.fs,'fsaverage','surf','rh.sulc'));

if isequal(hemi,'l') 
    g = gI_L;
    sulcal_labels = sulcal_labels_L;
    views_plot = {[40,-30],[-45,-10],[-90,20],[0,-90]};
elseif isequal(hemi,'r') 
    g = gI_R;
    sulcal_labels = sulcal_labels_R;
    views_plot = {[-40,-30],[45,-10],[90,20],[0,-90]};
end

% load the color map
% cmap = cmapsign4(200);
% cmap = lines(180);
% cmap = cmap.^.2;
cmap = brewermap(180,'Set3').^.5;
cmap = [cmap; [.2 .2 .2]];

vert_amp_corr = NaN(size(sulcal_labels));
corr_plot = corr_out;
corr_plot(abs(corr_plot)<.1) = NaN;
% vert_amp_corr(atlas_def) = corr_plot;
% vert_amp_corr(roi_inds) = corr_plot;

% make a plot with electrode dots
for vv = 4%:length(views_plot)
    v_d = [views_plot{vv}(1),views_plot{vv}(2)];
    
    % plot 
    figure
%     vert_amp = vert_amp_corr;
%     vert_amp = vert_amp/max(abs(vert_amp)); % set max to 1
%     vert_amp = 100+round(vert_amp*100); % multiply by 100 and add 100
    
    vert_amp = v_HCP;
    vert_amp(vert_amp>180) = 0;
    vert_amp(corr_vert_indices) = 181;

    tH = ieeg_RenderGiftiLabels(g,vert_amp,cmap,[],sulcal_labels);
    ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle
    hold on, 
    els = g.vertices(loc_mni305.vertex_fsaverage(loc_ind),:);
    
    a_offset=.1*max(abs(els(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
    els = els + repmat(a_offset,size(els,1),1);

    plot3(els(:,1),els(:,2),els(:,3),'o','MarkerFaceColor', [.5 .5 1],'MarkerEdgeColor',[0 0 .3],'LineWidth',2,'MarkerSize',10)

    set(gcf,'PaperPositionMode','auto')

    print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','VSS',...
        ['sub-' sub_label '_' el_name '_bb_corr_roi_render' int2str(v_d(1)) int2str(v_d(2))]))

end
%     close all


%%
%% Other ways to do the same:
%%
%% find correlated fMRI vertices within the area where the electrode is located
% then we plot over time when this selectivity arises

% el_name = 'RPO2';
% el_nr = find(ismember(all_channels.name,el_name));
% thisBB = squeeze(bb_1000(el_nr,:,:));
% 
% % get MNI1305 position of electrode and check position in atlas
% loc_mni305 = readtable(fullfile(localDataPath.input, ['sub-' sub_label],'ses-ieeg01','ieeg',['sub-' sub_label '_ses-ieeg01_space-MNI305_electrodes.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');
% loc_ind = find(ismember(loc_mni305.name,el_name));
% 
% hemi = lower(el_name(1));
% 
% % get fmri for correct hemisphere
% if isequal(hemi,'r')
%     fmri_beta_all = nsd_fmri.beta_avgR;
% elseif isequal(hemi,'l')
%     fmri_beta_all = nsd_fmri.beta_avgL;
% end
% 
% % check which visual area this electrode is positioned in
% [v_Rosenke, v_Rosenke_label, cmap] = loadAtlas( localDataPath.fs, hemi, 'visf');
% [v_BensonV, v_BensonV_label, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonV');
% [v_Wang, v_Wang_label, ~] = loadAtlas( localDataPath.fs, hemi, 'Wang');
% [v_HCP, v_HCP_label, ~] = loadAtlas( localDataPath.fs, hemi, 'HCP');
% [v_BensonE, ~, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonE');
% 
% atlas_ind = 4;
% if atlas_ind==1
%     el_atlasLabels = v_Rosenke(loc_mni305.vertex_fsaverage(loc_ind));
%     v_Rosenke_label(el_atlasLabels)
%     roi_inds = find(v_Rosenke==el_atlasLabels);
% elseif atlas_ind==2
%     el_atlasLabels = v_BensonV(loc_mni305.vertex_fsaverage(loc_ind));
%     roi_inds = find(v_BensonV==el_atlasLabels);
%     v_BensonV_label(el_atlasLabels)
% elseif atlas_ind==3
%     el_atlasLabels = v_Wang(loc_mni305.vertex_fsaverage(loc_ind));
%     roi_inds = find(v_Wang==el_atlasLabels);
%     v_Wang_label(el_atlasLabels)
% elseif atlas_ind==4    
%     el_atlasLabels = v_HCP(loc_mni305.vertex_fsaverage(loc_ind));
%     roi_inds = find(v_HCP==el_atlasLabels);    
%     v_HCP_label(el_atlasLabels)
% end
% 
% if el_atlasLabels==0
%     disp('no atlas label, return')
%     return
% end
% 
% % fMRI beta within this roi
% fmri_beta = fmri_beta_all(roi_inds,:);
% 
% % correlate fMRI/iEEG for one time-bin around maximum iEEG
% aa = mean(thisBB(:,1:100),2);
% [~,max_bb_tint] = max(aa(tt>0 & tt<.8));
% max_bb_tint = find(tt>0,1)+max_bb_tint;
% t_int = max_bb_tint-srate/40:max_bb_tint+srate/40; %+/- 25 ms
% av_bb = mean(thisBB(t_int,:),1)';
% 
% % use the first 100 images to establish vertices that correlate
% corr_out = zeros(1,size(fmri_beta,1));
% for kk = 1:size(fmri_beta,1)
%     if mod(kk,10000)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta,1))]), end
%     R = faster_corr_mtrx([av_bb(1:100) fmri_beta(kk,1:100)']);
%     corr_out(kk) = R(1,2);
% end
% 
% corr_th = .1;
% disp([int2str(length(find(corr_out>corr_th))) ' voxels with r>' num2str(corr_th) ])
% 
% % vertices that correlate somewhat in this roi
% corr_vert_indices = roi_inds(corr_out>corr_th);
% 
% % now we use the other 900 images to get image specificy
% fmri_avbeta = fmri_beta(corr_out>corr_th,101:1000); % average across correlated vertices
% bb_tempcorr = thisBB(:,101:1000);
% 
% % we plot the fMRI images
% stim_900 = stimuliOrig(:,:,101:1000,:);
% [fmri_sort,ii_sort] = sort(mean(fmri_avbeta,1),'descend');
% f = makeimagestack_wrapper(stim_900(:,:,ii_sort(1:80),:),0,1,[8 10]);
% 
% 
% figure('Position',[0 0 300 600])
% 
% subplot(2,2,1)
% % plot(tt,mean(thisBB,2),'k')
% ieeg_plotCurvConf(tt,thisBB',[0 0 1],.5)
% xlim([-0.1 1.1]),ylim([-.1 .3]),yline(0)
% 
% % now we show at which timepoint these two vector point in the same direction
% subplot(2,2,2)
% % cross_proj = (bb_tempcorr*fmri_avbeta');
% % plot(tt,cross_proj)
% % rr = corr(bb_tempcorr',mean(fmri_avbeta,1)');
% % plot(tt,rr,'k')
% 
% rr = corr(bb_tempcorr',fmri_avbeta');
% ieeg_plotCurvConf(tt,rr',[0 0 1],.5)
% xlim([-0.1 1.1]),ylim([-.1 .3]),yline(0)
% 
% subplot(2,1,2)
% imshow(uint8(f))
% 
% % print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','VSS','prelim',...
% %     ['sub-' sub_label '_' el_name '_bb_corr_roi']))
% % print('-depsc','-painters','-r300',fullfile(localDataPath.input,'derivatives','VSS','prelim',...
% %     ['sub-' sub_label '_' el_name '_bb_corr_roi']))
% 
% 
% 
%%

el_name = 'RT5';
el_nr = find(ismember(all_channels.name,el_name));
thisBB = squeeze(bb_1000(el_nr,:,:));

% get MNI1305 position of electrode and check position in atlas
loc_mni305 = readtable(fullfile(localDataPath.input, ['sub-' sub_label],'ses-ieeg01','ieeg',['sub-' sub_label '_ses-ieeg01_space-MNI305_electrodes.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');
loc_ind = find(ismember(loc_mni305.name,el_name));

% check which visual area this electrode is positioned in
hemi = lower(el_name(1));
[v_Rosenke, v_Rosenke_label, cmap] = loadAtlas( localDataPath.fs, hemi, 'visf');
[v_BensonE, ~, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonE');
[v_BensonV, v_BensonV_label, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonV');
[v_Wang, v_Wang_label, ~] = loadAtlas( localDataPath.fs, hemi, 'Wang');
[v_HCP, v_HCP_label, ~] = loadAtlas( localDataPath.fs, hemi, 'HCP');

% restrict to areas of atlas:
% atlas_def = find(v_Rosenke>0 & v_Rosenke<=9);
atlas_def = find(v_HCP>120 & v_HCP<139);
% atlas_def = find(v_Wang==18)
% atlas_def = find(v_BensonV>0);

% correlate fMRI/iEEG for one time-bin
t_int = tt>0.1 & tt<.8;
av_bb = mean(thisBB(t_int,:),1)';

% get fMRI for correct area
if isequal(hemi,'r')
    fmri_beta = nsd_fmri.beta_avgR(atlas_def,:);
elseif isequal(hemi,'l')
    fmri_beta = nsd_fmri.beta_avgL(atlas_def,:);
end

% use the first 100 images to establish vertices that correlate
corr_out = zeros(1,size(fmri_beta,1));
for kk = 1:size(fmri_beta,1)
    if mod(kk,10000)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta,1))]), end
    R = faster_corr_mtrx([av_bb(special100_idx==0) fmri_beta(kk,special100_idx==0)']);
    corr_out(kk) = R(1,2);
end

corr_th = .10; % .15
length(find(corr_out>corr_th))

% now we use the repeated  100 images to get image specificy
fmri_avbeta = mean(fmri_beta(corr_out>corr_th,special100_idx>0),1); % average across correlated vertices
bb_tempcorr = thisBB(:,special100_idx>0);

% we plot the fMRI images
stim_900 = stimuliOrig(:,:,special100_idx>0,:);
[fmri_sort,ii_sort] = sort(fmri_avbeta,'descend');
f = makeimagestack_wrapper(stim_900(:,:,ii_sort(1:90),:),0,1,[9 10]);

figure
% now we show at which timepoint these two vector point in the same direction
subplot(2,2,1)
plot(tt,mean(thisBB,2))

subplot(2,2,2)
% cross_proj = (bb_tempcorr*fmri_avbeta');
% plot(tt,cross_proj)
rr = corr(bb_tempcorr',fmri_avbeta');
plot(tt,rr)

subplot(2,1,2)
imshow(uint8(f))

%%
%%

