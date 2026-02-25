%% This is the preprocessing script for the clean NSD analysis.
% Requires each subject to have annotated channels, events, and electrodes
% with SOZ annotated

%% Set paths, get preprocessed data for one subject

localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('functions');

%% Load BB output to check

subjects = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19'};

ss = 19;
sub_label = subjects{ss};
ses_label = 'ieeg01';

outdir = fullfile(localDataPath.output,'derivatives','preproc_car',['sub-' sub_label]);

% load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_ieeg.mat']), 'tt', 'srate', 'Mdata', 'eventsST', 'all_channels');
load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCARBB_ieeg.mat']), 'tt', 'srate', 'Mbb', 'eventsST', 'all_channels');
readtable(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_events.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');

cm = ieeg_kjmloccolormap;

%% order NSD color images according to data
 
im_dims = [425,425,3];

% load all NSD images
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

% reorder images according to ieeg data
stimuliOrd = NaN(im_dims(1),im_dims(2),height(eventsST),im_dims(3)); % events X resXresXcol
for kk = 1:height(eventsST)% event number
    if ~contains(eventsST.stim_file{kk},{'z1_','z2_','z3_'}) % not test image
        % get name of image that was shown during event
        thisEventImName = extractBetween(eventsST.stim_file{kk},'shared','_prepped1.png');
        thisEventImName = ['shared' thisEventImName{1} '.png'];

        % which number if this image in the imNames (filtered NSD images)
        imListNr = find(ismember(imNames,thisEventImName));
        
        % now we can put this image number in reordered_stimuli
        stimuliOrd(:,:,kk,:) = stimuliOrig(:,:,imListNr,:);
    end
end


%% Normalize bb power per run

% Initialize normalized log power of BB
Mbb_norm = log10(Mbb); 

% Indicate the interval for baseline, used in normalization
norm_int = find(tt>-.2 & tt<0);

% Normalize per run
for run_idx = 1:max(eventsST.tasknumber)
    this_run = find(eventsST.tasknumber==run_idx); % out of 1500
    
    % find pre-stim events with 'good' status
    trials_norm = find(ismember(eventsST.pre_status,'good') & eventsST.tasknumber==run_idx);

    Mbb_norm(:,:,this_run) = minus(Mbb_norm(:,:,this_run),mean(Mbb_norm(:,norm_int,trials_norm),[2 3],'omitnan'));
end


%% rank broadband over time bins

% sub-05: RPS 2 - body
% sub-05: RPS 3 - food

el_nr = find(ismember(all_channels.name,'LOC2'));
thisBB = squeeze(Mbb_norm(el_nr,:,ismember(eventsST.trial_type,'stim')));
thisStim = stimuliOrd(:,:,ismember(eventsST.trial_type,'stim'),:);

figure
subplot(2,1,1)
plot(tt,mean(thisBB,2))

subplot(2,1,2)
[bb_sort,ii_sort] = sort(mean(thisBB(tt>.050 & tt<.080,:)),'descend');
f = makeimagestack_wrapper(thisStim(:,:,ii_sort(1:50),:),0,1,[5 10]);
imshow(uint8(f))

%% SVD repeats

el_nr = find(ismember(all_channels.name,'LOC4'));
thisBB = squeeze(Mbb_norm(el_nr,:,:));

t_int = [-0.05 1.2];

% [u,s,v] = svd(thisBB(tt>0 & tt<1.2,:), 'econ');

eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);
[shared_idx_repeats] = shared_idx(nsd_repeats==1); % 100 images

% get 100 images from repeats
stimuliRepeats = stimuliOrd(:,:,nsd_repeats==1,:);

% initialize matrix for bb
repeats_bb = NaN(length(shared_idx_repeats),size(thisBB,1),6); % 100 X time X 6

for kk = 1:length(shared_idx_repeats)
    these_trials = find(shared_idx==shared_idx_repeats(kk));    % for this repeat, find the correct 6 trial numbers out of the 1500 and get the image and the data
    repeats_bb(kk,:,1:length(these_trials)) = thisBB(:,these_trials); 
end

% split half to get SVD
train_set = mean(repeats_bb(:,tt>t_int(1) & tt<t_int(2),1:2:end),3);

[u,s,v] = svd(train_set, 'econ');

% test on other half with increasing number of components
test_set = mean(repeats_bb(:,tt>t_int(1) & tt<t_int(2),2:2:end),3);

nr_comp = 10;
r_pred = zeros(nr_comp,1);

for kk = 1:nr_comp
    pred = u(:,1:kk)*s(1:kk,1:kk)*v(:,1:kk)';
    
    % R^2
    SSres = sum((test_set(:)-pred(:)).^2);
    SStot = sum((test_set(:)-mean(test_set(:))).^2);
    r_pred(kk) = 1 - SSres/SStot;
end

figure,
subplot(3,1,1)
plot(tt,mean(thisBB,2))
subplot(3,1,2)
plot(r_pred)
% ylim([0 1])
subplot(3,1,3)
plot(tt(tt>t_int(1) & tt<t_int(2)),v(:,1:2))


%% SVD non-repeats and test

el_nr = find(ismember(all_channels.name,'LOC4'));
thisBB = squeeze(Mbb_norm(el_nr,:,:));

t_int = [-0.05 1.5];

eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);
[shared_idx_repeats] = shared_idx(nsd_repeats==1); % 100 images

% get 100 images from repeats
stimuliRepeats = stimuliOrd(:,:,nsd_repeats==1,:);

% initialize matrix for bb
repeats_bb = NaN(length(shared_idx_repeats),size(thisBB,1),6); % 100 X time X 6

for kk = 1:length(shared_idx_repeats)
    these_trials = find(shared_idx==shared_idx_repeats(kk));    % for this repeat, find the correct 6 trial numbers out of the 1500 and get the image and the data
    repeats_bb(kk,:,1:length(these_trials)) = thisBB(:,these_trials); 
end

nonrepeats_set = thisBB(tt>t_int(1) & tt<t_int(2),nsd_repeats==0)';
[u,s,v] = svd(nonrepeats_set, 'econ');

% training data for weights
train_set = mean(repeats_bb(:,tt>t_int(1) & tt<t_int(2),1:2:end),3);

% test on other half with increasing number of components
test_set = mean(repeats_bb(:,tt>t_int(1) & tt<t_int(2),2:2:end),3);

nr_comp = 10;
r_pred = zeros(nr_comp,1);

for kk = 1:nr_comp
    % prediction from training set
    pc_weights = train_set/v(:,1:kk)'; 
    pred = pc_weights*v(:,1:kk)';
    
    % R^2
    SSres = sum((test_set(:)-pred(:)).^2);
    SStot = sum((test_set(:)-mean(test_set(:))).^2);
    r_pred(kk) = 1 - SSres/SStot;

end

figure,
subplot(3,1,1)
plot(tt,mean(thisBB,2))
subplot(3,1,2)
plot(r_pred)
% ylim([0 1])
subplot(3,1,3)
plot(tt(tt>t_int(1) & tt<t_int(2)),v(:,1:3))


%%
%% Find repeated images, calculate SNR
%%

eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);
[shared_idx_repeats] = shared_idx(nsd_repeats==1); % 100 images

% get 100 images from repeats
stimuliRepeats = stimuliOrd(:,:,nsd_repeats==1,:);
    
all_chan_snr = NaN(size(Mbb_norm,1),1);
t_avg = tt>0.2 & tt<0.3;
% t_avg = tt>0.4 & tt<0.7;

% el_nr = find(ismember(all_channels.name,'ROC5'));
el_nr = find(ismember(all_channels.name,'RPS3'));
thisBB = squeeze(Mbb_norm(el_nr,:,:));

bb_strength = squeeze(mean(Mbb_norm(el_nr,t_avg==1,:),2));

repeats_bb_strength = cell(length(shared_idx_repeats),1);
av_bb_repeats = NaN(length(shared_idx_repeats),1);
for kk = 1:length(shared_idx_repeats)
    these_trials = find(shared_idx==shared_idx_repeats(kk));    % for this repeat, find the correct 6 trial numbers out of the 1500 and get the image and the data
    repeats_bb_strength{kk} = bb_strength(these_trials); 
    av_bb_repeats(kk) = mean(bb_strength(these_trials)); 
end

[NCSNR, p, NCSNRNull] = estimateNCSNR(repeats_bb_strength, 1000);

% stimuliRepeats
% av_bb_repeats

figure
subplot(1,2,1)
[bb_sort,ii_sort] = sort(av_bb_repeats,'descend');
f = makeimagestack_wrapper(stimuliRepeats(:,:,ii_sort(1:20),:));
imshow(uint8(f))
subplot(1,2,2)
[bb_sort,ii_sort] = sort(av_bb_repeats,'ascend');
f = makeimagestack_wrapper(stimuliRepeats(:,:,ii_sort(1:20),:));
imshow(uint8(f))


%% 
%% for just one electrode
%% Find repeated images, calculate SNR

eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);

all_chan_snr = NaN(size(Mbb_norm,1),1);
t_avg = tt>0.1 & tt<.5;

for el_nr = 1:size(Mbb_norm,1)
    
    if ismember(all_channels.type(el_nr),'SEEG') && all_channels.status(el_nr)==1
        bb_strength = squeeze(mean(Mbb_norm(el_nr,t_avg==1,:),2));
        
        all_repeats = find(nsd_repeats>0);
        shared_idx_repeats = unique(shared_idx(all_repeats)); % 100 images
        repeats_bb_strength = cell(length(shared_idx_repeats),1);
        for kk = 1:length(shared_idx_repeats)
            these_trials = find(shared_idx==shared_idx_repeats(kk));    % for this repeat, find the correct 6 trial numbers out of the 1500 and get the image and the data
            repeats_bb_strength{kk} = bb_strength(these_trials); 
        end
        
        [NCSNR, p, NCSNRNull] = estimateNCSNR(repeats_bb_strength, 1000);
        all_chan_snr(el_nr) = NCSNR;
    end
end



%%
%% Calculate spectograms
%%

el_nr = find(ismember(all_channels.name,'LOC4'));

% data_spect = minus(Mdata(el_nr,:,:),mean(Mdata(el_nr,:,:),3));
data_spect = minus(Mdata(el_nr,:,:)-Mdata(el_nr+1,:,:),mean(Mdata(el_nr,:,:)-Mdata(el_nr+1,:,:),3)); % try bipolar
[S, f] = ieeg_getWaveletSpectrogram(squeeze(data_spect), srate, [1 175]);
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

%% %% SVD no crossval

S_temp = S_norm(:,tt>-.2 & tt<1.4,:);

S_stack = reshape(S_temp,[length(f),size(S_temp,2)*size(S_temp,3)]);

[u,s,v] = svd(S_stack,'econ');

figure,
subplot(1,3,1)
plot(f,u(:,1:3))
xlabel('frequency (Hz)')
title([all_channels.name{el_nr}])

% Reshape pc1/pc2/etc 
% Project data into first pc:
pc_nr = 1;
aa = u(:,pc_nr)'*S_stack;

subplot(1,3,2)
% plot(tt(tt>-.2 & tt<1.4),-reshape(aa,size(S_temp,2),size(S_temp,3)));
plot(tt(tt>-.2 & tt<1.4),mean(-reshape(aa,size(S_temp,2),size(S_temp,3)),2));
title(['Data projected into PC ' int2str(pc_nr)])

% Reconstruct with first pc:
aa = u(:,pc_nr)*u(:,pc_nr)'*S_stack;
bb = reshape(aa,size(S_temp,1),size(S_temp,2),size(S_temp,3));

subplot(1,3,3)
uimagesc(tt(tt>-.2 & tt<1.4),f,mean(bb,3),[-.3 .3])
axis xy
xlim([-.5 1.5])
xlabel('Time (s)'),ylabel('Frequency (Hz)')
colormap(cm); a = colorbar;
ylabel(a,'log power change wrt baseline','FontSize',12,'Rotation',90)
title(['Data reconstructed with only PC ' int2str(pc_nr)])

%%
%% SVD spectrogram + crossval
%%

% get trial labels repeats/non-repeats
eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);
[shared_idx_repeats] = shared_idx(nsd_repeats==1); % 100 images
stimuliRepeats = stimuliOrd(:,:,nsd_repeats==1,:);% get 100 images from repeats

% initialize matrix for S
S_repeats = NaN(size(S_norm,1),size(S_norm,2),6,length(shared_idx_repeats)); % frequency X time X 6 X 100

for kk = 1:length(shared_idx_repeats)
    these_trials = find(shared_idx==shared_idx_repeats(kk));    % for this repeat, find the correct 6 trial numbers out of the 1500 and get the image and the data
    S_repeats(:,:,1:length(these_trials),kk) = S_norm(:,:,these_trials); 
end

S_train = squeeze(mean(S_repeats(:,tt>-.2 & tt<1.4,1:2:end,:),3)); % frequency X time X 100
S_trainStack = reshape(S_train,[length(f),size(S_train,2)*size(S_train,3)]);

[u,s,v] = svd(S_trainStack,'econ');

% cross-validate
S_test = squeeze(mean(S_repeats(:,tt>-.2 & tt<1.4,2:2:end,:),3)); % frequency X time X 100
S_testStack = reshape(S_test,[length(f),size(S_test,2)*size(S_test,3)]);
 
nr_pcs = 10;
r_pred = zeros(nr_pcs,1);
for kk = 1:nr_pcs
    pred = u(:,1:kk)*s(1:kk,1:kk)*v(:,1:kk)';
    
    % R^2
    SSres = sum((S_testStack(:)-pred(:)).^2);
    SStot = sum((S_testStack(:)-mean(S_testStack(:))).^2);
    r_pred(kk) = 1 - SSres/SStot;
end

figure,
subplot(2,1,1)
plot(f,u(:,1:3))
subplot(2,1,2)
plot(r_pred)

% R^2 train-test
SSres = sum((S_testStack(:)-S_trainStack(:)).^2);
SStot = sum((S_testStack(:)-mean(S_testStack(:))).^2);
% 1 - SSres/SStot


%% Reshape first PC 

% Project data into first pc:
pc_nr = 5;
aa = u(:,pc_nr)'*S_stack;

figure,
subplot(1,2,1)
% plot(tt(tt>-.2 & tt<1.4),-reshape(aa,size(S_temp,2),size(S_temp,3)));
plot(tt(tt>-.2 & tt<1.4),mean(-reshape(aa,size(S_temp,2),size(S_temp,3)),2));


% Reconstruct with first pc:
aa = u(:,pc_nr)*u(:,pc_nr)'*S_stack;
bb = reshape(aa,size(S_temp,1),size(S_temp,2),size(S_temp,3));

subplot(1,2,2)
uimagesc(tt(tt>-.2 & tt<1.4),f,mean(bb,3),[-.3 .3])
axis xy
xlim([-.5 1.5])
xlabel('Time (s)'),ylabel('Frequency (Hz)')
title([all_channels.name{el_nr}])
colormap(cm); a = colorbar;
ylabel(a,'log power change wrt baseline','FontSize',12,'Rotation',90)



%% Try for multiple electrodes
%% Calculate spectograms

el_nr = find(ismember(all_channels.name,{'ROC2','ROC3','ROC4','ROC5','ROC6','ROC7','ROC8','ROC9','ROC10'}));

% stack data for first 100 stimuli
stack_data = squeeze(Mdata(el_nr,:,1:100));
stack_data = shiftdim(stack_data,1);
stack_data = reshape(stack_data,size(stack_data,1),size(stack_data,2)*size(stack_data,3));

[S, f] = ieeg_getWaveletSpectrogram(stack_data, srate, [1 175]);
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

%% SVD

S_temp = S_norm(:,tt>-.2 & tt<1.4,:);

S_stack = reshape(S_temp,[length(f),size(S_temp,2)*size(S_temp,3)]);

[u,s,v] = svd(S_stack,'econ');

figure,plot(f,u(:,1:3))

%% Reshape first PC 

pc_nr = 4;

% Project data into first pc:
aa = u(:,pc_nr)'*S_stack;

figure,
% plot(tt(tt>-.2 & tt<1.4),-reshape(aa,size(S_temp,2),size(S_temp,3)));
plot(tt(tt>-.2 & tt<1.4),mean(-reshape(aa,size(S_temp,2),size(S_temp,3)),2));


% Reconstruct with first pc:
aa = u(:,pc_nr)*u(:,pc_nr)'*S_stack;
bb = reshape(aa,size(S_temp,1),size(S_temp,2),size(S_temp,3));

figure('Position',[0 0 300 300]),
uimagesc(tt(tt>-.2 & tt<1.4),f,mean(bb,3),[-.3 .3])
axis xy
xlim([-.5 1.5])
xlabel('Time (s)'),ylabel('Frequency (Hz)')
title([all_channels.name{el_nr}])
colormap(cm); a = colorbar;
ylabel(a,'log power change wrt baseline','FontSize',12,'Rotation',90)



%% imagesc broadband and render prelim

good_channel_nrs = find(all_channels.status==1);

figure
imagesc(tt,1:length(good_channel_nrs),mean(Mbb_norm(good_channel_nrs,:,:),3),[-.2 .2])
set(gca,'YTick',1:length(good_channel_nrs),'YTickLabels',all_channels.name(good_channel_nrs))

%% render and plot noise ceiling SNR

elecsPath = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']);
elecs = ieeg_readtableRmHyphens(elecsPath);

name = all_channels.name;
all_channels_table = table(name);
elecs = ieeg_sortElectrodes(elecs, all_channels_table, 0);

% load pial and inflated giftis
gL = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['white.L.surf.gii']));
gR = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['white.R.surf.gii']));
gL_infl = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['inflated.L.surf.gii']));
gR_infl = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['inflated.R.surf.gii']));

% snap electrodes to surface and then move to inflated
xyz_inflated = ieeg_snap2inflated(elecs,gR,gL,gR_infl,gL_infl,4);

% render with Wang labeling
Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};

for hh = 1:2

    if hh==1
        hemi = 'l';
        g = gL_infl;
        views_plot = {[40,-30],[-45,-10],[-90,20]};
    elseif hh==2
        hemi = 'r';
        g = gR_infl;
        views_plot = {[-40,-30],[45,-10],[90,20]};
    end

    % surface labels 
    surface_labels = MRIread(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],'surf',...
        [hemi 'h.wang15_mplbl.mgz']));
    vert_label = surface_labels.vol(:);

    % sulcal labels
    sulcal_labels = read_curv(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],'surf',...
        [hemi 'h.sulc']));

    % load the color map
    cmap = make_WangColormap();

    electrodes_thisHemi = find(ismember(elecs.hemisphere,upper(hemi)));

    % make a plot with electrode dots
    for vv = 1:length(views_plot)
        v_d = [views_plot{vv}(1),views_plot{vv}(2)];

        % get the inflated coordinates
        els = xyz_inflated;
        % calculate popout so we can better read labels and see amplitude
        a_offset = .1*max(abs(els(:,1)))*[cosd(v_d(1)-90)*cosd(v_d(2)) sind(v_d(1)-90)*cosd(v_d(2)) sind(v_d(2))];
        els_pop = els+repmat(a_offset,size(els,1),1);

        % no labels
        figure
        tH = ieeg_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names,sulcal_labels);
        ieeg_elAdd(els(electrodes_thisHemi,:),[.99 .99 .99],25) % add electrode positions
        ieeg_elAdd(els(electrodes_thisHemi,:),[.1 .1 .1],15) % add electrode positions
        ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','render',['sub-' sub_label],...
            ['inflated_dots_sub-' sub_label '_WangAreas_v' int2str(v_d(1)) '_' int2str(v_d(2)) '_' hemi]))
         
        % with labels 
        figure
        tH = ieeg_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names,sulcal_labels);
        ieeg_elAdd(els_pop(electrodes_thisHemi,:),[.1 .1 .1],15) % add electrode positions
        ieeg_label(els_pop(electrodes_thisHemi,:),20,6,elecs.name(electrodes_thisHemi)) % add electrode names
        ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
        print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','render',['sub-' sub_label],...
            ['inflated_labels_sub-' sub_label '_WangAreas_v' int2str(v_d(1)) '_' int2str(v_d(2)) '_' hemi]))
             
%         % with activity
%         figure
%         tH = ieeg_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names,sulcal_labels);
%         all_chan_snr_plot = all_chan_snr;
%         all_chan_snr_plot(all_chan_snr_plot<.2) = 0;
%         ieeg_elAdd_sizable(els_pop(electrodes_thisHemi,:),all_chan_snr_plot(electrodes_thisHemi),.8,40) % add electrode positions
%         ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
%         set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','render',['sub-' sub_label],...
%             ['NCSNR_sub-' sub_label '_WangAreas_v' int2str(v_d(1)) '_' int2str(v_d(2)) '_' hemi]))
    end
    close all

end



