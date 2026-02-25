%% This is the preprocessing script for the clean NSD analysis.
% Requires each subject to have annotated channels, events, and electrodes
% with SOZ annotated

%% Set paths, get preprocessed data for one subject

localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('functions');

%% Load BB output to check

subjects = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17'};

ss = 5;
sub_label = subjects{ss};
ses_label = 'ieeg01';

outdir = fullfile(localDataPath.iEEG,'derivatives','preproc_car',['sub-' sub_label]);

load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_ieeg.mat']), 'tt', 'srate', 'Mdata', 'eventsST', 'all_channels');
load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCARBB_ieeg.mat']), 'tt', 'srate', 'Mbb', 'eventsST', 'all_channels');
readtable(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_events.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');

cm = ieeg_kjmloccolormap;

Mdata = single(Mdata);

%% Normalize bb power per run and baseline subtract signal

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


%% load all NSD images 

im_dims = [425,425,3];

stimuliOrig = zeros(im_dims(1),im_dims(2),1000,im_dims(3));
imNames = cell(1000,1);
imPaths = dir(fullfile(localDataPath.iEEG,'stimuli','shared1000','*.png'));
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
figure,plot(aa-designFile.sharedix'), title('should be all zero')
clear aa a_temp

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

%% plot single trace raw data

% sub-05
% el_nr = find(ismember(all_channels.name,'ROC4'));
% el_nr = find(ismember(all_channels.name,'RPS3'));

% sub-06 
% ROC2 nice bb and offset response
% ROC6 nice alpha - trial_nr = 39, perhaps 566 

el_nr = find(ismember(all_channels.name,'RPS3')); 
thisBB = squeeze(bb_1000(el_nr,:,:));
thisSig = squeeze(sig_1000(el_nr,:,:));

% Make a nice colormap
cm1 = [repmat([0 0 0],100,1)];
cm1(1:40,1) = [0.7]';
cm1(1:40,2) = [0.7:-0.6/39:0.1]';
cm1(1:40,3) = [0.7:-0.7/39:0]';
cm1(40:100,1) = [0.7:(1-0.7)/60:1]';
cm1(40:100,2) = [0.1:.9/60:1]';
cm2 = [repmat([0 0 0],100,1)];
cm2(1:30,3) = [0.7]';
cm2(1:30,1) = [0.7:-0.7/29:0]';
cm2(1:30,2) = [0.7:-0.7/29:0]';
cm2(30:100,3) = [0.7:(1-0.7)/70:1]';
cm2(30:100,2) = [0:1/70:1]';
cm = [cm2(end:-1:1,:); cm1];

% figure,plot(tt,mean(thisBB,2))

%%
% trial_nr = 801; % 801, 689, 664, 656 (oranges), 416 (bananas), 298 (flowers), 270 (great alpha as well), 193 not bad
% electrode ROC6, trial 39: great alpha, trial 466: nice alpha
if isequal(sub_label,'06') && isequal(all_channels.name{el_nr},'ROC6') 
    trial_nr = 39; 
elseif isequal(sub_label,'06') && isequal(all_channels.name{el_nr},'ROC2') 
    trial_nr = 664; % 270 566 801
end

figure('Position',[0 0 400 400],'Color',[1 1 1])
subplot(2,1,1),hold on
mm = max(abs(thisSig(:,trial_nr)));
fill([0 .8 .8 0],[mm mm -mm -mm],[.8 .8 .8],'EdgeColor','none','FaceAlpha',.3)
fill([-1.6 -.8 -.8 -1.6],[mm mm -mm -mm],[.8 .8 .8],'EdgeColor','none','FaceAlpha',.3)
fill([1.6 max(tt) max(tt) 1.6],[mm mm -mm -mm],[.8 .8 .8],'EdgeColor','none','FaceAlpha',.3)
plot(tt,thisSig(:,trial_nr),'k','LineWidth',1)
box off
title([all_channels.name{el_nr}])
yline(0),xline([-1.6 0 1.6])
plot([1.9 1.9],[50 150],'k','LineWidth',2)
axis off

subplot(2,1,2)
[S, f] = ieeg_getWaveletSpectrogram(thisSig(:,trial_nr), srate, [1 175]);
S = log10(S);
% normalization wrt baseline is necessary (Power/baselinePower)
S_base = mean(S(:,tt>-0.5 & tt<0),2);
S_norm = S - repmat(S_base,1,size(S,2));
uimagesc(tt,f,mean(S_norm,3),[-1.5 1.5]);

axis xy
xline([-1.6 0 1.6])
xlabel('Time (s)'),ylabel('Frequency (Hz)')
colormap(cm); %a = colorbar;

% % plot(tt,smooth(mean(S_norm(f>10 & f<20,:),1),srate/10))
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','traces',['sub-' sub_label],...
    [sub_label '_' all_channels.name{el_nr} '_tr' int2str(trial_nr)]))
print('-depsc','-r300','-painters',fullfile(localDataPath.input,'derivatives','traces',['sub-' sub_label],...
    [sub_label '_' all_channels.name{el_nr} '_tr' int2str(trial_nr)]))

%% plot images shown
figure('Position',[0 0 700 300])
bb = squeeze(stimuliOrig(:,:,trial_nr,:));
aa = squeeze(stimuliOrig(:,:,shared_idx(find(shared_idx==trial_nr)-1,:),:));
cc = squeeze(stimuliOrig(:,:,shared_idx(find(shared_idx==trial_nr)+1,:),:));
subplot(1,3,1)
imshow(uint8(aa))
subplot(1,3,2)
imshow(uint8(bb))
subplot(1,3,3)
imshow(uint8(cc))
set(gcf,'PaperPositionMode','auto')
print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','traces',['sub-' sub_label],...
    [sub_label '_' all_channels.name{el_nr} '_tr' int2str(trial_nr) '_images']))


%% rank broadband over time bins

% sub-05: RPS 2 - body
% sub-05: RPS 3 - food

% el_nr = find(ismember(all_channels.name,'ROC4'));
el_nr = find(ismember(all_channels.name,'RPS4'));
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


%% How similar is this representation to the previous one?

% average bb per t_int window
thisBB_av1000 = NaN(size(t_int,1),1000);
for kk = 1:size(t_int,1)
    thisBB_av1000(kk,:) = mean(thisBB(tt>t_int(kk,1) & tt<t_int(kk,2),:),1);
end
R = faster_corr_mtrx([thisBB_av1000']);

%% try svd
[u,s,v] = svd(thisBB(tt>-0.1 & tt<1.05,:),'econ');
figure,
subplot(2,2,1)
plot(cumsum(diag(s).^2)./sum(diag(s).^2))
subplot(2,2,2),
plot(tt(tt>-0.1 & tt<1.05),u(:,1:3))
subplot(2,2,3),
imagesc(thisBB(tt>-0.1 & tt<1.05,:)',[-0.5 0.5])
colormap(cm)


%% representational similarity over time?
%% only show for repeated images?
%%
%% Calculate SNR over time for special 100
%%

% get 100 images from repeats
stimuliRepeats = stimuliOrig(:,:,special100_idx>0,:);
    
el_nr = find(ismember(all_channels.name,'RPS3'));
thisBB = squeeze(bb_100(el_nr,:,:,:));

t_int = [[-.1:.05:1]' [-.05:.05:1.05]'];
thisBB_av = NaN(size(t_int,1),100,6);

this_chan_snr = NaN(size(t_int,1),1);

repeats_bb_strength = cell(size(t_int,1),1);

for kk = 1:size(t_int,1)
    repeats_bb_strength = cell(100,1);
    for ii = 1:100
        thisBB_av(kk,ii,:) = mean(thisBB(tt>t_int(kk,1) & tt<t_int(kk,2),ii,:),1);
        repeats_bb_strength{ii} = squeeze(mean(thisBB(tt>t_int(kk,1) & tt<t_int(kk,2),ii,:),1));
    end
    [NCSNR, p, NCSNRNull] = estimateNCSNR(repeats_bb_strength, 1000);    

    this_chan_snr(kk) = NCSNR;
end
figure,
subplot(3,1,1) 
plot(t_int(:,1),mean(thisBB_av,[2 3]))
subplot(3,1,2) 
plot(t_int(:,1),this_chan_snr)

subplot(3,1,3) 
r = corr(mean(thisBB_av(:,:,1:2:6),3)',mean(thisBB_av(:,:,2:2:6),3)');

% R = faster_corr_mtrx([thisBB_av1000']);

%% stimuliRepeats
 
av_bb_repeats = mean(thisBB_av,3);

figure
subplot(1,2,1)
[bb_sort,ii_sort] = sort(av_bb_repeats(7,:),'descend');
f = makeimagestack_wrapper(stimuliRepeats(:,:,ii_sort(1:50),:),0,1,[5 10]);
imshow(uint8(f))

subplot(1,2,2)
[bb_sort,ii_sort] = sort(av_bb_repeats(7,:),'ascend');
f = makeimagestack_wrapper(stimuliRepeats(:,:,ii_sort(1:50),:),0,1,[5 10]);
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





%% Calculate spectograms

el_nr = find(ismember(all_channels.name,'ROC4'));

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

%% SVD

S_temp = S_norm(:,tt>-.2 & tt<1.4,:);

S_stack = reshape(S_temp,[length(f),size(S_temp,2)*size(S_temp,3)]);

[u,s,v] = svd(S_stack,'econ');

figure,plot(f,u(:,1:3))

%% Reshape first PC 

% Project data into first pc:
aa = u(:,1)'*S_stack;

figure,
% plot(tt(tt>-.2 & tt<1.4),-reshape(aa,size(S_temp,2),size(S_temp,3)));
plot(tt(tt>-.2 & tt<1.4),mean(-reshape(aa,size(S_temp,2),size(S_temp,3)),2));


% Reconstruct with first pc:
aa = u(:,1)*u(:,1)'*S_stack;
bb = reshape(aa,size(S_temp,1),size(S_temp,2),size(S_temp,3));

figure('Position',[0 0 300 300]),
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



