%% This is the preprocessing script for the clean NSD analysis.
% Requires each subject to have annotated channels, events, and electrodes
% with SOZ annotated

%% Set paths, get preprocessed data for one subject

localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('testing');

%% Load BB output to check

subjects = {'02','04','06','07','17','18'};
sub_elecs = {{'LOC1','LOC2','LOC3','LOC4','LOC5'},... % 1:sub-02
    {'LPI1','LPI2','LPI2','LPI4'},... %2: sub-04 - nystagmus
    {'ROC1','ROC2','ROC3','ROC4','ROC5'},...%3: sub-06
    {'ROC1','ROC2','ROC3','ROC4'},... %4: sub-07
    {'LOC1','LOC2','LOC3'},...%5: sub-17
    {'ROc1','ROc2'}}; %6: sub-18

% todo: sub-18

for ss = 6%1:5
    sub_label = subjects{ss};
    ses_label = 'ieeg01';
    
    outdir = fullfile(localDataPath.output,'derivatives','preproc_car',['sub-' sub_label]);
    
    load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCARBB_ieeg.mat']), 'tt', 'srate', 'Mbb', 'eventsST', 'all_channels');
    readtable(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_events.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');
    
    cm = ieeg_kjmloccolormap;
    
    
    % Normalize bb power and signal per run
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
    
    % load all NSD images and reorder data according to images
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
    
    eventsST.status_description = cellstr(string(eventsST.status_description));
    [events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);
    
    bb_1000 = NaN(size(Mbb_norm,1),size(Mbb_norm,2),1000);
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
            bb_1000(:,:,kk) = mean(Mbb_norm(:,:,these_trials),3,'omitnan');
        end
        if length(these_trials)>1
            count_100 = count_100+1;
            special100_idx(kk) = count_100;
            bb_100(:,:,count_100,1:length(these_trials)) = Mbb_norm(:,:,these_trials);
        end
    end
    
    % save bb1000 and bb100 for relevant electrodes
    for ii_el = 1:length(sub_elecs{ss})
        el_name = sub_elecs{ss}{ii_el};
        el_nr = find(ismember(all_channels.name,el_name));
        bb_el_1000 = squeeze(bb_1000(el_nr,:,:));
        bb_el_100 = squeeze(bb_100(el_nr,:,:,:));
        out_file = fullfile(localDataPath.output,'derivatives','zeeshan',['sub-' sub_label],...
            ['sub-' sub_label '_el' sub_elecs{ss}{ii_el} '_bb1000']);
    
        save(out_file,'special100_idx','tt','srate','bb_el_1000','bb_el_100')
    end
end

%%
%% now we can load small files for all subjects
%%

localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)

subjects = {'02','06','07','17','18'};
% sub_elecs = {{'LOC2','LOC3','LOC4','LOC5'},... % sub-02
%     {'LPI1','LPI2','LPI2','LPI4'},... % sub-04 - nystagmus
%     {'ROC1','ROC2','ROC3','ROC4','ROC5'},...% sub-06
%     {'ROC1','ROC2','ROC3','ROC4'},... % sub-07
%     {'LOC1','LOC2','LOC3'},... % sub-17
%     {'ROC1','ROC2'}};% sub-18
sub_elecs = {{'LOC1'},... % sub-02 
    {'ROC2','ROC3','ROC4','ROC5'},...% sub-06
    {'ROC1','ROC2','ROC3'},... % sub-07
    {'LOC1','LOC3'},...% sub-17
    {'ROc1','ROc2'}};% sub-18

% load all NSD images
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

% Load fmri data
% fMRI shared NSD (average across all subjects)
nsd_fmri = load(fullfile(localDataPath.fmri,'betaAvg_shared_allsub.mat'));
designFile = load(fullfile(localDataPath.fmri,'nsd_expdesign.mat'));
% make sure order is the same as iEEG
aa = zeros(1000,1);
for kk = 1:1000
    a_temp = extractBetween(imNames{kk},'_nsd','.png');
    aa(kk) = str2num(a_temp{1});
end

% load atlases 
[v_Rosenke_r, v_Rosenke_label, ~] = loadAtlas( localDataPath.fs, 'r', 'visf');
[v_BensonV_r, v_BensonV_label, ~] = loadAtlas( localDataPath.fs, 'r', 'BensonV');
[v_Wang_r, v_Wang_label, ~] = loadAtlas( localDataPath.fs, 'r', 'Wang');
[v_HCP_r, v_HCP_label, ~] = loadAtlas( localDataPath.fs, 'r', 'HCP');
% [v_BensonE, ~, ~] = loadAtlas( localDataPath.fs, hemi, 'BensonE');

[v_Rosenke_l, ~, ~] = loadAtlas( localDataPath.fs, 'l', 'visf');
[v_BensonV_l, ~, ~] = loadAtlas( localDataPath.fs, 'l', 'BensonV');
[v_Wang_l, ~, ~] = loadAtlas( localDataPath.fs, 'l', 'Wang');
[v_HCP_l, ~, ~] = loadAtlas( localDataPath.fs, 'l', 'HCP');

%% now we can run through each V1-3 electrode

out = [];

for ss = 1:5
    disp(['sub ' int2str(ss)])
    sub_label = subjects{ss};
    ses_label = 'ieeg01';
    for ii_el = 1:length(sub_elecs{ss})
        el_name = sub_elecs{ss}{ii_el};
        bb1000_filename = fullfile(localDataPath.output,'derivatives','zeeshan',['sub-' sub_label],...
            ['sub-' sub_label '_el' sub_elecs{ss}{ii_el} '_bb1000']);
        load(bb1000_filename,'special100_idx','tt','srate','bb_el_1000','bb_el_100')
        hemi = lower(el_name(1));
        
        if isequal(hemi,'r')
            fmri_beta = nsd_fmri.beta_avgR; 
            v_Rosenke = v_Rosenke_r;
            v_BensonV = v_BensonV_r;
            v_Wang = v_Wang_r;
            v_HCP = v_HCP_r;
        elseif isequal(hemi,'l')
            fmri_beta = nsd_fmri.beta_avgL; 
            v_Rosenke = v_Rosenke_l;
            v_BensonV = v_BensonV_l;
            v_Wang = v_Wang_l;
            v_HCP = v_HCP_l;        
        end

        clear vertmask
        vertmask(1).idx = v_BensonV==1;
        vertmask(1).nr = length(find(vertmask(1).idx==1));
        vertmask(2).idx = v_BensonV==2;
        vertmask(2).nr = length(find(vertmask(2).idx==1));
        vertmask(3).idx = v_BensonV==3;
        vertmask(3).nr = length(find(vertmask(3).idx==1));
        vertmask(4).idx = v_BensonV==4;
        vertmask(4).nr = length(find(vertmask(4).idx==1));
%         vertmask(5).idx = v_Rosenke==1 | v_Rosenke==2; % faces
        % % vertmask(5).idx = v_Rosenke==4 | v_Rosenke==5 | v_Rosenke==6 | v_Rosenke==7; % faces
        % % vertmask(5).idx = v_Rosenke>0 & v_Rosenke<10; % VTC
%         vertmask(5).idx = v_Rosenke==8 | v_Rosenke==9; % places
        % % vertmask(5).idx = (v_HCP>120 & v_HCP<139);
        % % vertmask(5).idx = (v_Wang==16 | v_Wang==17); %Wang V3B
%         vertmask(5).nr = length(find(vertmask(5).idx==1));
%         vertmask(6).idx = (v_Wang==18 | v_Wang==19 | v_Wang==20 | v_Wang==21); %Wang IPS
%         vertmask(6).nr = length(find(vertmask(6).idx==1));

        % used for VSS:
        vertmask(5).idx = (v_BensonV==5 | v_BensonV==6); %V01 V02
        vertmask(5).nr = length(find(vertmask(5).idx==1));
        vertmask(6).idx = (v_Wang==10 | v_Wang==11); %PHC1 PHC2
        vertmask(6).nr = length(find(vertmask(6).idx==1));
        
        %%%% find vertices with shared computation
        t_ints = [.1 .7];
        bb = mean(bb_el_1000(tt>t_ints(1) & tt<t_ints(2),special100_idx==0),1)';
        vert_shared_comp = [];
        for jj = 1:length(vertmask)
            fmri_beta_mask = fmri_beta(vertmask(jj).idx==1,special100_idx==0); 
            corr_out = zeros(1,size(fmri_beta_mask,1));
            for kk = 1:size(fmri_beta_mask,1)
                if mod(kk,10000)==0, disp([int2str(kk) ' of ' int2str(size(fmri_beta_mask,1))]), end
                R = faster_corr_mtrx([bb(~isnan(bb)) fmri_beta_mask(kk,~isnan(bb))']);
                corr_out(kk) = R(1,2);
            end
            vert_shared_comp(jj).corr = corr_out';
        end       

        % test when shared computation happens
        t_ints = [[-.420:0.020:1.000]' [-.400:0.020:1.020]'];
        sum_stats = zeros(4,size(t_ints,1),length(vertmask));
        
        for jj = 1:length(vertmask)
            fmri_beta_mask = fmri_beta(vertmask(jj).idx==1,special100_idx>0); 
        
            % betas from vertices with shared representation
        %     beta_test = fmri_beta_mask(vert_shared_comp(jj).corr>r_th,:); 
            beta_test = fmri_beta_mask(vert_shared_comp(jj).corr>quantile(vert_shared_comp(jj).corr,.95),:); 
        
            corr_out = zeros(size(t_ints,1),size(beta_test,1));
            for ii = 1:size(t_ints,1)
                % disp(['starting ' int2str(ii) ' of ' int2str(size(t_ints,1))])
        
                bb_test = mean(bb_el_1000(tt>t_ints(ii,1) & tt<t_ints(ii,2),special100_idx>0),1)';
            
                for kk = 1:size(beta_test,1)
                    if mod(kk,10000)==0, disp([int2str(kk) ' of ' int2str(size(beta_test,1))]), end
                    R = faster_corr_mtrx([bb_test(~isnan(bb_test)) beta_test(kk,~isnan(bb_test))']);
                    corr_out(ii,kk) = R(1,2);
                end
            end
            
            out(ss).el(ii_el).area(jj).corr = corr_out; 
        end       

    end
end


%%

% initialize
for ii_area = 1:length(vertmask)
    r_stack(ii_area).corr = [];
end
r_subInd = [];

for ss = 1:length(out)
    for ii_el = 1:length(out(ss).el)
        r_subInd = [r_subInd; ss];
        for ii_area = 1:length(vertmask)
            r_stack(ii_area).corr = [r_stack(ii_area).corr mean(out(ss).el(ii_el).area(ii_area).corr,2)];
        end
    end
end

% cm = lines(length(vertmask));
cm = brewermap(length(vertmask),'Set1');
cm = parula(length(vertmask)+1).^2;

figure('Position',[0 0 500 300])
subplot(1,2,1), hold on
for ii_area = [1 3:6]%1:length(vertmask)
    % plot(t_ints(:,1),mean(r_stack(ii_area).corr,2),'Color',cm(ii_area,:))
    ieeg_plotCurvConf(t_ints(:,1),r_stack(ii_area).corr',cm(ii_area,:),.3,10000)
    plot(t_ints(:,1),mean(r_stack(ii_area).corr,2),'Color',cm(ii_area,:),'LineWidth',2)
    
end
yline(0)
xlim([-0.2 t_ints(end,1)]),ylim([-0.1 .8])
ylabel('Pearsons r'),xlabel('Time(s)')
title('iEEG-fMRI similarity')

subplot(1,2,2), hold on
conf_int = [];
for ii_area = [1 3:6]
    r_plot = r_stack(ii_area).corr./sqrt(sum(r_stack(ii_area).corr.^2,1));
    plot(t_ints(:,1),mean(r_plot,2),'Color',cm(ii_area,:))
%     plot(t_ints(:,1),r_plot,'Color',cm(ii_area,:))
    [h,conf] = ieeg_plotCurvConf(t_ints(:,1),r_plot',cm(ii_area,:),.3,10000);
    conf_int(ii_area,:) = conf;
    plot(t_ints(:,1),mean(r_plot,2),'Color',cm(ii_area,:),'LineWidth',2)
    
end
yline(0)
xlim([-0.2 t_ints(end,1)])
ylim([-0.05 .3])
ylabel('Unit length normalized'),xlabel('Time(s)')
title('Normalized iEEG-fMRI similarity')
set(gca,'YTick',[0:.1:.3])


% V1>V3
v1_lowerbound = conf_int(1,1:size(t_ints,1));
v2_upperbound = conf_int(2,end:-1:size(t_ints,1)+1);
v3_upperbound = conf_int(3,end:-1:size(t_ints,1)+1);
v4_upperbound = conf_int(4,end:-1:size(t_ints,1)+1);
VO_upperbound = conf_int(5,end:-1:size(t_ints,1)+1);
PHC_upperbound = conf_int(6,end:-1:size(t_ints,1)+1);
if ~isempty(find(v1_lowerbound>v3_upperbound))
    plot(t_ints(find(v1_lowerbound>v3_upperbound),1),.25,'.','Color',cm(3,:))
end
if ~isempty(find(v1_lowerbound>v4_upperbound))
    plot(t_ints(find(v1_lowerbound>v4_upperbound),1),.255,'.','Color',cm(4,:))
end
if ~isempty(find(v1_lowerbound>VO_upperbound))
    plot(t_ints(find(v1_lowerbound>VO_upperbound),1),.26,'.','Color',cm(5,:))
end
if ~isempty(find(v1_lowerbound>PHC_upperbound))
    plot(t_ints(find(v1_lowerbound>PHC_upperbound),1),.265,'.','Color',cm(6,:))
end
set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','zeeshan',...
%         ['foursub_similarity']))
% print('-depsc','-painters','-r300',fullfile(localDataPath.input,'derivatives','zeeshan',...
%         ['foursub_similarity']))


%% per electrode effect 

figure
for ii_area = 1:length(vertmask)
    for kk = 1:size(r_stack(ii_area).corr,2)
        subplot(5,3,kk), hold on
        r_plot = r_stack(ii_area).corr(:,kk)./sqrt(sum(r_stack(ii_area).corr(:,kk).^2,1));
        plot(t_ints(:,1),r_plot,'Color',cm(ii_area,:))
    end
end

%% show for one electrode in nice figure

figure('Position',[0 0 500 300])
subplot(1,2,1), hold on
for ii_area = [1 3:6]
    for kk = 4
        plot(t_ints(:,1),r_stack(ii_area).corr(:,kk),'Color',cm(ii_area,:))
    end
end
yline(0)
xlim([-0.2 t_ints(end,1)]),ylim([-0.1 .8])
ylabel('Pearsons r'),xlabel('Time(s)')
title('iEEG-fMRI similarity')
legend({'V1','V3','V4','VO','PHC'})


subplot(1,2,2), hold on
for ii_area = [1 3:6]
    for kk = 4
        r_plot = r_stack(ii_area).corr(:,kk)./sqrt(sum(r_stack(ii_area).corr(:,kk).^2,1));
        plot(t_ints(:,1),r_plot,'Color',cm(ii_area,:))
    end
end
yline(0)
xlim([-0.2 t_ints(end,1)])
ylim([-0.05 .3])
ylabel('Unit length normalized'),xlabel('Time(s)')
title('Normalized iEEG-fMRI similarity')
set(gca,'YTick',[0:.1:.3])
set(gcf,'PaperPositionMode','auto')
% print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','zeeshan',...
%         ['single_electrode_similarity']))
% print('-depsc','-painters','-r300',fullfile(localDataPath.input,'derivatives','zeeshan',...
%         ['single_electrode_similarity']))


%% overall correlation with V1 check - this should be good

for ii_area = 1%:length(vertmask)
    figure
    for ii_el = 1:size(r_stack(ii_area).corr,2)
        subplot(4,5,ii_el)
        plot(r_stack(ii_area).corr(:,ii_el))
        ylim([-0.2 1])

    end
end


