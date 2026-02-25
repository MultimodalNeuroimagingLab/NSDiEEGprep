%% This is the preprocessing script for the clean NSD analysis.
% Requires each subject to have annotated channels, events, and electrodes
% with SOZ annotated

%% Set paths, get filenames and tables for one subject
% generates data_info with information for all runs
% generates all_channels with good channels across runs marked

localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('functions');

% subject to preprocess
ss = 5;
sub_label = sprintf('%02d', ss);

ses_label = 'ieeg01';
% list of task labels that should be preprocessed and concatenated
task_labels = {'NSDspecial01','NSDspecial02','NSDspecial03','NSDspecial04','NSDspecial05',...
    'NSDspecial06','NSDspecial07','NSDspecial08','NSDspecial09','NSDspecial10'};
task_list = readtable(fullfile(localDataPath.input,['sub-' sub_label],['ses-' ses_label],['sub-' sub_label '_ses-' ses_label '_scans.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');

% electrodes path (to exclude electrodes on the SOZ)
elecsPath = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']);
elecs = ieeg_readtableRmHyphens(elecsPath);

% For all possible runs, set data_info
data_info = [];
set_counter = 0;
for ii = 1:height(task_list.filename)
    % check if task is an accepted task label
    if ismember(extractBetween(task_list.filename{ii},'task-','_run-'),task_labels)
        set_counter = set_counter+1;
        % get core of the iEEG filename
        coreName = extractBetween(task_list.filename{ii},'ieeg/','_ieeg.mefd');
        % data iEEG filename
        data_info(set_counter).filename = [coreName{1} '_ieeg.mefd'];
        % events name
        data_info(set_counter).eventsname = [coreName{1} '_events.tsv'];
        % channels name
        data_info(set_counter).channelsname = [coreName{1} '_channels.tsv'];
        % track good channels throughout all runs
        channels_table = ieeg_readtableRmHyphens(fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', data_info(set_counter).channelsname));
        if ii==1
            good_channel_bool = ismember(channels_table.type,{'SEEG','ECOG'}) & ismember(channels_table.status,'good');
        else
            good_channel_bool = ismember(channels_table.type,{'SEEG','ECOG'}) & ismember(channels_table.status,'good') & good_channel_bool==1;
        end
        data_info(set_counter).general = task_list(ii,:);
    end
end

% Channel info across runs
all_channels.name = channels_table.name;
all_channels.type = channels_table.type;
all_channels.status = good_channel_bool; % 1 is good, zero is bad, -1 is SOZ
% SOZ channels
elecsSOZ = elecs.name(contains(elecs.seizure_zone, 'SOZ')); % first set SOZ channels to bad
all_channels.soz = ismember(all_channels.name, elecsSOZ);


%% Loop through NSD01 - NSD10 and all runs individually, concatenate output at the end

Mdata = []; % data matrix
Mbb = []; % data matrix for broadband filtered data
eventsST = []; % events

for ii = 1:length(data_info)
    fprintf('Loading %s\n', data_info(ii).filename);

    %%% ------------------------------------------------------------%%% 
    %%% Load events
    try
        events =  readtable(fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ...
            data_info(ii).eventsname), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');
    catch
        warning('Error loading events file for %02d, skipping', ii); % e.g., no 2nd run for a certain NSD session
        continue
    end

    % make sure status description is a cell array to concatenate
    if ~iscell(events.status_description)
        events.status_description = cellstr(string(events.status_description));
    end

    %%% ------------------------------------------------------------%%% 
    %%% Load mef data
    [metadata, data] = readMef3(fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', data_info(ii).filename));
    srate = metadata.time_series_metadata.section_2.sampling_frequency;

    %%% ------------------------------------------------------------%%% 
    %%% Reorganize and clean up events so that adjacent rest and stim events are combined
    samprange = -2*srate:2*srate-1; % -2:+2 seconds
    tt = samprange/srate;

    % first event should be rest. Sometimes this pdiode onset is not acquired bc pdiode was on during presequence gray
    if strcmp(events.trial_type(1), 'stim')
        warning('Appending artificial rest trial to start of events');
        T1 = cell2table({events.onset(1) - 2, 2, events.sample_start(1) - 2*srate, 'rest', '', 'good', 'n/a'}, ...
                        'VariableNames', events.Properties.VariableNames);
        events = [T1; events];
    end

    % reassemble events table to be stim/test only. 
    idxST = find(ismember(events.trial_type, {'stim', 'test'}));
    eventsSTCurr = events(idxST, :); % only keep stim and test events
    % add pre (rest) status and durations
    assert(all(strcmp(events.trial_type(idxST-1), 'rest')), 'Error: trials preceding some stim/test trials are not rest');
    eventsSTCurr.pre_duration = events.duration(idxST-1, :);
    eventsSTCurr.pre_status = events.status(idxST-1, :);
    eventsSTCurr.pre_status_description = events.status_description(idxST-1, :);

    % add task (NSD1-10), run (1 2...), session (ordered blocks of time) info
    eventsSTCurr.tasknumber = repmat(data_info(ii).general.nsd_special_nr, height(eventsSTCurr), 1);
    eventsSTCurr.run = repmat(data_info(ii).general.run, height(eventsSTCurr), 1);
    eventsSTCurr.recording_set = repmat(data_info(ii).general.recording_set, height(eventsSTCurr), 1);

    % remove any events missing the full sample range (e.g. due to crash)
    eventsSTCurr(eventsSTCurr.sample_start+samprange(end)+1 > size(data, 2), :) = [];

    %%% ------------------------------------------------------------%%% 
    %%% Preprocess: highpass and CAR by lowest variance
    
    % highpass the SEEG channels
    data(strcmp(all_channels.type, 'SEEG'), :) = ieeg_highpass(data(strcmp(all_channels.type, 'SEEG'), :)', srate)';

    % re-reference: CAR by lowest variance (across data from first to last events)
    opts = struct();
    opts.vartype = 'var';
    opts.pctThresh = 25; % leave at default of 25%
    opts.winResp = [eventsSTCurr.onset(1) eventsSTCurr.onset(end)]; 
    tt_all = (1:size(data,2))/srate;
    badChs = find(all_channels.status==0 | all_channels.soz==1); % bad channels or soz
    [car_tmp, chsUsed] = ccep_CARVariance(tt_all, data, srate, badChs, [], opts);
    % only add CAR-reref data for sEEG channels
    data_car = data;
    data_car(strcmp(all_channels.type, 'SEEG'), :) = car_tmp(strcmp(all_channels.type, 'SEEG'), :);
    clear car_tmp

    %%% ------------------------------------------------------------%%% 
    %%% Epoch data
    M = zeros(size(data_car, 1), length(samprange), height(eventsSTCurr));
    for jj = 1:height(eventsSTCurr)
        M(:, :, jj) = data_car(:, (round(eventsSTCurr.onset(jj)*srate)+samprange(1)+1) : (round(eventsSTCurr.onset(jj)*srate)+samprange(end)+1)); % convert to 1 indexing
    end

    % after epoching - check the DI1 channel to make sure that the first sample of 1 occurs at t=0
    %di1 = squeeze(M(strcmp(all_channels.name, 'DigitalInput1'), :, :));
    %figure; plot(tt, di1); xlim([-0.01, 0.01]);
    
    %%% ------------------------------------------------------------%%% 
    %%% Extract broadband and epoch data
    bands = [70  80;80 90; 90 100; 100 110; 130 140; 140 150; 150 160; 160 170];  % 10 hz bins, avoiding 60 and 120
    bb_power = data_car; % initizalize with all channels
    bb_power(strcmp(all_channels.type, 'SEEG'), :) = ieeg_getHilbert(data_car(strcmp(all_channels.type, 'SEEG'), :)', bands, srate, 'power')';
    % epoch broadband
    BB = zeros(size(data_car, 1), length(samprange), height(eventsSTCurr));
    for jj = 1:height(eventsSTCurr)
        BB(:,:,jj) = bb_power(:, (round(eventsSTCurr.onset(jj)*srate)+samprange(1)+1) : (round(eventsSTCurr.onset(jj)*srate)+samprange(end)+1)); % convert to 1 indexing
    end
    
    %%% ------------------------------------------------------------%%% 
    %%% Concatenate data and events
    Mdata = cat(3, Mdata, M);
    Mbb = cat(3, Mbb, BB);
    eventsST = [eventsST; eventsSTCurr];

end

% save as single to save space
Mbb = single(Mbb);
Mdata = single(Mdata);

% Save all outputs
outdir = fullfile(localDataPath.output,'derivatives','preproc_car',['sub-' sub_label]);
if ~exist(outdir,'dir')
    sprintf(['making output directory: ' outdir])
    mkdir(outdir);
end
save(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_ieeg.mat']), 'tt', 'srate', 'Mdata', 'eventsST', 'all_channels', '-v7.3');
save(fullfile(outdir, ['sub-' sub_label '_desc-preprocCARBB_ieeg.mat']), 'tt', 'srate', 'Mbb', 'eventsST', 'all_channels', '-v7.3');
writetable(eventsST, fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_events.tsv']), 'FileType', 'text', 'Delimiter', '\t');


%% Load BB output to check

clear all
localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('functions');

subjects = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17','18','19'};

ss = 5;
sub_label = subjects{ss};
ses_label = 'ieeg01';

outdir = fullfile(localDataPath.output,'derivatives','preproc_car',['sub-' sub_label]);

% load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_ieeg.mat']), 'tt', 'srate', 'Mdata', 'eventsST', 'all_channels');
load(fullfile(outdir, ['sub-' sub_label '_desc-preprocCARBB_ieeg.mat']), 'tt', 'srate', 'Mbb', 'eventsST', 'all_channels');
readtable(fullfile(outdir, ['sub-' sub_label '_desc-preprocCAR_events.tsv']), 'FileType', 'text', 'Delimiter', '\t', 'TreatAsEmpty', 'n/a');

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

%% Find repeated images, calculate SNR

eventsST.status_description = cellstr(string(eventsST.status_description));
[events_status,nsd_idx,shared_idx,nsd_repeats] = ieeg_nsdParseEvents(eventsST);

all_chan_snr = NaN(size(Mbb_norm,1),1);
t_avg = tt>0.1 & tt<.5;

for el_nr = 1:size(Mbb_norm,1)
    
    if ismember(all_channels.type(el_nr),'SEEG') && all_channels.status(el_nr)==1
        
        bb_strength = squeeze(mean(Mbb_norm(el_nr,t_avg==1,:),2));
        
        if sum(~isnan(bb_strength))~=0

            all_repeats = find(nsd_repeats>0);
            shared_idx_repeats = unique(shared_idx(all_repeats)); % 100 images
            repeats_bb_strength = cell(length(shared_idx_repeats),1);
            for kk = 1:length(shared_idx_repeats)
                these_trials = find(shared_idx==shared_idx_repeats(kk));    % for this repeat, find the correct 6 trial numbers out of the 1500 and get the image and the data
                repeats_bb_strength{kk} = bb_strength(these_trials); 
            end
            
            % [NCSNR, p, NCSNRNull] = estimateNCSNR(repeats_bb_strength, 1000);
            [NCSNR] = estimateNCSNR(repeats_bb_strength);    % fast, no permutations

            all_chan_snr(el_nr) = NCSNR;
        else
            disp(['channel ' all_channels.name{el_nr} ' should be set to bad'])
        end
    end
end

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

for hh = 2%1:2

    if hh==1
        hemi = 'l';
        g = gL_infl;
        views_plot = {[40,-30],[-45,-10],[-90,20],[90 0]};
    elseif hh==2
        hemi = 'r';
        g = gR_infl;
        views_plot = {[-40,-30],[45,-10],[90,20],[270 0]};
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
    for vv = 1%1:length(views_plot)
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
%         print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','render',['sub-' sub_label],...
%             ['inflated_dots_sub-' sub_label '_WangAreas_v' int2str(v_d(1)) '_' int2str(v_d(2)) '_' hemi]))
         
        % with labels 
        figure
        tH = ieeg_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names,sulcal_labels);
        ieeg_elAdd(els_pop(electrodes_thisHemi,:),[.1 .1 .1],15) % add electrode positions
        ieeg_label(els_pop(electrodes_thisHemi,:),20,6,elecs.name(electrodes_thisHemi)) % add electrode names
        ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','render',['sub-' sub_label],...
%             ['inflated_labels_sub-' sub_label '_WangAreas_v' int2str(v_d(1)) '_' int2str(v_d(2)) '_' hemi]))
             
        % with activity
        figure
        tH = ieeg_RenderGiftiLabels(g,vert_label,cmap,Wang_ROI_Names,sulcal_labels);
        all_chan_snr_plot = all_chan_snr;
        all_chan_snr_plot(all_chan_snr_plot<.2) = 0;
        ieeg_elAdd_sizable(els_pop(electrodes_thisHemi,:),all_chan_snr_plot(electrodes_thisHemi),.8,40) % add electrode positions
        ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
        set(gcf,'PaperPositionMode','auto')
%         print('-dpng','-r300',fullfile(localDataPath.input,'derivatives','render',['sub-' sub_label],...
%             ['NCSNR_sub-' sub_label '_WangAreas_v' int2str(v_d(1)) '_' int2str(v_d(2)) '_' hemi]))
    end
%     close all

end



