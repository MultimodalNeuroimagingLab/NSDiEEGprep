%% This is the preprocessing script for the clean NSD analysis.
% Requires each subject to have annotated channels, events, and electrodes
% with SOZ annotated

%% Set paths, get filenames and tables for one subject
% generates data_info with information for all runs
% generates all_channels with good channels across runs marked

clear all
localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
addpath('functions');

subjects = {'01','02','03','04','05','06','07','08','09','10','11','12','13','14','15','16','17'};

ss = 15;
sub_label = subjects{ss};
ses_label = 'ieeg01';

%% render plain for testing before movie

elecsPath = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']);
elecs = ieeg_readtableRmHyphens(elecsPath);

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

for hh = 2

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
    for vv = 2%:length(views_plot)
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
         
    end
end


%% movie for VSS presentation

% load pial and inflated giftis
gL = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['white.L.surf.gii']));
gR = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['white.R.surf.gii']));
gL_infl = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['inflated.L.surf.gii']));
gR_infl = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['inflated.R.surf.gii']));
gL_pial = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['pial.L.surf.gii']));
gR_pial = gifti(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],['pial.R.surf.gii']));

% surface labels 
surface_labels = MRIread(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],'surf',...
    [hemi 'h.wang15_mplbl.mgz']));
vert_label = surface_labels.vol(:);

% sulcal labels
sulcal_labels = read_curv(fullfile(localDataPath.input,'sourcedata','freesurfer',['sub-' sub_label],'surf',...
    [hemi 'h.sulc']));


hemi = 'r';
if isequal(hemi,'l')
    g = gL;
    gI = gL_infl;
    gP = gL_pial;
elseif isequal(hemi,'r')
    g = gR;
    gI = gR_infl;
    gP = gR_pial;
end

% snap electrodes to surface and then move to inflated
xyz_inflated = ieeg_snap2inflated(elecs,gR,gL,gR_infl,gL_infl,4);
xyz = [elecs.x elecs.y elecs.z];


Wang_ROI_Names = {...
    'V1v' 'V1d' 'V2v' 'V2d' 'V3v' 'V3d' 'hV4' 'VO1' 'VO2' 'PHC1' 'PHC2' ...
    'TO2' 'TO1' 'LO2' 'LO1' 'V3B' 'V3A' 'IPS0' 'IPS1' 'IPS2' 'IPS3' 'IPS4' ...
    'IPS5' 'SPL1' 'FEF'};
views_plot = {[70 -10]};
vv = 1; % view 1

videoName = fullfile(localDataPath.input,'derivatives','render',['sub-' sub_label],...
        ['inflating_sub-' sub_label '_WangAreas_v' int2str(views_plot{vv}(1)) '_' int2str(views_plot{vv}(2)) '_' hemi '_movie']);

fid = figure('Color',[1 1 1]);

vidObj = VideoWriter(videoName,'MPEG-4'); %

open(vidObj); 


% load the color map
cmap = make_WangColormap();

electrodes_thisHemi = find(ismember(elecs.hemisphere,upper(hemi)));

% make the brain transparent
g_plot = gP;
v_d = [views_plot{vv}(1),views_plot{vv}(2)];

% get the coordinates
els = xyz;

tH = ieeg_RenderGiftiLabels(g_plot,vert_label,cmap,Wang_ROI_Names,sulcal_labels);
hold on, 
plot3(els(electrodes_thisHemi,1),els(electrodes_thisHemi,2),els(electrodes_thisHemi,3),'o','MarkerFaceColor', [.1 .1 .1],'MarkerEdgeColor',[.99 .99 .99],'LineWidth',1)

ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
for kk = 2:.05:9
    tH.FaceAlpha = kk/10;
    writeVideo(vidObj,getframe(fid));
end
clf
% go from pial to inflated
% gP --> gI
perc_infl = [0:.01:1];

for kk = 1:length(perc_infl)
    
    g_plot = gP;
    g_plot.vertices = (1-perc_infl(kk))*gP.vertices + perc_infl(kk)*gI.vertices;
    v_d = [views_plot{vv}(1),views_plot{vv}(2)];
    
    % get the coordinates
    els = (1-perc_infl(kk))*xyz + perc_infl(kk)*xyz_inflated;

    tH = ieeg_RenderGiftiLabels(g_plot,vert_label,cmap,Wang_ROI_Names,sulcal_labels);
    plot3(els(electrodes_thisHemi,1),els(electrodes_thisHemi,2),els(electrodes_thisHemi,3),'o','MarkerFaceColor', [.1 .1 .1],'MarkerEdgeColor',[.99 .99 .99],'LineWidth',1)
    ieeg_viewLight(v_d(1),v_d(2)) % change viewing angle   
    tH.FaceAlpha = 0.9;
    
    writeVideo(vidObj,getframe(fid));
    
    clf

end

close(vidObj);