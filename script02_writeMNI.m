%% This script calculates MNI152 and MNI305 positions each subject and saves them

%localDataPath = setLocalDataPath(1); % runs local PersonalDataPath (gitignored)
localDataPath = personalDataPath();
addpath('functions');

ss = 17;
sub_label = sprintf('%02d', ss);

ses_label = 'ieeg01';
outdir = fullfile(localDataPath.output,'derivatives','preproc_car',['sub-' sub_label]);

% load electrodes
elecsPath = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']);
elecs = ieeg_readtableRmHyphens(elecsPath);
elecmatrix = [elecs.x, elecs.y, elecs.z];

% FS dir at root level, then of subject
FSRootdir = fullfile(localDataPath.input, 'sourcedata', 'freesurfer');
FSdir = fullfile(FSRootdir, sprintf('sub-%s', sub_label));

%% 1) Get and save MNI152 positions to electrodesMni152.tsv (volumetric, SPM12)

% locate forward deformation field from SPM. There are variabilities in session name, so we use dir to find a matching one
niiPath = dir(fullfile(localDataPath.input, 'sourcedata', 'spm_forward_deformation_fields', sprintf('sub-%s_ses-*_T1w_acpc.nii', sub_label)));
assert(length(niiPath) == 1, 'Error: did not find exactly one match in sourcedata T1w MRI'); % check for only one unique match
niiPath = fullfile(niiPath.folder, niiPath.name);

% create a location in derivatives to save the transformed electrode images
rootdirMni = fullfile(localDataPath.input, 'derivatives', 'MNI152_electrode_transformations', sprintf('sub-%s', sub_label));
mkdir(rootdirMni);

% calculate MNI152 coordinates for electrodes
xyzsMni152 = ieeg_getXyzMni(elecmatrix, niiPath, rootdirMni);

% save as separate MNI 152 electrodes table
elecsMni152Path = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_space-' 'MNI152NLin2009' '_electrodes.tsv']);
elecsMni152 = elecs;
elecsMni152.x = xyzsMni152(:, 1); elecsMni152.y = xyzsBipMni152(:, 2); elecsMni152.z = xyzsBipMni152(:, 3);
writetable(elecsMni152, elecsMni152Path, 'FileType', 'text', 'Delimiter', '\t');

fprintf('Saved to %s\n', elecsMni152Path);

%% 2) Get and save MNI305 positions (through fsaverage)

% calculate MNI305 coordinates for electrodes
[xyzsMni305, vertIdxFsavg, minDists, surfUsed] = ieeg_mni305ThroughFsSphere(elecmatrix, elecs.hemisphere, FSdir, FSRootdir, 'closest', 5);

% save as separate MNI 305 electrodes table
elecsBipMni305Path = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_space-' 'MNI305' '_electrodes.tsv']);
elecsBipMni305 = elecs;
elecsBipMni305.x = xyzsMni305(:, 1); elecsBipMni305.y = xyzsMni305(:, 2); elecsBipMni305.z = xyzsMni305(:, 3);
elecsBipMni305.vertex_fsaverage = vertIdxFsavg; % also add a column to indicate vertex on fsavg, so we can easily get position for inflated brain
writetable(elecsBipMni305, elecsBipMni305Path, 'FileType', 'text', 'Delimiter', '\t');

fprintf('Saved to %s\n', elecsBipMni305Path);

%% 3) Calculate bipolar electrodes in native space, save

% note: we decide not to sort electrodes by channels before bipolar deriv, to avoid many n/a rows.
% Channels.tsv contains non-ephys channels anyway, which we would need to remove before sorting and before matching to elecs when loading

% load segmented data to get segs for this subject
participants = readtable(fullfile(localDataPath.input, 'participants.tsv'), 'Filetype', 'text', 'Delimiter', '\t');
seg5 = participants.seg5{strcmp(participants.participant_id, sprintf('sub-%s', sub_label))};
seg5 = strip(split(seg5, ','));
seg6 = participants.seg6{strcmp(participants.participant_id, sprintf('sub-%s', sub_label))};
seg6 = strip(split(seg6, ','));
if strcmp(seg5, 'n/a'), seg5 = {}; end
if strcmp(seg6, 'n/a'), seg6 = {}; end

% create dummy data to get bipolar names/channels from
Mdummy = zeros(1, height(elecs)); % note warnings here are ok, because electrodes.tsv doesn't contain all channels
[~, bipolarNames, bipolarChans] = ieeg_bipolarSEEG(Mdummy, elecs.name, [], seg5, seg6);

% get hemisphere info, checking that bipolar pair have the same hemi (or is outside of brain '')
hemiBip = elecs.hemisphere(bipolarChans);
assert(all(strcmp(hemiBip(:, 1), hemiBip(:, 2)) | any(strcmp(hemiBip, ''), 2)), ...
    'Error: hemisphere mismatch found between electrodes in bipolar pairs');
hemiBip = hemiBip(:, 1); % keep just one col

% get bipolar electrode positions
xyzsBip = nan(length(bipolarNames), 3);
for ii = 1:length(bipolarNames)
    xyzsBip(ii, :) = mean(elecmatrix(bipolarChans(ii, :), :)); % center of the pair xyzs
end

% get interpolated destrieux position (requires vistasoft
[labs, labs_val] = ieeg_getLabelXyzDestrieux(xyzsBip, FSdir, 3);

% keep seizure_zone columns for each original electrode
sozs = elecs.seizure_zone;
sozs(cellfun(@isempty, sozs)) = {'n/a'};
soz_labels = sozs(bipolarChans);

% assemble table and save
elecsBip = table(bipolarNames, xyzsBip(:, 1), xyzsBip(:, 2), xyzsBip(:, 3), hemiBip, soz_labels(:, 1), soz_labels(:, 2), labs_val, labs, ...
    'VariableNames', {'name', 'x', 'y', 'z', 'hemisphere', 'seizure_zone_1', 'seizure_zone_2', 'Destrieux_label', 'Destrieux_label_text'});
bipolarDir = fullfile(localDataPath.output, 'derivatives', 'pipeline', sprintf('sub-%s', sub_label));
mkdir(bipolarDir);
writetable(elecsBip, fullfile(bipolarDir, sprintf('sub-%s_desc-bipolar_electrodes.tsv', sub_label)), 'Filetype', 'text', 'Delimiter', '\t');

%% 4) Bipolar - Get and save MNI152 positions (run 3 first)

% locate forward deformation field from SPM. There are variabilities in session name, so we use dir to find a matching one
niiPath = dir(fullfile(localDataPath.input, 'sourcedata', 'spm_forward_deformation_fields', sprintf('sub-%s_ses-*_T1w_acpc.nii', sub_label)));
assert(length(niiPath) == 1, 'Error: did not find exactly one match in sourcedata T1w MRI'); % check for only one unique match
niiPath = fullfile(niiPath.folder, niiPath.name);

% create a location in derivatives to save the transformed electrode images. It will be hard to tell which were for bipolar, which were CAR elecs. Is this important?
rootdirMni = fullfile(localDataPath.input, 'derivatives', 'MNI152_electrode_transformations', sprintf('sub-%s', sub_label));
mkdir(rootdirMni);

% calculate MNI152 coordinates for bipolar electrodes
xyzsBipMni152 = ieeg_getXyzMni(xyzsBip, niiPath, rootdirMni);

% save as separate MNI 152 electrodes table
elecsBipMni152Path = fullfile(bipolarDir, sprintf('sub-%s_space-MNI152NLin2009_desc-bipolar_electrodes.tsv', sub_label));
elecsBipMni152 = elecsBip;
elecsBipMni152.x = xyzsBipMni152(:, 1); elecsBipMni152.y = xyzsBipMni152(:, 2); elecsBipMni152.z = xyzsBipMni152(:, 3);
writetable(elecsBipMni152, elecsBipMni152Path, 'FileType', 'text', 'Delimiter', '\t');

fprintf('Saved to %s\n', elecsBipMni152Path);

%% 5) Bipolar - Get and save MNI305 positions (run 3 first)

% calculate MNI305 coordinates for electrodes
[xyzsBipMni305, vertIdxFsavg] = ieeg_mni305ThroughFsSphere(xyzsBip, elecsBip.hemisphere, FSdir, FSRootdir, 'closest', 5);

% save as separate MNI 305 electrodes table
elecsBipMni305Path = fullfile(bipolarDir, sprintf('sub-%s_space-MNI305_desc-bipolar_electrodes.tsv', sub_label));
elecsBipMni305 = elecsBip;
elecsBipMni305.x = xyzsBipMni305(:, 1); elecsBipMni305.y = xyzsBipMni305(:, 2); elecsBipMni305.z = xyzsBipMni305(:, 3);
elecsBipMni305.vertex_fsaverage = vertIdxFsavg; % also add a column to indicate vertex on fsavg, so we can easily get position for inflated brain
writetable(elecsBipMni305, elecsBipMni305Path, 'FileType', 'text', 'Delimiter', '\t');

fprintf('Saved to %s\n', elecsBipMni305Path);

return

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

%% render and plot noise ceiling SNR

elecsPath = fullfile(localDataPath.input, ['sub-' sub_label], ['ses-' ses_label], 'ieeg', ['sub-' sub_label '_ses-' ses_label '_electrodes.tsv']);
elecs = ieeg_readtableRmHyphens(elecsPath);

name = all_channels.name;
all_channels_table = table(name);
elecs = ieeg_sortElectrodes(elecs, all_channels_table, 0);

% load pial and inflated giftis
gL = gifti(fullfile(localDataPath.input,'derivatives','freesurfer',['sub-' sub_label],['white.L.surf.gii']));
gR = gifti(fullfile(localDataPath.input,'derivatives','freesurfer',['sub-' sub_label],['white.R.surf.gii']));
gL_infl = gifti(fullfile(localDataPath.input,'derivatives','freesurfer',['sub-' sub_label],['inflated.L.surf.gii']));
gR_infl = gifti(fullfile(localDataPath.input,'derivatives','freesurfer',['sub-' sub_label],['inflated.R.surf.gii']));

% snap electrodes to surface and then move to inflated
xyz_inflated = ieeg_snap2inflated(elecs,gR,gL,gR_infl,gL_infl,4);