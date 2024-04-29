%% Load CAR trial data, bipolar-rereference, and calculate broadband to save

localDataPath = personalDataPath();

ss = 1;
% for ss = 13:17
%      if ss == 4, continue; end

sub = sprintf('sub-%02d', ss);

%% 1) Load CAR data, bipolar-rereference

fprintf('Loading CAR evoked potentials for %s\n', sub);
load(fullfile(localDataPath.input, 'derivatives', 'preproc_car', sub, sprintf('%s_desc-preprocCAR_ieeg.mat', sub)));

% load segmented data to get segs for this subject
participants = readtable(fullfile(localDataPath.input, 'participants.tsv'), 'Filetype', 'text', 'Delimiter', '\t');
seg5 = participants.seg5{strcmp(participants.participant_id, sub)};
seg5 = strip(split(seg5, ','));
seg6 = participants.seg6{strcmp(participants.participant_id, sub)};
seg6 = strip(split(seg6, ','));
if strcmp(seg5, 'n/a'), seg5 = {}; end
if strcmp(seg6, 'n/a'), seg6 = {}; end

% Get SEEG channels
ephysIdxes = strcmp(all_channels.type, 'SEEG');
MdataEphys = Mdata(ephysIdxes, :, :);
nEphys = sum(ephysIdxes);
chNamesEphys = all_channels.name(ephysIdxes);

% First run twice, to get how many bipolar channels, and to get bipolar bad channel numbers, then bipolar SOZ channels
fprintf('Performing bipolar rereference per trial\n');
[~, bipolarNames, bipolarChans, badChansBip] = ieeg_bipolarSEEG(MdataEphys(:, :, 1)', chNamesEphys, find(all_channels.status(ephysIdxes) == 0), seg5, seg6);
[~, ~, ~, sozBip] = ieeg_bipolarSEEG(MdataEphys(:, :, 1)', chNamesEphys, find(all_channels.soz(ephysIdxes) == 1), seg5, seg6, false);

% Get all bipolar channels
MdataEphysBip = nan(length(bipolarNames), length(tt), height(eventsST));
for ii = 1:height(eventsST)
    MdataEphysBip(:, :, ii) = ieeg_bipolarSEEG(MdataEphys(:, :, ii)', chNamesEphys, [], seg5, seg6, false)';
end

%% 2) Calculate broadband on trial data and save
fprintf('Filtering for broadband per channel\n');
bands = [70, 80; 80, 90; 90, 100; 100, 110; 130, 140; 140, 150; 150, 160; 160, 170];  % 10 hz bins, avoiding 60 and 120, matches NSDieegPrep CAR bands
MbbBip = nan(size(MdataEphysBip));
for ii = 1:length(bipolarNames)
    fprintf('.');
    MbbBip(ii, :, :) = ieeg_getHilbert(squeeze(MdataEphysBip(ii, :, :)), bands, srate, 'power');
end
fprintf('\n');

fprintf('Saving bipolar broadband data...');
MbbBip = single(MbbBip); % single for less storage
all_channels_bipolar = struct();
all_channels_bipolar.name = bipolarNames;
all_channels_bipolar.originalChannels = bipolarChans;
all_channels_bipolar.badChannels = badChansBip;
all_channels_bipolar.soz = sozBip;
save(fullfile(localDataPath.output, 'derivatives', 'pipeline', sub, sprintf('%s_desc-preprocBipolarBB_ieeg.mat', sub)), ...
    'MbbBip', 'tt', 'srate', 'all_channels_bipolar', 'eventsST', '-v7.3');
fprintf('.\n');

%end
