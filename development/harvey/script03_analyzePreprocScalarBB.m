%% Script to perform first-step ERP, broadband analysis on preprocessed Mef data
% First, load preprocessed data and all giftis

% cd here
setPathsNSD;
addpath('functions');
cmkjm = ieeg_kjmloccolormap();

sub = 'X';

% Determine which task/runs are bad to ignore, if any
fprintf('analyzePreprocScalarBB01.m on sub-%s\n', sub);
[~, ~, badtasks] = subjectConfig(sub);

% Load electrodes
elecPath = fullfile(localRoot, sprintf('sub-%s', sub), 'ses-ieeg01', 'ieeg', sprintf('sub-%s_ses-ieeg01_electrodes.tsv', sub));
elecs = ieeg_readtableRmHyphens(elecPath);

% load preprocessed mef data from local copy
outdir = fullfile('output', 'fromMef', sub);
tic; load(fullfile(localRoot, sprintf('sub-%s', sub), 'NSDpreproc', sprintf('preprocessed_%s.mat', sub))); toc;

% sort electrodes
elecs = ieeg_sortElectrodes(elecs, channelsEphys);

% extract hemisphere information from electrodes
hemis = unique(cellfun(@(s) {s(1)}, lower(elecs.name)));
disp('Hemispheres:');
disp(hemis);

% Remove bad trials from events and taskruns table based on task/runs to ignore
badTrs = zeros(height(eventsST), 1, 'logical');
for ii = 1:size(badtasks, 1)
    badTrs(eventsST.tasknumber == badtasks(ii, 1) & eventsST.run == badtasks(ii, 2)) = true;
    taskruns(taskruns.tasknumber == badtasks(ii, 1) & taskruns.run == badtasks(ii, 2), :) = [];
end
fprintf('Removing %d trials from bad tasks\n', sum(badTrs));
eventsST(badTrs, :) = [];
Mephys(:, :, badTrs) = [];
Meye(:, :, badTrs) = [];

% Make relevant output directories
mkdir(fullfile(outdir, 'ERP'));
mkdir(fullfile(outdir, 'broadband'));
mkdir(fullfile(outdir, 'giftiLabeled'));

% Load giftis (pial, inflated, and white (used to transform to inflated))
giis = cellfun(@(hemi) gifti(fullfile(localRoot, sprintf('sub-%s', sub), sprintf('pial.%s.surf.gii', hemi))), upper(hemis), 'UniformOutput', false);
giiInfs = cellfun(@(hemi) gifti(fullfile(localRoot, sprintf('sub-%s', sub), sprintf('inflated.%s.surf.gii', hemi))), upper(hemis), 'UniformOutput', false);
whites = cellfun(@(hemi) gifti(fullfile(localRoot, sprintf('sub-%s', sub), sprintf('white.%s.surf.gii', hemi))), upper(hemis), 'UniformOutput', false);
sulcs = cellfun(@(hemi) read_curv(fullfile(localRoot, sprintf('sub-%s', sub), sprintf('%sh.sulc', hemi))), lower(hemis), 'UniformOutput', false);
for ii = 1:length(hemis)
    assert(length(giis{ii}.faces) == length(giiInfs{ii}.faces) && length(giis{ii}.faces) == length(whites{ii}.faces), 'Error: face mismatch between %dth gifti/giiInf/white', ii);
end

% Create bipolar-rereferenced data, for broadband analysis. Load manually-written lead-segment metadata
segmentedLeads = readtable(fullfile('data', 'metadata', sprintf('%s_segmentedLeads.txt', sub)), 'FileType', 'text', 'Delimiter', '\t');
seg5 = {}; seg6 = {};
if ~isnumeric(segmentedLeads.seg5), seg5 = strtrim(split(upper(segmentedLeads.seg5), {' ', ','})); seg5(cellfun(@isempty, seg5)) = []; end
if ~isnumeric(segmentedLeads.seg6), seg6 = strtrim(split(upper(segmentedLeads.seg6), {' ', ','})); seg6(cellfun(@isempty, seg6)) = []; end
[~, bipolarNames, bipolarChans] = ieeg_bipolarSEEG(Mephys(:, :, 1)', channelsEphys.name, [], seg5, seg6); % initialize to get number of bipolar channels
MephysBip = nan(length(bipolarNames), length(tt), height(eventsST));
for ii = 1:height(eventsST)
    MephysBip(:, :, ii) = ieeg_bipolarSEEG(Mephys(:, :, ii)', channelsEphys.name, [], seg5, seg6, false)'; % don't consider bad channels here
end

% Determine run-specific bipolar bad channels
for ii = 1:height(taskruns)
    [~, ~, ~, badChs] = ieeg_bipolarSEEG(Mephys(:, :, 1)', channelsEphys.name, taskruns.badChs{ii}, seg5, seg6, false);
    taskruns.badChsBip{ii} = badChs;
end

% get interpolated bipolar pair coordinates
elecsBip = ieeg_Els2BipLocs(elecs, bipolarNames);

% Make masks indicating for each channel with trials are bad (i.e., which task-runs are bad). % can add to this later other bad trials
badTrsMask = zeros(height(channelsEphys), height(eventsST), 'logical');
badTrsMaskBip = zeros(length(bipolarNames), height(eventsST), 'logical');
for ii = 1:height(taskruns)
    badTrsMask(taskruns.badChs{ii}, eventsST.tasknumber == taskruns.tasknumber(ii) & eventsST.run == taskruns.run(ii)) = true;
    badTrsMaskBip(taskruns.badChsBip{ii}, eventsST.tasknumber == taskruns.tasknumber(ii) & eventsST.run == taskruns.run(ii)) = true;
end

%% Calculate inflated pial positions for elecs and bipolar-interpolated elecs, add to elecs and elecsBip tables

fprintf('Plotting labeled pial and inflated pials\n');

% Calculate inflated bipolar electrode positions
elecs.xInf = nan(height(elecs), 1);
elecs.yInf = nan(height(elecs), 1);
elecs.zInf = nan(height(elecs), 1);
elecsBip.xInf = nan(height(elecsBip), 1);
elecsBip.yInf = nan(height(elecsBip), 1);
elecsBip.zInf = nan(height(elecsBip), 1);

% check each hemisphere separately
for ii = 1:length(hemis)

    % Calculate inflated positions for real electrodes
    xyzInf = getXyzsInf(elecs, hemis{ii}, whites{ii}, giiInfs{ii}, 4, false); % must be within 4 mm of Gray/white boundary

    % also check 2 mm of pial surface, add any additional discovered in case of thick GM
    xyzInfPial = getXyzsInf(elecs, hemis{ii}, giis{ii}, giiInfs{ii}, 2, false);
    extraElecsPial = isnan(xyzInf(:, 1)) & ~isnan(xyzInfPial(:, 1));
    xyzInf(extraElecsPial, :) = xyzInfPial(extraElecsPial, :);

    idxes = ~isnan(xyzInf(:, 1)); % inflated positions found for current hemisphere
    elecs.xInf(idxes) = xyzInf(idxes, 1); elecs.yInf(idxes) = xyzInf(idxes, 2); elecs.zInf(idxes) = xyzInf(idxes, 3);

    % Calculate inflated positions for bipolar-interpolated electrode pairs
    xyzInf = getXyzsInf(elecsBip, hemis{ii}, whites{ii}, giiInfs{ii}, 4, false); % must be within 4 mm of Gray/white boundary

    % also check 2 mm of pial surface, add any additional discovered in case of thick GM
    xyzInfPial = getXyzsInf(elecsBip, hemis{ii}, giis{ii}, giiInfs{ii}, 2, false);
    extraElecsPial = isnan(xyzInf(:, 1)) & ~isnan(xyzInfPial(:, 1));
    xyzInf(extraElecsPial, :) = xyzInfPial(extraElecsPial, :);

    idxes = ~isnan(xyzInf(:, 1)); % inflated positions found for current hemisphere
    elecsBip.xInf(idxes) = xyzInf(idxes, 1); elecsBip.yInf(idxes) = xyzInf(idxes, 2); elecsBip.zInf(idxes) = xyzInf(idxes, 3);
end

% find which channels are good in at least one run
goodChsOn1 = goodChsNRuns(taskruns.badChs, height(elecs), 1);
goodChsBipOn1 = goodChsNRuns(taskruns.badChsBip, height(elecsBip), 1);

% Plot electrode positions to gifti for reference (good elecs only)
plotElectrodesToGifti(elecs(goodChsOn1, :), elecsBip(goodChsBipOn1, :), giis, giiInfs, sulcs, hemis, fullfile(outdir, 'giftiLabeled'), 'off');

%% Calculate PSDs and average (scalar) broadband power for CAR electrodes

fprintf('CAR electrodes\n');

% include channels that are good on 5+ runs
[goodChsOn5, badChsOn5] = goodChsNRuns(taskruns.badChs, height(elecs), 5);

% Calculate PSD for each trial
intervalStim = [0, 0.5]; intervalRest = [-0.5, 0]; % stim = first 500ms after stimulus, rest = last 500 ms before stimulus
pwelchWin = hann(srate/4);
psdStim = nan(height(channelsEphys), srate/2+1, height(eventsST));
psdRest = nan(height(channelsEphys), srate/2+1, height(eventsST));
for ii = 1:height(channelsEphys)
    [psdStim(ii, :, :), f] = pwelch(squeeze(Mephys(ii, tt >= intervalStim(1) & tt < intervalStim(2), :)), pwelchWin, length(pwelchWin)/2, srate, srate);
    psdRest(ii, :, :) = pwelch(squeeze(Mephys(ii, tt >= intervalRest(1) & tt < intervalRest(2), :)), pwelchWin, length(pwelchWin)/2, srate, srate);
end

% Plot average stationary PSD for all channels and save
toplot = nan(height(channelsEphys), length(f));
for ii = 1:height(channelsEphys)
    toplot(ii, :) = mean(log10(psdStim(ii, :, ~badTrsMask(ii, :))), 3, 'omitnan') - mean(log10(psdRest(ii, :, ~badTrsMask(ii, :))), 3, 'omitnan');
end
toplot(badChsOn5, :) = nan;
figure('Position', [200, 200, 600, 1200], 'Visible', 'off');
imagescNaN(f, 1:height(channelsEphys), toplot, 'CLim', [-0.2, 0.2]); colormap(cmkjm);
xlim([0, 200]);
set(gca, 'YTick', 1:2:height(channelsEphys), 'YTickLabels', channelsEphys.name(1:2:end));
xlabel('Frequency (Hz)');
ylabel('channels');
colorbar;
saveas(gcf, fullfile(outdir, sprintf('PSDsummary_%s.png', sub)));

% Calculate average broadband power for each trial (geomean across frequency bins as in Yuasa et al., 2023)
BBpowerStim = squeeze(geomean(psdStim(:, f >= 70 & f <= 115 | f >= 125 & f <= 170, :), 2));
BBpowerRest = squeeze(geomean(psdRest(:, f >= 70 & f <= 115 | f >= 125 & f <= 170, :), 2));

%% Plot Rsq of BB power stim vs rest as activations to brain

% Calculate rsq for each channel, across trials that are good for each channel
rsqBB = nan(height(channelsEphys), 2);
for ii = 1:height(channelsEphys)
    [rsq, p] = mnl_rsq(BBpowerStim(ii, ~badTrsMask(ii, :)), BBpowerRest(ii, ~badTrsMask(ii, :)));
    rsqBB(ii, :) = [rsq, p];
end
% set bad channels (based on N runs) to nan
rsqBB(badChsOn5, :) = nan;

% Add to bipolar electrodes table and write as output
elecs.BBrsq = rsqBB(:, 1);
elecs.BBrsqP = rsqBB(:, 2);
writetable(elecs, fullfile(outdir, 'elecsBB1.tsv'), 'FileType', 'text', 'Delimiter', '\t');

% Plot broadband Rsq activations on gifti and inflated gifti
plotActivationsToGiftis([], elecs, elecs.BBrsq, giis, giiInfs, sulcs, hemis, 0.5, fullfile(outdir, 'broadband'), 'BBrsq', 'off');

% save a colorbar
fCbar = plotCtxWtsAdjCbar(0.5, 'SigFraction', 0.05, 'SizeJump', 2);
saveas(fCbar, fullfile(outdir, 'broadband', 'cbar_BBrsq.png'));

%% Calculate noise ceiling SNR for each CAR electrode using the special100 images

% find indices of the special100 images
special100Idxes = getSpecial100Idxes(eventsST.stim_file);

% copy BBpowerStim, and put nans in bad chxtrs. These are removed in estimateNCSNR function
BBpowerStimNan = BBpowerStim;
BBpowerStimNan(badTrsMask) = nan;

rng('default');

NCSNRs = nan(height(channelsEphys), 2); % cols = NCSNR, p
for ch = 1:height(channelsEphys)

    fprintf('.');
    if any(ch == badChsOn5), continue; end

    % pull out broadband power, normalize by mean rest power and convert to log power for normality
    BBpowerRestMean = geomean(BBpowerRest(ch, ~badTrsMask(ch, :)));
    BBlogpowerStimSpecial100 = cell(100, 1);
    for ss = 1:100
        BBlogpowerStimSpecial100{ss} = log10(BBpowerStimNan(ch, special100Idxes{ss, 2}) / BBpowerRestMean);
    end
    
    % calculate SNR and pvalue
    [NCSNR, p] = estimateNCSNR(BBlogpowerStimSpecial100);
    NCSNRs(ch, 1) = NCSNR; NCSNRs(ch, 2) = p;
    
end
fprintf('\n');

elecs.BBncsnr = NCSNRs(:, 1);
elecs.BBncsnrP = NCSNRs(:, 2);
writetable(elecs, fullfile(outdir, 'elecsBB2.tsv'), 'FileType', 'text', 'Delimiter', '\t');

NCSNRsSig = NCSNRs(:, 1);
NCSNRsSig(NCSNRs(:, 2) > 0.05) = 0;
plotActivationsToGiftis([], elecs, NCSNRsSig, giis, giiInfs, sulcs, hemis, 0.6, fullfile(outdir, 'broadband'), 'BBncsnr', 'off');

% save a colorbar
fCbar = plotCtxWtsAdjCbar(0.6, 'SigFraction', 0.05, 'SizeJump', 2, 'InsigText', 'Insig.');
saveas(fCbar, fullfile(outdir, 'broadband', 'cbar_BBncsnr.png'));



%% BIPOLAR
%% Scalar average broadband power, calculated for bipolar-reref electrodes
% Add BB rsq and p-value to bipolar electrodes table, save, and also plot mean BB rsq on giftis and inflated giftis

fprintf('Bipolar electrodes analysis\n');

% include bipolar channels that are good on 5+ runs
[goodChsBipOn5, badChsBipOn5] = goodChsNRuns(taskruns.badChsBip, height(elecsBip), 5);

% Calculate PSD for each trial
intervalStim = [0, 0.5]; intervalRest = [-0.5, 0]; % stim = first 500ms after stimulus, rest = last 500 ms before stimulus
pwelchWin = hann(srate/4);
psdStimBip = nan(length(bipolarNames), srate/2+1, height(eventsST));
psdRestBip = nan(length(bipolarNames), srate/2+1, height(eventsST));
for ii = 1:length(bipolarNames)
    [psdStimBip(ii, :, :), f] = pwelch(squeeze(MephysBip(ii, tt >= intervalStim(1) & tt < intervalStim(2), :)), pwelchWin, length(pwelchWin)/2, srate, srate);
    psdRestBip(ii, :, :) = pwelch(squeeze(MephysBip(ii, tt >= intervalRest(1) & tt < intervalRest(2), :)), pwelchWin, length(pwelchWin)/2, srate, srate);
end

% Plot average stationary PSD for all bipolar channels and save
toplot = nan(length(bipolarNames), length(f));
for ii = 1:length(bipolarNames)
    toplot(ii, :) = mean(log10(psdStimBip(ii, :, ~badTrsMaskBip(ii, :))), 3, 'omitnan') - mean(log10(psdRestBip(ii, :, ~badTrsMaskBip(ii, :))), 3, 'omitnan');
end
toplot(badChsBipOn5, :) = nan;
figure('Position', [200, 200, 600, 1200], 'Visible', 'off');
imagescNaN(f, 1:length(bipolarNames), toplot, 'CLim', [-0.2, 0.2]); colormap(cmkjm);
xlim([0, 200]);
set(gca, 'YTick', 1:length(bipolarNames), 'YTickLabels', bipolarNames);
xlabel('Frequency (Hz)');
ylabel('channels');
colorbar;
saveas(gcf, fullfile(outdir, sprintf('PSDsummaryBipolar_%s.png', sub)));

% Calculate average broadband power for each trial (geomean across frequency bins as in Yuasa et al., 2023)
BBpowerStimBip = squeeze(geomean(psdStimBip(:, f >= 70 & f <= 115 | f >= 125 & f <= 170, :), 2));
BBpowerRestBip = squeeze(geomean(psdRestBip(:, f >= 70 & f <= 115 | f >= 125 & f <= 170, :), 2));

% consider normalizing by sessions. But then again, R^2 should still work in the absence of adjusting for any other effects (univariate)
% % Normalize stim and rest power by mean rest power in each session, to account for large impedance differences over time
% sessions = unique(eventsST.session);
% for ii = 1:length(sessions)
%     sessionBaseline = geomean(BBpowerRestBip(:, eventsST.session == sessions(ii)), 2);
%     temp(:, eventsST.session == sessions(ii)) = temp(:, eventsST.session == sessions(ii)) ./ sessionBaseline;
% end

%% Calculate rsq for each channel, across trials good for those channels

rsqBBBip = nan(length(bipolarNames), 2);
for ii = 1:length(bipolarNames)
    [rsq, p] = mnl_rsq(BBpowerStimBip(ii, ~badTrsMaskBip(ii, :)), BBpowerRestBip(ii, ~badTrsMaskBip(ii, :)));
    rsqBBBip(ii, :) = [rsq, p];
end
% set bad channels (in any run) to nan
rsqBBBip(badChsBipOn5, :) = nan;

% Add to bipolar electrodes table and write as output
elecsBip.BBrsq = rsqBBBip(:, 1);
elecsBip.BBrsqP = rsqBBBip(:, 2);
writetable(elecsBip, fullfile(outdir, 'elecsBipBB1.tsv'), 'FileType', 'text', 'Delimiter', '\t');

% Plot broadband Rsq activations on gifti and inflated gifti
plotActivationsToGiftis(elecs, elecsBip, elecsBip.BBrsq, giis, giiInfs, sulcs, hemis, 0.5, fullfile(outdir, 'broadband'), 'BBrsqBip', 'off');

%% Calculate noise ceiling SNR for each bipolar electrode using the special100 images

% find indices of the special100 images
special100Idxes = getSpecial100Idxes(eventsST.stim_file);

% copy BBpowerStimBip, and put nans in bad chxtrs. These nans are removed in estimateNCSNR function
BBpowerStimBipNan = BBpowerStimBip;
BBpowerStimBipNan(badTrsMaskBip) = nan;

rng('default');

NCSNRsBip = nan(length(bipolarNames), 2); % cols = NCSNR, p
for ch = 1:length(bipolarNames)

    fprintf('.');
    if any(ch == badChsBipOn5), continue; end

    % pull out broadband power, normalize by mean rest power and convert to log power for normality
    BBpowerRestMean = geomean(BBpowerRestBip(ch, ~badTrsMaskBip(ch, :)));
    BBlogpowerStimSpecial100 = cell(100, 1);
    for ss = 1:100
        BBlogpowerStimSpecial100{ss} = log10(BBpowerStimBipNan(ch, special100Idxes{ss, 2}) / BBpowerRestMean);
    end
    
    % calculate SNR and pvalue
    [NCSNR, p] = estimateNCSNR(BBlogpowerStimSpecial100);
    NCSNRsBip(ch, 1) = NCSNR; NCSNRsBip(ch, 2) = p;
    
end
fprintf('\n');

elecsBip.BBncsnr = NCSNRsBip(:, 1);
elecsBip.BBncsnrP = NCSNRsBip(:, 2);
writetable(elecsBip, fullfile(outdir, 'elecsBipBB2.tsv'), 'FileType', 'text', 'Delimiter', '\t');

NCSNRsBipSig = NCSNRsBip(:, 1);
NCSNRsBipSig(NCSNRsBip(:, 2) > 0.05) = 0;
plotActivationsToGiftis(elecs, elecsBip, NCSNRsBipSig, giis, giiInfs, sulcs, hemis, 0.6, fullfile(outdir, 'broadband'), 'BBncsnrBip', 'off');



