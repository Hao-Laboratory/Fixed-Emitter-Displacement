function results = run_template(params)
%RUN_TEMPLATE Five-frame templates, MAD displacement estimates and global bridging.
tAllStart = tic;

[filePaths, fileNames] = read_tiff_sequence(params.dataDir, params.ext);
nF = numel(filePaths);

if ~isempty(params.maxFramesOverride)
    nKeep = min(max(1, round(params.maxFramesOverride)), nF);
    filePaths = filePaths(1:nKeep);
    fileNames = fileNames(1:nKeep);
    nF = numel(filePaths);
end

if nF < 5
    error('Need at least 5 frames to build initial template.');
end

nPairs = nF - 1;
templateLen = 5;
templateAnchorInSeg = 3;
templateUpdateInterval = params.templateUpdateInterval;
templateUpdateMode = lower(string(params.templateUpdateMode));

% Timing arrays
% per-frame columns: [imread, detect, localize]
timeFrame = zeros(nF, 3);
time_total = zeros(nF, 1);
time_read = zeros(nF, 1);
time_detect = zeros(nF, 1);
time_localize = zeros(nF, 1);
time_match = zeros(nF, 1);
time_update_total = zeros(nF, 1);
time_update_build = zeros(nF, 1);
time_update_bridge = zeros(nF, 1);
time_update_switch = zeros(nF, 1);
did_update = false(nF, 1);
nCandidates_raw = zeros(nF, 1);
nCandidates_afterBorder = zeros(nF, 1);
nCandidates_final = zeros(nF, 1);
timePair = zeros(nPairs, 2);

logmsg(params, 'Found %d frames.\n', nF);

update_events = cell(0,1);
templates = cell(0,1);
frame_match_ok = false(nF,1);
locs = cell(nF, 1);
frameSize = [NaN NaN];
recentLocs = cell(templateLen, 1);
recentFrameIdx = nan(templateLen, 1);

%% ====== INITIAL 5 FRAMES ======
for t = 1:templateLen
    [xy, ~, frameSize, timeFrame(t,:), detectInfo] = process_frame_localization( ...
        filePaths{t}, t, nF, params, frameSize);
    time_read(t) = timeFrame(t,1);
    time_detect(t) = timeFrame(t,2);
    time_localize(t) = timeFrame(t,3);
    time_total(t) = time_read(t) + time_detect(t) + time_localize(t);
    nCandidates_raw(t) = detectInfo.nCandidates_raw;
    nCandidates_afterBorder(t) = detectInfo.nCandidates_afterBorder;
    nCandidates_final(t) = detectInfo.nCandidates_final;
    locs{t} = xy;
    recentLocs{t} = xy;
    recentFrameIdx(t) = t;
end

%% ====== INITIAL TEMPLATE ======
tt = tic;
[template_cur, initBuildOk, initInfo] = build_template_from_5frames( ...
    recentLocs, recentFrameIdx(:).', params, frameSize, templateAnchorInSeg);
time_init_template_build = toc(tt);
if ~initBuildOk
    error('Initial template build failed: %s', initInfo.reason);
end

template_cur_idx = 1;
framesSinceTemplate = 0;

traj_global_px = nan(nF, 2);
rel_to_template_px = nan(nF, 2);
template_idx_of_frame = nan(nF, 1);
drift_cont = nan(nPairs, 2);
se_cont = nan(nPairs, 2);
n_used = zeros(nPairs, 1);

[local1, ~, nValid1, ok1] = estimate_local_pos_from_template(template_cur.xy_local, locs{1}, params);
if ~ok1
    local1 = [0 0];
end
template_cur.global_xy = -local1;
templates{1,1} = template_cur;
frame_match_ok(1) = ok1;
traj_global_px(1,:) = [0 0];
rel_to_template_px(1,:) = local1;
template_idx_of_frame(1) = template_cur_idx;

for t = 2:templateLen
    idxPair = t - 1;
    tt = tic;
    [localPos, sePos, nValid, okPos] = estimate_local_pos_from_template(template_cur.xy_local, locs{t}, params);
    timePair(idxPair,1) = toc(tt);

    frame_match_ok(t) = okPos;
    if okPos
        rel_to_template_px(t,:) = localPos;
        traj_global_px(t,:) = template_cur.global_xy + localPos;
        se_cont(idxPair,:) = sePos;
    else
        traj_global_px(t,:) = traj_global_px(t-1,:);
        rel_to_template_px(t,:) = traj_global_px(t,:) - template_cur.global_xy;
    end

    template_idx_of_frame(t) = template_cur_idx;
    n_used(idxPair) = nValid;
    drift_cont(idxPair,:) = traj_global_px(t,:);
    time_match(t) = timePair(idxPair,1);
    time_total(t) = time_total(t) + time_match(t);
    logmsg(params, 'Template(init)->Frame %d: matches=%d\n', t, nValid);
end

initialMatchCount = max(1, round(median(n_used(1:templateLen-1), 'omitnan')));
templateMinMatchDynamic = max(1, floor(initialMatchCount * params.templateUpdateMatchFrac));

%% ====== ONLINE MAIN LOOP ======
for t = (templateLen + 1):nF
    [xy, ~, frameSize, timeFrame(t,:), detectInfo] = process_frame_localization( ...
        filePaths{t}, t, nF, params, frameSize);
    time_read(t) = timeFrame(t,1);
    time_detect(t) = timeFrame(t,2);
    time_localize(t) = timeFrame(t,3);
    nCandidates_raw(t) = detectInfo.nCandidates_raw;
    nCandidates_afterBorder(t) = detectInfo.nCandidates_afterBorder;
    nCandidates_final(t) = detectInfo.nCandidates_final;
    locs{t} = xy;
    recentLocs_before = recentLocs;
    recentFrameIdx_before = recentFrameIdx;
    recentLocs = [recentLocs(2:end); {xy}];
    recentFrameIdx = [recentFrameIdx(2:end); t];

    idxPair = t - 1;
    tt = tic;
    [localPos, sePos, nValid, okPos] = estimate_local_pos_from_template(template_cur.xy_local, locs{t}, params);
    timePair(idxPair,1) = toc(tt);

    frame_match_ok(t) = okPos;
    if okPos
        rel_to_template_px(t,:) = localPos;
        traj_global_px(t,:) = template_cur.global_xy + localPos;
        se_cont(idxPair,:) = sePos;
    else
        traj_global_px(t,:) = traj_global_px(t-1,:);
        rel_to_template_px(t,:) = traj_global_px(t,:) - template_cur.global_xy;
    end

    template_idx_of_frame(t) = template_cur_idx;
    n_used(idxPair) = nValid;
    drift_cont(idxPair,:) = traj_global_px(t,:);
    time_match(t) = timePair(idxPair,1);
    time_total(t) = time_read(t) + time_detect(t) + time_localize(t) + time_match(t);
    logmsg(params, 'Template(active)->Frame %d: matches=%d\n', t, nValid);

    framesSinceTemplate = framesSinceTemplate + 1;
    needUpdate = (nValid < templateMinMatchDynamic) || (framesSinceTemplate >= templateUpdateInterval);
    if ~needUpdate
        continue;
    end

    switch templateUpdateMode
        case "prev5_then_current_frame"
            updateLocs = recentLocs_before;
            updateFrameIdx = recentFrameIdx_before(:).';
        case "recent5_then_next_frame"
            updateLocs = recentLocs;
            updateFrameIdx = recentFrameIdx(:).';
    end
    logmsg(params, 'Template update candidate at frame %d (mode=%s, frames=%s).\n', ...
        t, char(templateUpdateMode), mat2str(updateFrameIdx));

    event = struct('frame',t,'low_match',nValid < templateMinMatchDynamic, ...
        'interval_reached',framesSinceTemplate >= templateUpdateInterval, ...
        'matches_before_update',nValid,'minimum_matches',templateMinMatchDynamic, ...
        'frames_since_template',framesSinceTemplate,'frames_used',updateFrameIdx, ...
        'mode',templateUpdateMode,'build_info',struct(),'bridge_info',struct(), ...
        'delta_template_px',[nan nan],'success',false);
    update_events{end+1,1} = event;
    ttUpdate = tic;
    tt = tic;
    [template_new, updOk, updInfo] = build_template_from_5frames( ...
        updateLocs, updateFrameIdx, params, frameSize, templateAnchorInSeg);
    time_update_build(t) = toc(tt);
    update_events{end}.build_info = updInfo;
    if ~updOk || template_new.n_points < 1
        time_update_total(t) = toc(ttUpdate);
        time_total(t) = time_read(t) + time_detect(t) + time_localize(t) + time_match(t) + time_update_total(t);
        logmsg(params, 'Template update skipped at frame %d: %s\n', t, char(string(updInfo.reason)));
        continue;
    end

    tt = tic;
    [delta_template_px, shiftOk, shiftInfo] = estimate_template_shift(template_cur, template_new, params);
    time_update_bridge(t) = toc(tt);
    update_events{end}.bridge_info = shiftInfo;
    update_events{end}.delta_template_px = delta_template_px;
    if ~shiftOk
        time_update_total(t) = toc(ttUpdate);
        time_total(t) = time_read(t) + time_detect(t) + time_localize(t) + time_match(t) + time_update_total(t);
        logmsg(params, 'Template update skipped at frame %d: %s\n', t, char(string(shiftInfo.reason)));
        continue;
    end

    tt = tic;
    template_new.global_xy = template_cur.global_xy + delta_template_px;
    template_cur_idx = template_cur_idx + 1;
    template_cur = template_new;
    templates{end+1,1} = template_cur;
    update_events{end}.success = true;
    framesSinceTemplate = 0;
    if templateUpdateMode == "prev5_then_current_frame"
        ttRematch = tic;
        [localPosNew, sePosNew, nValidNew, okPosNew] = estimate_local_pos_from_template( ...
            template_cur.xy_local, locs{t}, params);
        timePair(idxPair,2) = toc(ttRematch);
        time_match(t) = timePair(idxPair,1) + timePair(idxPair,2);

        frame_match_ok(t) = okPosNew;
        if okPosNew
            rel_to_template_px(t,:) = localPosNew;
            traj_global_px(t,:) = template_cur.global_xy + localPosNew;
            se_cont(idxPair,:) = sePosNew;
            n_used(idxPair) = nValidNew;
        else
            traj_global_px(t,:) = traj_global_px(t-1,:);
            rel_to_template_px(t,:) = traj_global_px(t,:) - template_cur.global_xy;
            se_cont(idxPair,:) = [nan nan];
            n_used(idxPair) = 0;
        end

        template_idx_of_frame(t) = template_cur_idx;
        drift_cont(idxPair,:) = traj_global_px(t,:);
        logmsg(params, ['Template updated at frame %d: nCommon=%d, delta=[%.4f %.4f], ' ...
            'rematch=%d, frames=%s.\n'], ...
            t, shiftInfo.nCommonUsed, delta_template_px(1), delta_template_px(2), ...
            n_used(idxPair), mat2str(updateFrameIdx));
    else
        logmsg(params, 'Template updated at frame %d: nCommon=%d, delta=[%.4f %.4f], frames=%s.\n', ...
            t, shiftInfo.nCommonUsed, delta_template_px(1), delta_template_px(2), mat2str(updateFrameIdx));
    end
    time_update_switch(t) = toc(tt);
    time_update_total(t) = toc(ttUpdate);
    time_total(t) = time_read(t) + time_detect(t) + time_localize(t) + time_match(t) + time_update_total(t);
    did_update(t) = true;
end

%% ====== UNIT CONVERSION ======
scale = params.pixelSizeNm;
drift_cont = drift_cont * scale;
se_cont = se_cont * scale;
cum_cont = [0 0; drift_cont];

%% ====== PACKAGE RESULTS ======
results = struct();
results.params = params;
results.fileNames = fileNames;
results.locs = locs;
results.drift_cont = drift_cont;
results.se_cont = se_cont;
results.n_used = n_used;
results.cum_cont = cum_cont;
results.traj_global = cum_cont;
results.local_pos = rel_to_template_px * scale;
results.rel_to_template_px = rel_to_template_px;
results.template_idx_of_frame = template_idx_of_frame;
results.step_pair_nm = diff(cum_cont,1,1);
results.frame_n_used = [nValid1; n_used];
results.frame_match_ok = frame_match_ok;
results.initialTemplateBuildInfo = initInfo;
results.initialTemplateMatchCount = initialMatchCount;
results.templateMinMatchDynamic = templateMinMatchDynamic;
results.did_update = did_update;
results.updateFrames = find(did_update);
results.update_events = update_events;
results.templates = templates;

summary = report_online_timing_summary( ...
    time_total, time_read, time_detect, time_localize, time_match, ...
    time_update_total, time_update_build, time_update_bridge, time_update_switch, ...
    did_update, time_init_template_build, toc(tAllStart), params);
detectSummary = report_detect_timing_summary( ...
    time_detect, nCandidates_raw, nCandidates_afterBorder, nCandidates_final, params);

results.timing = struct();
results.timing.summary = summary;
results.timing.detectSummary = detectSummary;

if params.printCompletionMessage
    disp('All done.');
end
end

function logmsg(params, varargin)
if params.verbose
    fprintf(varargin{:});
end
end

function [xy, q, frameSize, timeRow, detectInfo] = process_frame_localization(filePath, frameIdx, nF, params, frameSize)
timeRow = zeros(1, 3);
tt = tic;
Iraw = imread(char(filePath));
if lower(string(params.numericType)) == "double"
    I = im2double(Iraw);
else
    I = im2single(Iraw);
end
timeRow(1) = toc(tt);

if frameIdx == 1
    frameSize = size(I);
    frameSize = frameSize(1:2);
end

paramsDetect = params;
paramsDetect.detectFrameIdx = frameIdx;
[candXY, Wdog, ~, ~, detectInfo] = detect_dog_candidates(I, paramsDetect);
timeRow(2) = detectInfo.time_total;

tt = tic;
[xy, q] = localize_weighted_centroid(I, Wdog, candXY, params);
timeRow(3) = toc(tt);

logmsg(params, 'Frame %d/%d: cand=%d, final=%d\n', frameIdx, nF, size(candXY,1), size(xy,1));
end

function [template, ok, info] = build_template_from_5frames(segLocs, frameIdx, params, frameSize, anchorIdx)
[xy_local, ok, infoSeg] = build_template_from_segment(segLocs, params, frameSize, anchorIdx);

template = struct();
template.xy_local = xy_local;
template.ids = (1:size(xy_local,1)).';
template.frames_used = frameIdx(:).';
template.global_xy = [nan nan];
template.n_points = size(xy_local,1);

info = infoSeg;
info.frames_used = template.frames_used;

fs = info.filterStats;
fStart = frameIdx(1);
fEnd = frameIdx(end);
logmsg(params, ['Template build stats [frames %d-%d]: anchor=%d -> presence=%d (-%d) -> ' ...
    'dNN=%d (-%d) -> sigma=%d (-%d) -> edge=%d (-%d), reason=%s\n'], ...
    fStart, fEnd, ...
    fs.anchor_count, fs.after_presence_count, fs.drop_presence, ...
    fs.after_dnn_count, fs.drop_dnn, ...
    fs.after_sigma_count, fs.drop_sigma, ...
    fs.after_edge_count, fs.drop_edge, ...
    char(string(info.reason)));
end

function [delta_template_px, ok, info] = estimate_template_shift(template_old, template_new, params)
delta_template_px = [nan nan];
ok = false;
info = struct('reason', "", 'nCommonRaw', 0, 'nCommonUsed', 0);

if isempty(template_old.xy_local) || isempty(template_new.xy_local)
    info.reason = "empty_template";
    return;
end

bridgeMatchDist = get_template_assoc_dist(params);
if ~isempty(params.templateBridgeMatchDist)
    bridgeMatchDist = max(0, params.templateBridgeMatchDist);
end

[mOld, mNew] = match_points_nn(template_old.xy_local, template_new.xy_local, bridgeMatchDist);
if isempty(mOld)
    info.reason = "no_common_points";
    return;
end

deltaSamples = mNew - mOld;
[dEst, ~, keep] = estimate_drift_continuous(deltaSamples, params.madTau);

info.nCommonRaw = size(deltaSamples, 1);
info.nCommonUsed = sum(keep);
if info.nCommonUsed < params.templateBridgeMinCommonPoints
    info.reason = "insufficient_common_points";
    return;
end
if any(~isfinite(dEst))
    info.reason = "invalid_template_shift";
    return;
end

delta_template_px = dEst;
ok = true;
info.reason = "ok";
end

function [template_locs, ok, info] = build_template_from_segment(segLocs, params, frameSize, anchorIdx)
L = numel(segLocs);
ok = false;
filterStats = struct( ...
    'anchor_count', 0, ...
    'after_presence_count', 0, ...
    'after_sigma_count', 0, ...
    'after_dnn_count', 0, ...
    'after_edge_count', 0, ...
    'drop_presence', 0, ...
    'drop_sigma', 0, ...
    'drop_dnn', 0, ...
    'drop_edge', 0);
info = struct('reason', "", 'nTemplate', 0, 'filterStats', filterStats);
template_locs = zeros(0, 2);

anchor = segLocs{anchorIdx};
if isempty(anchor)
    info.reason = "empty_anchor";
    return;
end

shifted = cell(L, 1);
shifted{anchorIdx} = anchor;
for k = 1:L
    if k == anchorIdx
        continue;
    end
    A = segLocs{k};
    if isempty(A)
        shifted{k} = zeros(0, 2);
        continue;
    end
    [m1, m2] = match_points_nn(A, anchor, get_template_assoc_dist(params));
    if isempty(m1)
        shifted{k} = A;
        continue;
    end
    dispXY = m2 - m1;
    [d_est, ~, keep] = estimate_drift_continuous(dispXY, params.madTau);
    if ~any(keep) || any(~isfinite(d_est))
        d_est = [0 0];
    end
    shifted{k} = A + d_est;
end

anchorPts = shifted{anchorIdx};
nA = size(anchorPts, 1);
filterStats.anchor_count = nA;
obsStack = zeros(L, 2, nA, 'like', anchorPts);
obsCount = ones(nA, 1);
obsFillPos = 2 * ones(nA, 1);
obsStack(1, :, :) = reshape(anchorPts.', 1, 2, nA);
matchedIdxByFrame = cell(L, 1);
matchedPtsByFrame = cell(L, 1);

for k = 1:L
    if k == anchorIdx || isempty(shifted{k})
        continue;
    end
    [~, mK, idxA] = match_points_nn(anchorPts, shifted{k}, get_template_assoc_dist(params));
    matchedIdxByFrame{k} = idxA;
    matchedPtsByFrame{k} = mK;
    if isempty(idxA)
        continue;
    end
    obsCount = obsCount + accumarray(idxA, 1, [nA, 1]);
end

for k = 1:L
    idxA = matchedIdxByFrame{k};
    if isempty(idxA)
        continue;
    end
    mK = matchedPtsByFrame{k};
    for j = 1:numel(idxA)
        ii = idxA(j);
        obsStack(obsFillPos(ii), :, ii) = mK(j,:);
        obsFillPos(ii) = obsFillPos(ii) + 1;
    end
end

mu = zeros(nA, 2, 'like', anchorPts);
sigma_r = zeros(nA, 1, 'like', anchorPts);
keepPresence = false(nA, 1);
minAppearFrames = max(1, min(5, round(params.templateMinAppearFrames)));

for i = 1:nA
    nObs = obsCount(i);
    if nObs < minAppearFrames
        continue;
    end
    pts = reshape(obsStack(1:nObs, :, i), nObs, 2);
    m = mean(pts, 1);
    r = pts - m;
    sr = sqrt(mean(sum(r.^2, 2)));
    mu(i,:) = m;
    sigma_r(i,1) = sr;
    keepPresence(i) = true;
end
mu = mu(keepPresence, :);
sigma_r = sigma_r(keepPresence);
filterStats.after_presence_count = size(mu, 1);
filterStats.drop_presence = max(filterStats.anchor_count - filterStats.after_presence_count, 0);

if isempty(mu)
    info.filterStats = filterStats;
    info.reason = "no_bead_after_presence_filter";
    return;
end

if size(mu, 1) >= 2
    dNN = nearest_neighbor_distance(mu);
    keepDnn = dNN > params.templateDnnMinPxWhenMinSepOff;
    mu = mu(keepDnn, :);
    sigma_r = sigma_r(keepDnn);
end
filterStats.after_dnn_count = size(mu, 1);
filterStats.drop_dnn = max(filterStats.after_presence_count - filterStats.after_dnn_count, 0);

if isempty(mu)
    info.filterStats = filterStats;
    info.reason = "no_bead_after_quality_filters";
    return;
end

[~, ord] = sort(sigma_r, 'ascend');
sigmaKeepFrac = min(max(params.templateSigmaKeepFrac, 0), 1);
keepN = max(1, floor(sigmaKeepFrac * numel(ord)));
mu = mu(ord(1:keepN), :);
filterStats.after_sigma_count = size(mu, 1);
filterStats.drop_sigma = max(filterStats.after_dnn_count - filterStats.after_sigma_count, 0);

edgePad = params.R;
H = frameSize(1);
W = frameSize(2);
keepEdge = mu(:,1) > edgePad & mu(:,1) < (W - edgePad) & ...
           mu(:,2) > edgePad & mu(:,2) < (H - edgePad);
mu = mu(keepEdge, :);
filterStats.after_edge_count = size(mu, 1);
filterStats.drop_edge = max(filterStats.after_sigma_count - filterStats.after_edge_count, 0);

if isempty(mu)
    info.filterStats = filterStats;
    info.reason = "no_bead_after_quality_filters";
    return;
end

template_locs = mu;
ok = true;
info.reason = "ok";
info.nTemplate = size(template_locs, 1);
info.filterStats = filterStats;
end

function [localPos, sePos, nValid, ok] = estimate_local_pos_from_template(template_locs, frame_locs, params)
localPos = [nan nan];
sePos = [nan nan];
nValid = 0;
ok = false;

[m1, m2] = match_points_nn(template_locs, frame_locs, params.maxMatchDist);
if isempty(m1)
    return;
end

dispXY = m2 - m1;
[localPos, sePos, keep] = estimate_drift_continuous(dispXY, params.madTau);
nValid = sum(keep);
ok = all(isfinite(localPos));
end

function d = get_template_assoc_dist(params)
if ~isempty(params.maxMatchDist)
    d = max(2.0, 20 * params.maxMatchDist);
else
    d = 2.0;
end
end

function timing = report_online_timing_summary( ...
    time_total, time_read, time_detect, time_localize, time_match, ...
    time_update_total, time_update_build, time_update_bridge, time_update_switch, ...
    did_update, time_init_template_build, totalWall, params)

validFrame = time_total > 0;
if ~any(validFrame)
    validFrame = true(size(time_total));
end

moduleNames = {'read', 'detect', 'localize', 'match', 'update'};
moduleTotals = [ ...
    sum(time_read, 'omitnan'), ...
    sum(time_detect, 'omitnan'), ...
    sum(time_localize, 'omitnan'), ...
    sum(time_match, 'omitnan'), ...
    sum(time_update_total, 'omitnan')];
moduleTotalSum = sum(moduleTotals);
if moduleTotalSum <= 0
    moduleTotalSum = 1;
end

updateFrames = find(did_update);
nUpdate = numel(updateFrames);

timing = struct();
timing.totalWall = totalWall;
timing.timeInitTemplateBuild = time_init_template_build;
timing.moduleNames = moduleNames;
timing.moduleTotals = moduleTotals;
timing.moduleFractions = moduleTotals / moduleTotalSum;
timing.frameMean = mean(time_total(validFrame), 'omitnan');
timing.frameMedian = median(time_total(validFrame), 'omitnan');
timing.frameMax = max(time_total(validFrame));
timing.frameMaxIdx = find(time_total == timing.frameMax, 1, 'first');
timing.nUpdate = nUpdate;
timing.updateFrames = updateFrames;

if nUpdate > 0
    timing.updateMean = mean(time_update_total(updateFrames), 'omitnan');
    timing.updateMedian = median(time_update_total(updateFrames), 'omitnan');
    timing.updateBuildMean = mean(time_update_build(updateFrames), 'omitnan');
    timing.updateBridgeMean = mean(time_update_bridge(updateFrames), 'omitnan');
    timing.updateSwitchMean = mean(time_update_switch(updateFrames), 'omitnan');
else
    timing.updateMean = 0;
    timing.updateMedian = 0;
    timing.updateBuildMean = 0;
    timing.updateBridgeMean = 0;
    timing.updateSwitchMean = 0;
end

if params.printTimingSummary
    fprintf('\n========== Online Timing Summary ==========\n');
    fprintf('Total wall time        : %.4f s\n', totalWall);
    fprintf('Initial template build : %.4f s\n', time_init_template_build);
    fprintf('\n[Module totals / fractions]\n');
    for i = 1:numel(moduleNames)
        fprintf('  %-10s : %8.4f s   (%6.2f%%)\n', ...
            moduleNames{i}, moduleTotals(i), 100 * moduleTotals(i) / moduleTotalSum);
    end
    fprintf('\n[Per-frame total runtime]\n');
    fprintf('  mean   : %.4f s\n', timing.frameMean);
    fprintf('  median : %.4f s\n', timing.frameMedian);
    fprintf('  max    : %.4f s  (frame %d)\n', timing.frameMax, timing.frameMaxIdx);
    fprintf('\n[Template updates]\n');
    fprintf('  count        : %d\n', nUpdate);
    if nUpdate > 0
        fprintf('  mean total   : %.4f s\n', timing.updateMean);
        fprintf('  median total : %.4f s\n', timing.updateMedian);
        fprintf('  mean build   : %.4f s\n', timing.updateBuildMean);
        fprintf('  mean bridge  : %.4f s\n', timing.updateBridgeMean);
        fprintf('  mean switch  : %.4f s\n', timing.updateSwitchMean);
    end
end
end

function detectTiming = report_detect_timing_summary(time_detect, nCandidates_raw, nCandidates_afterBorder, nCandidates_final, params)
validFrame = time_detect > 0;
if ~any(validFrame)
    validFrame = true(size(time_detect));
end

detectTiming = struct();
detectTiming.detectTotal = sum(time_detect, 'omitnan');
detectTiming.detectMean = mean(time_detect(validFrame), 'omitnan');
detectTiming.detectMax = max(time_detect(validFrame));
detectTiming.nCandidatesRawMean = mean(nCandidates_raw(validFrame), 'omitnan');
detectTiming.nCandidatesRawMedian = median(nCandidates_raw(validFrame), 'omitnan');
detectTiming.nCandidatesRawMax = max(nCandidates_raw(validFrame));
detectTiming.nCandidatesAfterBorderMean = mean(nCandidates_afterBorder(validFrame), 'omitnan');
detectTiming.nCandidatesAfterBorderMedian = median(nCandidates_afterBorder(validFrame), 'omitnan');
detectTiming.nCandidatesAfterBorderMax = max(nCandidates_afterBorder(validFrame));
detectTiming.nCandidatesFinalMean = mean(nCandidates_final(validFrame), 'omitnan');
detectTiming.nCandidatesFinalMedian = median(nCandidates_final(validFrame), 'omitnan');
detectTiming.nCandidatesFinalMax = max(nCandidates_final(validFrame));

if params.printTimingSummary
    fprintf('\n[Detection summary]\n');
    fprintf('  detect_total     : %.4f s | mean %.4f s | max %.4f s\n', ...
        detectTiming.detectTotal, detectTiming.detectMean, detectTiming.detectMax);
    fprintf('  candidates_final : mean %.2f | median %.2f | max %d\n', ...
        detectTiming.nCandidatesFinalMean, detectTiming.nCandidatesFinalMedian, detectTiming.nCandidatesFinalMax);
end
end
