function [xy, q] = localize_weighted_centroid(I, Wdog, candXY, params)
%LOCALIZE_WEIGHTED_CENTROID ROI-border background subtraction and weighted centroid.
% I: converted image; candXY: integer [x y] peaks; xy: localized coordinates.

profileEnabled = logical(params.enableLocalizeProfile);
time_loc_roi_extract = 0;
time_loc_border_bg = 0;
time_loc_bg_subtract = 0;
time_loc_weight_build = 0;
time_loc_centroid_sum = 0;
time_loc_validity_check = 0;
time_loc_output_write = 0;
if profileEnabled
    tLocTotal = tic;
end

R = params.R;
[H, W] = size(I);
win = 2 * R + 1;
nCand = size(candXY, 1);

returnQuality = logical(params.returnQuality);

qTemplate = struct('sumW',NaN,'bg',NaN,'peak',NaN,'dogPeak',NaN,'weightSource',"raw");

localizeMethod = lower(string(params.localizeMethod));

if nCand == 0
    xy = zeros(0, 2, 'like', I);
    q = repmat(qTemplate, 1, 0);
    if profileEnabled
        update_localize_profile_accum(struct( ...
            'nCalls', 1, ...
            'nCandidatesTotal', 0, ...
            'nLocalizedTotal', 0, ...
            'time_loc_roi_extract', 0, ...
            'time_loc_border_bg', 0, ...
            'time_loc_bg_subtract', 0, ...
            'time_loc_weight_build', 0, ...
            'time_loc_centroid_sum', 0, ...
            'time_loc_validity_check', 0, ...
            'time_loc_output_write', 0, ...
            'time_loc_total', toc(tLocTotal)));
    end
    return;
end

rel = cast(-R:R, 'like', I);
[relY, relX] = ndgrid(rel, rel);
borderMask = false(win, win);
borderMask([1, end], :) = true;
borderMask(:, [1, end]) = true;
borderIdx = find(borderMask);
zeroLikeI = zeros(1, 'like', I);

% Use only the MEX built alongside this function.
localMexFile = fullfile(fileparts(mfilename('fullpath')), ...
    ['localize_weighted_centroid_mex.' mexext]);
if localizeMethod == "mex" && ~returnQuality && ~profileEnabled && ...
        isfile(localMexFile) && ...
        exist('localize_weighted_centroid_mex', 'file') == 3
    [xyAll, nValid] = localize_weighted_centroid_mex(I, candXY, R);
    if nValid <= 0
        xy = zeros(0, 2, 'like', I);
        q = repmat(qTemplate, 1, 0);
    else
        xy = xyAll(1:nValid, :);
        q = repmat(qTemplate, 1, nValid);
    end
    return;
end

xy = zeros(nCand, 2, 'like', I);
if returnQuality
    q = repmat(qTemplate, 1, nCand);
else
    q = repmat(qTemplate, 1, 0);
end
nValid = 0;

for i = 1:size(candXY,1)
    x0 = candXY(i,1);
    y0 = candXY(i,2);
    x0f = cast(x0, 'like', I);
    y0f = cast(y0, 'like', I);

    if profileEnabled
        tt = tic;
    end
    if x0-R < 1 || x0+R > W || y0-R < 1 || y0+R > H
        if profileEnabled
            time_loc_validity_check = time_loc_validity_check + toc(tt);
        end
        continue;
    end
    if profileEnabled
        time_loc_validity_check = time_loc_validity_check + toc(tt);
    end

    if profileEnabled
        tt = tic;
    end
    Iroi = I(y0-R:y0+R, x0-R:x0+R);
    if profileEnabled
        time_loc_roi_extract = time_loc_roi_extract + toc(tt);
    end

    if profileEnabled
        tt = tic;
    end
    bg = median(Iroi(borderIdx));
    if profileEnabled
        time_loc_border_bg = time_loc_border_bg + toc(tt);
    end

    if profileEnabled
        tt = tic;
    end
    IroiBg = Iroi - bg;
    if profileEnabled
        time_loc_bg_subtract = time_loc_bg_subtract + toc(tt);
    end

    if profileEnabled
        tt = tic;
    end
    w = max(IroiBg, zeroLikeI);
    if profileEnabled
        time_loc_weight_build = time_loc_weight_build + toc(tt);
    end

    if profileEnabled
        tt = tic;
    end
    sumW = sum(w(:));
    if sumW <= 0
        if profileEnabled
            time_loc_validity_check = time_loc_validity_check + toc(tt);
        end
        continue;
    end
    if profileEnabled
        time_loc_validity_check = time_loc_validity_check + toc(tt);
    end

    if profileEnabled
        tt = tic;
    end
    xhat = x0f + sum(relX(:).*w(:)) / sumW;
    yhat = y0f + sum(relY(:).*w(:)) / sumW;
    if profileEnabled
        time_loc_centroid_sum = time_loc_centroid_sum + toc(tt);
    end

    if profileEnabled
        tt = tic;
    end
    nValid = nValid + 1;
    xy(nValid,:) = [xhat, yhat];

    if returnQuality
        q(nValid).sumW = sumW;
        q(nValid).bg = bg;
        q(nValid).peak = max(Iroi(:));
        q(nValid).dogPeak = Wdog(y0,x0);
        q(nValid).weightSource = "raw";
    end
    if profileEnabled
        time_loc_output_write = time_loc_output_write + toc(tt);
    end
end

xy = xy(1:nValid, :);
if returnQuality
    q = q(1:nValid);
else
    q = repmat(qTemplate, 1, nValid);
end

if profileEnabled
    update_localize_profile_accum(struct( ...
        'nCalls', 1, ...
        'nCandidatesTotal', nCand, ...
        'nLocalizedTotal', nValid, ...
        'time_loc_roi_extract', time_loc_roi_extract, ...
        'time_loc_border_bg', time_loc_border_bg, ...
        'time_loc_bg_subtract', time_loc_bg_subtract, ...
        'time_loc_weight_build', time_loc_weight_build, ...
        'time_loc_centroid_sum', time_loc_centroid_sum, ...
        'time_loc_validity_check', time_loc_validity_check, ...
        'time_loc_output_write', time_loc_output_write, ...
        'time_loc_total', toc(tLocTotal)));
end
end

function update_localize_profile_accum(delta)
global LOCALIZE_PROFILE_ACCUM
if isempty(LOCALIZE_PROFILE_ACCUM)
    LOCALIZE_PROFILE_ACCUM = init_localize_profile_accum();
end
fields = fieldnames(delta);
for k = 1:numel(fields)
    fn = fields{k};
    LOCALIZE_PROFILE_ACCUM.(fn) = LOCALIZE_PROFILE_ACCUM.(fn) + delta.(fn);
end
end

function accum = init_localize_profile_accum()
accum = struct( ...
    'nCalls', 0, ...
    'nCandidatesTotal', 0, ...
    'nLocalizedTotal', 0, ...
    'time_loc_roi_extract', 0, ...
    'time_loc_border_bg', 0, ...
    'time_loc_bg_subtract', 0, ...
    'time_loc_weight_build', 0, ...
    'time_loc_centroid_sum', 0, ...
    'time_loc_validity_check', 0, ...
    'time_loc_output_write', 0, ...
    'time_loc_total', 0);
end
