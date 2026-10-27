function [candXY, Wdog, threshUsed, G1, detectInfo] = detect_dog_candidates(I, params)
%DETECT_DOG_CANDIDATES Raw-image threshold and local maxima; historical name.
% params contains the configured detection fields and the current frame index.
tDetectTotal = tic;
persistent threshState

kSigma = double(params.kSigma);
manualThresh = params.detectThresh;
detectThreshRefresh = max(1, round(double(params.detectThreshRefresh)));
frameIdx = max(1, round(double(params.detectFrameIdx)));
thresholdMethod = lower(string(params.thresholdMethod));
thresholdSampleStride = max(1, round(params.thresholdSampleStride));
localMaxWindowPx = params.localMaxWindowPx;

Idetect = double(I);
Wdog = Idetect;
G1 = Idetect;
bg = nan;
sigma_n = nan;
thresholdSampleCount = 0;
thresholdSampleFraction = 0;

signature = build_state_signature(kSigma, thresholdMethod, thresholdSampleStride, ...
    detectThreshRefresh, manualThresh);
if isempty(threshState) || frameIdx == 1 || ...
        ~strcmp(threshState.signature, signature) || ...
        ~isequal(threshState.imageSize, size(Idetect))
    threshState = struct('signature',signature,'imageSize',size(Idetect), ...
        'hasAutoStats',false,'bg',nan,'sigma_n',nan,'threshUsed',nan, ...
        'thresholdSampleCount',0,'thresholdSampleFraction',0,'lastFrameIdx',0);
end
if frameIdx <= threshState.lastFrameIdx
    threshState.hasAutoStats = false;
end

if isempty(manualThresh)
    didThreshRefresh = ~threshState.hasAutoStats || ...
        mod(frameIdx - 1, detectThreshRefresh) == 0 || ~isfinite(threshState.threshUsed);
    if didThreshRefresh
        switch thresholdMethod
            case "exact"
                detectVec = Idetect(:);
                thresholdSampleCount = numel(detectVec);
                thresholdSampleFraction = 1;
            case "sampled"
                detectSample = Idetect(1:thresholdSampleStride:end, 1:thresholdSampleStride:end);
                detectVec = detectSample(:);
                thresholdSampleCount = numel(detectVec);
                thresholdSampleFraction = thresholdSampleCount / numel(Idetect);
        end
        bg = median(detectVec);
        sigma_n = 1.4826 * mad(detectVec, 1);
        threshUsed = bg + kSigma * sigma_n;
        threshState.bg = bg;
        threshState.sigma_n = sigma_n;
        threshState.threshUsed = threshUsed;
        threshState.thresholdSampleCount = thresholdSampleCount;
        threshState.thresholdSampleFraction = thresholdSampleFraction;
        threshState.hasAutoStats = true;
    else
        bg = threshState.bg;
        sigma_n = threshState.sigma_n;
        threshUsed = threshState.threshUsed;
        thresholdSampleCount = threshState.thresholdSampleCount;
        thresholdSampleFraction = threshState.thresholdSampleFraction;
    end
else
    didThreshRefresh = false;
    threshUsed = double(manualThresh);
    threshState.hasAutoStats = false;
end
threshState.lastFrameIdx = frameIdx;

mx = imdilate(Idetect, true(localMaxWindowPx));
pk = (Idetect == mx);
cand = pk & (Idetect > threshUsed);
[r, c] = find(cand);
candXY = [c, r];
nCandidates = size(candXY, 1);
detectInfo = struct('method',"raw_thresh_localmax",'kSigma',kSigma, ...
    'thresholdMethod',thresholdMethod,'thresholdSampleStride',thresholdSampleStride, ...
    'localMaxWindowPx',localMaxWindowPx,'detectThreshRefresh',detectThreshRefresh, ...
    'frameIdx',frameIdx,'didThreshRefresh',didThreshRefresh,'bg',bg,'sigma_n',sigma_n, ...
    'thresholdSampleCount',thresholdSampleCount,'thresholdSampleFraction',thresholdSampleFraction, ...
    'time_total',toc(tDetectTotal),'nCandidates_raw',nCandidates, ...
    'nCandidates_afterBorder',nCandidates,'nCandidates_final',nCandidates);
end

function sig = build_state_signature(kSigma, thresholdMethod, stride, refresh, manualThresh)
% Cache identity retains the original precision and reset rules.
if isempty(manualThresh)
    manualSig = 'auto';
else
    manualSig = sprintf('manual:%.16g', double(manualThresh(1)));
end
sig = sprintf('k=%.16g|m=%s|s=%d|r=%d|%s', ...
    kSigma, char(thresholdMethod), stride, refresh, manualSig);
end
