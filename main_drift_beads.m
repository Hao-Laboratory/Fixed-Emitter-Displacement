%% main_drift_beads.m
% Set the image directory and parameters, then run this script.
params = struct();

%% Input images
params.dataDir = fullfile(pwd, 'images'); % Set your TIFF directory.
params.ext = '*.tif';                    % Use zero-padded frame numbers.
params.pixelSizeNm = 61.9;               % nm/px
params.maxFramesOverride = [];          % [] processes all frames.

%% Emitter detection and localization
params.numericType = 'double';          % 'double' or 'single'
params.localizeMethod = 'mex';           % 'mex' or 'matlab'; MATLAB fallback if MEX is absent.
params.returnQuality = false;           % true uses the MATLAB localizer.
params.enableLocalizeProfile = false;   % true uses the MATLAB localizer.
params.detectThresh = [];               % [] is automatic; manual values use converted intensities.
params.kSigma = 6;
params.thresholdMethod = 'sampled';      % 'sampled' or 'exact'
params.thresholdSampleStride = 4;
params.detectThreshRefresh = 1;
params.localMaxWindowPx = 3;
params.R = 3;                           % Centroid ROI half-width, px.

%% Emitter matching and outlier rejection
params.maxMatchDist = 0.2;               % Frame-to-template radius, px.
params.madTau = 3;                       % <=0 disables MAD rejection.

%% Template construction, update and bridging
params.templateMinAppearFrames = 3;
params.templateSigmaKeepFrac = 0.85;
params.templateDnnMinPxWhenMinSepOff = 6; % d_sep: strict DNN threshold on template means, px.
params.templateUpdateMatchFrac = 0.83;
params.templateUpdateInterval = 150;
params.templateBridgeMatchDist = 2;
params.templateBridgeMinCommonPoints = 50;
params.templateUpdateMode = "prev5_then_current_frame";
% Alternative: "recent5_then_next_frame".

% Console output.
params.verbose = false;
params.printTimingSummary = false;
params.printCompletionMessage = true;

%% Run and display
assert(isfolder(params.dataDir), 'Data directory not found: %s', params.dataDir);
assert(all(cellfun(@(f) exist(f, 'file') ~= 0, ...
    {'imdilate','knnsearch','mad','im2double','im2single'})), ...
    'Requires Image Processing Toolbox and Statistics and Machine Learning Toolbox.');
results = run_template(params);
trajectoryFig = plot_trajectory_2d(results);
