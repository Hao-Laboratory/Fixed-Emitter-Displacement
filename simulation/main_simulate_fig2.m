%% Select one Fig. 2 condition.
config = struct();
config.nPoints = 500;        % 100, 300, 500, 1000, 1500, 2000, 3000
labelNm = 5;                 % 1, 3, 5, 7, 8, 10: nominal 2D precision label
config.rep = 1;              % 1 to 5; selects the original emitter coordinates

%% Images and displacement.
config.nFrames = 50;
config.imageSizePx = [1024, 1024];
config.pixelSizeNm = 61.9;
config.driftStepNm = [1, 2];  % nm/frame
config.trueDriftNm = (0:config.nFrames-1).' * config.driftStepNm;
% For a custom trajectory, replace trueDriftNm with an nFrames-by-2 array in nm.

%% PSF and EMCCD camera.
config.psfSigmaPx = 2.18;
config.patchRadiusSigma = 8;
config.backgroundADU = 36.25;
config.offsetADU = 508;
config.emGain = 50;
config.electronsPerADU = 17.3;
config.readNoiseOutputElectron = 130;
config.excessNoiseFactorSquared = 2; % Implemented by Poisson/Gamma sampling.
config.jitterMeanNm = 0;
config.jitterStdAxisNm = 0;          % No independent position jitter in Fig. 2.

%% Historical signal mapping and random seeds.
pointCounts = [100, 300, 500, 1000, 1500, 2000, 3000];
precisionLabels = [1, 3, 5, 7, 8, 10];
signalsADU = [267000, 39400, 18700, 12000, 10200, 7750];
pointCountIndex = find(pointCounts == config.nPoints, 1);
conditionId = find(precisionLabels == labelNm, 1);
assert(~isempty(pointCountIndex) && ~isempty(conditionId) && ismember(config.rep, 1:5), ...
    'Choose one of the listed Fig. 2 conditions and repetitions.');
config.target2DLabelNm = labelNm;
config.signalADU = signalsADU(conditionId); % Expected integrated ADU per emitter.
config.rngAlgorithm = 'twister';
config.sourceRepSeed = 1000 + config.rep;
config.jitterSeed = config.sourceRepSeed;
config.noiseSeedBase = 52000000 + conditionId*100000 + ...
    (config.rep-1)*1000 + (pointCountIndex-1)*1000000;

%% Package data and output.
packageDir = fileparts(mfilename('fullpath'));
config.sourceDataFile = fullfile(packageDir, 'data', 'fig2_source_data.mat');
assert(isfile(config.sourceDataFile), 'Missing data/fig2_source_data.mat.');
assert(exist('poissrnd', 'file') ~= 0 && exist('gamrnd', 'file') ~= 0, ...
    'Requires Statistics and Machine Learning Toolbox.');
sourceData = load(config.sourceDataFile, 'referenceRandom2DNm');
config.referenceSinglePSFRandom2DNm = sourceData.referenceRandom2DNm(conditionId);
config.saveStack = false;       % Optional additional multipage TIFF.
config.makePreviews = false;    % Optional display images; raw TIFFs are unchanged.
config.outputDir = fullfile(packageDir, 'outputs', ...
    sprintf('N%04d_label%02d_rep%03d_%s', config.nPoints, labelNm, ...
    config.rep, datestr(now, 'yyyymmdd_HHMMSS_FFF')));
assert(~isfolder(config.outputDir), 'Output directory already exists.');
mkdir(config.outputDir);

%% Generate images and ground truth.
result = generate_fig2_sequence(config);
fprintf('Images: %s\n', fullfile(config.outputDir, 'frames'));
fprintf('Drift truth: %s\n', fullfile(config.outputDir, 'true_drift.csv'));
