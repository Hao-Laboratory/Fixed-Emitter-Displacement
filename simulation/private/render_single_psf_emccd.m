function [image, diagnostic] = render_single_psf_emccd(config, trueXY, conditionId, positionId, repeatId)
%RENDER_SINGLE_PSF_EMCCD Single-frame EMCCD reference used by generation checks.

validateattributes(trueXY, {'numeric'}, {'real','finite','numel',2});
validateattributes(conditionId, {'numeric'}, {'scalar','integer','>=',1,'<=',numel(config.signalADU)});
validateattributes(positionId, {'numeric'}, {'scalar','integer','>=',1,'<=',config.nPositions});
validateattributes(repeatId, {'numeric'}, {'scalar','integer','>=',1,'<=',config.nRepeats});
assert(config.emGain > 0 && config.electronsPerADU > 0);
assert(config.excessNoiseFactorSquared == 2, ...
    'Poisson input and Gamma multiplication implement excess-noise factor squared 2.');
assert(all(trueXY >= 1) && trueXY(1) <= config.imageSizePx(2) && ...
    trueXY(2) <= config.imageSizePx(1), 'True position outside image.');

seed = config.noiseSeedBase + (conditionId-1)*config.conditionSeedStride + ...
    (positionId-1)*config.positionSeedStride + repeatId;
assert(seed < 2^32 && repeatId < config.positionSeedStride && ...
    config.nPositions*config.positionSeedStride < config.conditionSeedStride);
priorRng = rng;
restoreRng = onCleanup(@() rng(priorRng)); %#ok<NASGU>
rng(seed, config.rngAlgorithm);

[X,Y] = meshgrid(1:config.imageSizePx(2), 1:config.imageSizePx(1));
psf = exp(-((X-trueXY(1)).^2+(Y-trueXY(2)).^2)/(2*config.psfSigmaPx^2));
psf = psf / sum(psf(:));
signalExpectedADU = config.signalADU(conditionId) * psf;
lambdaInput = (signalExpectedADU + config.backgroundADU) * ...
    config.electronsPerADU / config.emGain;
inputElectron = poissrnd(lambdaInput);
outputElectron = zeros(size(inputElectron));
positive = inputElectron > 0;
outputElectron(positive) = gamrnd(inputElectron(positive), config.emGain);
outputElectron = outputElectron + config.readNoiseOutputElectron * randn(size(outputElectron));
aduBeforeQuantization = outputElectron / config.electronsPerADU + config.offsetADU;
image = uint16(round(min(max(aduBeforeQuantization,0),65535)));

if nargout > 1
    diagnostic = struct();
    diagnostic.seed = seed;
    diagnostic.psf = psf;
    diagnostic.expectedSignalADU = signalExpectedADU;
    diagnostic.lambdaInputElectron = lambdaInput;
    diagnostic.aduBeforeQuantization = aduBeforeQuantization;
    diagnostic.saturatedPixelCount = nnz(image == 65535);
    diagnostic.zeroPixelCount = nnz(image == 0);
end
end
