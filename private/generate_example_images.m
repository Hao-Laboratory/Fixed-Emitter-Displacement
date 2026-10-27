function truth = generate_example_images(outputDir)
%GENERATE_EXAMPLE_IMAGES Small deterministic demonstration, not a camera model.
% outputDir is a new directory supplied by example_run.
mkdir(outputDir);
oldRng = rng;
cleanup = onCleanup(@() rng(oldRng));
rng(104729, 'twister');
H = 176; W = 176; nFrames = 25;
[x, y] = meshgrid(1:W, 1:H);
[cx, cy] = meshgrid(16:16:160, 16:16:160);
centres = [cx(:), cy(:)] + 0.2;
stepPx = [0.018, 0.008];
truth = struct('pixelSizeNm',61.9,'cumDriftPx',(0:nFrames-1).' * stepPx);
truth.cumDriftNm = truth.cumDriftPx * truth.pixelSizeNm;
truth.emitterXY = centres;
truth.description = '100 Gaussian emitters, sigma 0.95 px, additive Gaussian noise; demonstration only.';
truth.filenames = cell(nFrames,1);
for t = 1:nFrames
    I = 200 + 3 * randn(H, W);
    xy = centres + truth.cumDriftPx(t,:);
    for k = 1:size(xy,1)
        I = I + 2200 * exp(-((x-xy(k,1)).^2 + (y-xy(k,2)).^2) / (2*0.95^2));
    end
    truth.filenames{t} = sprintf('frame_%04d.tif', t);
    imwrite(uint16(max(0, min(65535, round(I)))), ...
        fullfile(outputDir, truth.filenames{t}), 'Compression', 'none');
end
save(fullfile(outputDir, 'truth.mat'), 'truth');
end
