function result=generate_fig2_sequence(config)
%GENERATE_FIG2_SEQUENCE One historical EMCCD sequence and its known displacement.
% Use the first N historical points; add drift and optional independent jitter.
% Sum emitter PSFs before applying camera noise, background and offset once.
priorRng=rng; restoreRng=onCleanup(@() rng(priorRng)); %#ok<NASGU>
outDir=config.outputDir;
assert(config.nPoints>=1 && config.nPoints<=3000 && config.nPoints==fix(config.nPoints) && ...
    config.nFrames>=2 && config.nFrames==fix(config.nFrames) && config.excessNoiseFactorSquared==2);
signalTolerance=max(1e-7,1e-12*config.nPoints*config.signalADU);
assert(isfolder(outDir) && ~isfile(fullfile(outDir,'sequence_noisy.tif')));
source=load(config.sourceDataFile);
old=struct('basePositionsPx',source.basePositionsPx(:,:,config.rep), ...
    'imageSize',source.imageSize,'pixelSizeNm',source.pixelSizeNm, ...
    'repSeed',source.repSeeds(config.rep),'trueDriftNm',source.historicalDriftNm);
assert(isequal(old.imageSize(:).',[1024 1024]) && old.pixelSizeNm==61.9);
assert(size(old.basePositionsPx,1)==3000 && isfinite(old.repSeed));
assert(size(old.basePositionsPx,1)>=config.nPoints);
assert(old.repSeed==config.sourceRepSeed);
assert(isequal(config.imageSizePx*config.pixelSizeNm,old.imageSize(:).'*old.pixelSizeNm));
basePositionsPxFull=(old.basePositionsPx-0.5)*old.pixelSizeNm/config.pixelSizeNm+0.5;
basePositionsPx=basePositionsPxFull(1:config.nPoints,:);
if isfield(config,'trueDriftNm')
    % Explicit trajectories also support custom motion and frame counts.
    validateattributes(config.trueDriftNm,{'numeric'},{'real','finite','size',[config.nFrames,2]});
    trueDriftNm=config.trueDriftNm;
else
    trueDriftNm=(0:config.nFrames-1).'*config.driftStepNm;
    assert(config.nFrames==50 && isequal(trueDriftNm,old.trueDriftNm));
end
rng(config.jitterSeed,config.rngAlgorithm);
jitterNm=config.jitterMeanNm+config.jitterStdAxisNm.*randn(config.nPoints,2,config.nFrames);
truePositionsPx=zeros(config.nPoints,2,config.nFrames);
for t=1:config.nFrames
    truePositionsPx(:,:,t)=basePositionsPx+trueDriftNm(t,:)/config.pixelSizeNm+ ...
        jitterNm(:,:,t)/config.pixelSizeNm;
end
assert(all(reshape(truePositionsPx(:,1,:),[],1)>=1) && ...
    all(reshape(truePositionsPx(:,1,:),[],1)<=config.imageSizePx(2)) && ...
    all(reshape(truePositionsPx(:,2,:),[],1)>=1) && ...
    all(reshape(truePositionsPx(:,2,:),[],1)<=config.imageSizePx(1)));

packageDir=fileparts(fileparts(mfilename('fullpath')));
sourceNames={'data/fig2_source_data.mat';'main_simulate_fig2.m'; ...
    'private/generate_fig2_sequence.m';'private/render_single_psf_emccd.m'; ...
    'private/save_single_frame_tif.m';'private/emccd_thompson_precision.m'};
sourcePaths=cellfun(@(name) fullfile(packageDir,name),sourceNames,'UniformOutput',false);
sourceHashes=cellfun(@local_file_hash,sourcePaths,'UniformOutput',false);
provenance=table(string(sourceNames),string(sourceHashes),'VariableNames',{'path','sha256'});
writetable(provenance,fullfile(outDir,'source_provenance.csv'));

% Check the single-PSF reference and two probe frames before generation.
validation=local_validate(config,basePositionsPx,basePositionsPxFull,old,jitterNm,truePositionsPx,trueDriftNm);
probeImages=cell(2,1); probeDiagnostics=cell(2,1);
for t=1:2
    [probeImages{t},probeDiagnostics{t}]=local_render(config,truePositionsPx(:,:,t),config.noiseSeedBase+t);
end
baselineADU=config.backgroundADU+config.offsetADU;
backgroundVarianceADU2=2*config.backgroundADU*config.emGain/config.electronsPerADU+ ...
    (config.readNoiseOutputElectron/config.electronsPerADU)^2;
bgValues=[];
for t=1:2
    d=probeDiagnostics{t}; mask=d.signalExpectedADU<1e-10;
    pixels=double(probeImages{t}); bgValues=[bgValues;pixels(mask)]; %#ok<AGROW>
    validation=local_check(validation,sprintf('probe_%d_integrated_signal',t), ...
        abs(sum(d.signalExpectedADU(:))-config.nPoints*config.signalADU)<signalTolerance);
    validation=local_check(validation,sprintf('probe_%d_no_saturation_or_clipping',t), ...
        ~any(d.aduBeforeQuantization(:)<0 | d.aduBeforeQuantization(:)>65535));
end
validation=local_check(validation,'background_mean_once',abs(mean(bgValues)-baselineADU)<0.3);
validation=local_check(validation,'background_variance_EM_model', ...
    abs(var(bgValues,0)/backgroundVarianceADU2-1)<0.03);
writetable(validation,fullfile(outDir,'validation_results.csv'));
assert(all(validation.passed),'Small-sample validation failed; sequence generation not started.');
fprintf('Small-sample checks passed: %d/%d. Generating one %d-frame sequence.\n', ...
    nnz(validation.passed),height(validation),config.nFrames);

frameDir=fullfile(outDir,'frames'); mkdir(frameDir);
stackPath=fullfile(outDir,'sequence_noisy.tif');
rows=repmat(struct('frame_id',0,'noise_seed',0,'raw_total_adu',0,'net_total_adu',0, ...
    'expected_signal_total_adu',0,'peak_adu',0,'saturated_pixels',0,'zero_pixels',0, ...
    'clipped_before_quantization',0,'sha256',""),config.nFrames,1);
for t=1:config.nFrames
    if t<=2
        image=probeImages{t}; d=probeDiagnostics{t};
    else
        [image,d]=local_render(config,truePositionsPx(:,:,t),config.noiseSeedBase+t);
    end
    assert(abs(sum(d.signalExpectedADU(:))-config.nPoints*config.signalADU)<signalTolerance);
    framePath=fullfile(frameDir,sprintf('frame_%03d.tif',t));
    save_single_frame_tif(image,framePath);
    assert(isequal(imread(framePath),image),'Individual TIFF pixel mismatch.');
    if t==1
        if config.saveStack, imwrite(image,stackPath,'tif','Compression','none'); end
        firstNoisy=image; firstExpectedADU=d.signalExpectedADU+baselineADU;
        firstNoiseless=uint16(round(firstExpectedADU));
        if config.makePreviews
            save_single_frame_tif(firstNoiseless,fullfile(outDir,'first_frame_noiseless.tif'));
            assert(isequal(imread(fullfile(outDir,'first_frame_noiseless.tif')),firstNoiseless));
        end
    else
        if config.saveStack, imwrite(image,stackPath,'tif','WriteMode','append','Compression','none'); end
    end
    values=double(image(:));
    rows(t)=struct('frame_id',t,'noise_seed',config.noiseSeedBase+t, ...
        'raw_total_adu',sum(values),'net_total_adu',sum(values)-numel(image)*baselineADU, ...
        'expected_signal_total_adu',sum(d.signalExpectedADU(:)), ...
        'peak_adu',max(values),'saturated_pixels',nnz(image==65535),'zero_pixels',nnz(image==0), ...
        'clipped_before_quantization',nnz(d.aduBeforeQuantization<0 | d.aduBeforeQuantization>65535), ...
        'sha256',string(local_hash(typecast(image(:),'uint8'))));
    if mod(t,10)==0, fprintf('Generated and checked %d/%d frames.\n',t,config.nFrames); end
end
if config.saveStack
    info=imfinfo(stackPath); assert(numel(info)==config.nFrames);
    for t=1:config.nFrames
        page=imread(stackPath,t,'Info',info);
        assert(isa(page,'uint16') && isequal(size(page),config.imageSizePx) && info(t).BitDepth==16);
        assert(strcmp(local_hash(typecast(page(:),'uint8')),rows(t).sha256),'Stack pixel mismatch.');
    end
else
    assert(~isfile(stackPath));
end
frameSummary=struct2table(rows);
assert(all(frameSummary.saturated_pixels==0 & frameSummary.clipped_before_quantization==0));
writetable(frameSummary,fullfile(outDir,'frame_summary.csv'));
validation=local_check(validation,sprintf('all_%d_individual_TIFF_pixels_verified_and_stack_policy_observed',config.nFrames),true);
validation=local_check(validation,sprintf('all_%d_no_saturation_or_clipping',config.nFrames),true);

% Physical coordinates use the upper-left image edge as the nm origin.
frameId=kron((1:config.nFrames).',ones(config.nPoints,1));
particleId=repmat((1:config.nPoints).',config.nFrames,1);
x=reshape(truePositionsPx(:,1,:),[],1); y=reshape(truePositionsPx(:,2,:),[],1);
jx=reshape(jitterNm(:,1,:),[],1); jy=reshape(jitterNm(:,2,:),[],1);
truth=table(frameId,particleId,basePositionsPx(particleId,1),basePositionsPx(particleId,2), ...
    trueDriftNm(frameId,1),trueDriftNm(frameId,2),jx,jy,x,y, ...
    (x-0.5)*config.pixelSizeNm,(y-0.5)*config.pixelSizeNm, ...
    repmat(config.signalADU,numel(x),1),repmat(config.psfSigmaPx,numel(x),1), ...
    'VariableNames',{'frame_id','particle_id','base_x_px','base_y_px','drift_x_nm','drift_y_nm', ...
    'jitter_x_nm','jitter_y_nm','true_x_px','true_y_px','true_x_nm','true_y_nm','signal_adu','psf_sigma_px'});
writetable(truth,fullfile(outDir,'ground_truth.csv'));
writetable(table((1:config.nFrames).',trueDriftNm(:,1),trueDriftNm(:,2), ...
    'VariableNames',{'frame_id','drift_x_nm','drift_y_nm'}),fullfile(outDir,'true_drift.csv'));

theoryConfig=config; theoryConfig.target2DNm=config.target2DLabelNm;
theoryConfig.theoryVerificationToleranceNm=0.001;
theory=emccd_thompson_precision(theoryConfig);
writetable(theory,fullfile(outDir,'reference_single_psf_theory.csv'));
metadata=struct('config',config,'pointCount',config.nPoints,'nFrames',config.nFrames, ...
    'imageSize',config.imageSizePx,'pixelSizeNm',config.pixelSizeNm, ...
    'basePositionsPxFull',basePositionsPxFull,'basePositionsPx',basePositionsPx, ...
    'selectedPointIndices',(1:config.nPoints).','sourceRepSeed',old.repSeed, ...
    'trueDriftNm',trueDriftNm,'jitterNm',jitterNm,'truePositionsPx',truePositionsPx, ...
    'truePositionsNm',(truePositionsPx-0.5).*config.pixelSizeNm, ...
    'signalADU',repmat(config.signalADU,config.nPoints,1), ...
    'psfSigmaPx',repmat(config.psfSigmaPx,config.nPoints,1), ...
    'firstExpectedADU',firstExpectedADU,'backgroundVarianceADU2',backgroundVarianceADU2, ...
    'localizationPerformed',false,'coordinateConvention','1-based pixel centers; nm origin at image edge');
if config.makePreviews
    local_preview(config,firstNoisy,firstNoiseless,truePositionsPx(:,:,1),backgroundVarianceADU2);
end
for k=1:numel(sourcePaths)
    assert(strcmp(local_file_hash(sourcePaths{k}),sourceHashes{k}),'Source changed during generation.');
end
validation=local_check(validation,'all_source_hashes_unchanged',true);
writetable(validation,fullfile(outDir,'validation_results.csv'));
save(fullfile(outDir,'metadata.mat'),'metadata','frameSummary','validation','theory','provenance','-v7.3');
save(fullfile(outDir,'simulation_config.mat'),'config','-v7.3');
local_report(config,frameSummary,validation,theory,jitterNm);
result=struct('outputDir',outDir,'nFrames',config.nFrames,'nPoints',config.nPoints, ...
    'checksPassed',height(validation),'saturatedPixels',sum(frameSummary.saturated_pixels), ...
    'localizationPerformed',false);
end

function [image,diagnostic]=local_render(config,trueXY,seed)
signalExpectedADU=zeros(config.imageSizePx);
% Preserve local-grid sampling and per-emitter normalization in source order.
radius=ceil(config.patchRadiusSigma*config.psfSigmaPx);
for i=1:size(trueXY,1)
    x1=max(1,floor(trueXY(i,1))-radius); x2=min(config.imageSizePx(2),floor(trueXY(i,1))+radius);
    y1=max(1,floor(trueXY(i,2))-radius); y2=min(config.imageSizePx(1),floor(trueXY(i,2))+radius);
    [X,Y]=meshgrid(x1:x2,y1:y2);
    psf=exp(-((X-trueXY(i,1)).^2+(Y-trueXY(i,2)).^2)/(2*config.psfSigmaPx^2));
    psf=psf/sum(psf(:));
    signalExpectedADU(y1:y2,x1:x2)=signalExpectedADU(y1:y2,x1:x2)+config.signalADU*psf;
end
rng(seed,config.rngAlgorithm);
% Poisson input, Gamma EM multiplication and Gaussian output read noise.
lambdaInput = (signalExpectedADU + config.backgroundADU) * ...
    config.electronsPerADU / config.emGain;
inputElectron = poissrnd(lambdaInput);
outputElectron = zeros(size(inputElectron));
positive = inputElectron > 0;
outputElectron(positive) = gamrnd(inputElectron(positive), config.emGain);
outputElectron = outputElectron + config.readNoiseOutputElectron * randn(size(outputElectron));
aduBeforeQuantization = outputElectron / config.electronsPerADU + config.offsetADU;
image = uint16(round(min(max(aduBeforeQuantization,0),65535)));
diagnostic=struct('signalExpectedADU',signalExpectedADU,'lambdaInputElectron',lambdaInput, ...
    'aduBeforeQuantization',aduBeforeQuantization);
end

function checks=local_validate(config,base,baseFull,old,jitter,positions,drift)
checks=table(strings(0,1),false(0,1),'VariableNames',{'test','passed'});
checks=local_check(checks,'same_physical_FOV', ...
    isequal(config.imageSizePx*config.pixelSizeNm,old.imageSize(:).'*old.pixelSizeNm));
checks=local_check(checks,'all_3000_source_coordinates_preserved_physically', ...
    max(abs((baseFull(:)-0.5)*config.pixelSizeNm-(old.basePositionsPx(:)-0.5)*old.pixelSizeNm))<1e-10);
checks=local_check(checks,'first_N_source_points_selected',isequal(base,baseFull(1:config.nPoints,:)));
reconstructed=zeros(size(positions));
for t=1:config.nFrames
    reconstructed(:,:,t)=base+drift(t,:)/config.pixelSizeNm+jitter(:,:,t)/config.pixelSizeNm;
end
checks=local_check(checks,'jitter_applied_to_true_coordinates_not_cumulative',isequal(positions,reconstructed));
for axisId=1:2
    v=reshape(jitter(:,axisId,:),[],1);
    checks=local_check(checks,sprintf('jitter_axis%d_sample_mean',axisId),abs(mean(v))<0.08);
    checks=local_check(checks,sprintf('jitter_axis%d_sample_std',axisId),abs(std(v,0)-config.jitterStdAxisNm)<0.08);
end
probe=config; probe.imageSizePx=[16 16]; probe.nPositions=1; probe.nRepeats=1;
probe.noiseSeedBase=1233; probe.conditionSeedStride=1000000; probe.positionSeedStride=1000;
xy=[8.35 8.72];
[reference,referenceD]=render_single_psf_emccd(probe,xy,1,1,1);
[actual,d]=local_render(probe,xy,1234);
checks=local_check(checks,'single_PSF_exact_pixel_match_original_renderer',isequal(actual,reference));
checks=local_check(checks,'single_PSF_exact_expected_signal_match', ...
    isequal(d.signalExpectedADU,referenceD.expectedSignalADU));
checks=local_check(checks,'single_PSF_exact_prequantization_noise_match', ...
    isequal(d.aduBeforeQuantization,referenceD.aduBeforeQuantization));
% Compare the local PSF patch with full-canvas normalization.
probe=config; xy=positions(1,:,1);
[~,d]=local_render(probe,xy,1234);
[X,Y]=meshgrid(1:config.imageSizePx(2),1:config.imageSizePx(1));
g=exp(-((X-xy(1)).^2+(Y-xy(2)).^2)/(2*config.psfSigmaPx^2));
g=config.signalADU*g/sum(g(:));
checks=local_check(checks,'8sigma_patch_matches_full_canvas_PSF',max(abs(g(:)-d.signalExpectedADU(:)))<1e-9);
end

function checks=local_check(checks,name,passed)
checks=[checks;table(string(name),logical(passed),'VariableNames',{'test','passed'})];
end

function digest=local_hash(bytes)
md=java.security.MessageDigest.getInstance('SHA-256'); md.update(bytes);
raw=typecast(md.digest(),'uint8'); digest=lower(reshape(dec2hex(raw,2).',1,[]));
end

function digest=local_file_hash(path)
fid=fopen(path,'r'); assert(fid>=0); cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
digest=local_hash(fread(fid,Inf,'*uint8'));
end

function local_preview(config,noisy,clean,xy,bgVar)
outDir=config.outputDir; baseline=config.backgroundADU+config.offsetADU;
limits=[floor((baseline-3*sqrt(bgVar))/10)*10,ceil(double(max(noisy(:)))/50)*50];
for k=1:2
    if k==1, pixels=noisy; name='first_frame_noisy'; else, pixels=clean; name='first_frame_noiseless'; end
    display=uint8(round(255*min(max((double(pixels)-limits(1))/diff(limits),0),1)));
    imwrite(display,fullfile(outDir,[name,'_display.png']));
    imwrite(repelem(display,2,2),fullfile(outDir,[name,'_display_2x.png']));
end
dist=hypot(xy(:,1)-xy(:,1).',xy(:,2)-xy(:,2).');
dist(1:config.nPoints+1:end)=Inf; [~,id]=max(min(dist,[],2));
xc=round(xy(id,1)); yc=round(xy(id,2));
xs=max(1,xc-16):min(config.imageSizePx(2),xc+15);
ys=max(1,yc-16):min(config.imageSizePx(1),yc+15);
save_single_frame_tif(noisy(ys,xs),fullfile(outDir,'first_frame_isolated_crop_raw.tif'));
crop=uint8(round(255*min(max((double(noisy(ys,xs))-limits(1))/diff(limits),0),1)));
imwrite(repelem(crop,16,16),fullfile(outDir,'first_frame_isolated_crop_display_16x.png'));
fig=figure('Visible','off','Color','w','Position',[100 100 1150 900]);
data={clean,noisy,clean(ys,xs),noisy(ys,xs)};
names={sprintf('Noiseless: N=%d',config.nPoints),sprintf('EMCCD noisy: N=%d',config.nPoints), ...
    'Isolated crop: noiseless','Isolated crop: noisy'};
for k=1:4
    ax=subplot(2,2,k,'Parent',fig); imagesc(ax,data{k},limits); axis(ax,'image');
    colormap(ax,gray(256)); cb=colorbar(ax); cb.Label.String='ADU';
    title(ax,names{k}); xlabel(ax,'x (px)'); ylabel(ax,'y (px)');
end
exportgraphics(fig,fullfile(outDir,'first_frame_overview.png'),'Resolution',150,'BackgroundColor','white');
close(fig);
displayInfo=struct('limitsADU',limits,'isolatedParticleId',id,'cropXIndices',xs,'cropYIndices',ys, ...
    'minimumPairDistancePx',min(dist(:)));
save(fullfile(outDir,'display_settings.mat'),'displayInfo');
end

function local_report(config,summary,checks,theory,jitter)
fid=fopen(fullfile(config.outputDir,'README.md'),'w','n','UTF-8');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'# Fig. 2 simulation sequence\n\n');
fprintf(fid,'%d emitters; label %g nm; repetition %d; %d frames.\n\n', ...
    config.nPoints,config.target2DLabelNm,config.rep,config.nFrames);
fprintf(fid,'%d x %d pixels; %g nm/pixel; PSF sigma %g pixels.\n\n', ...
    config.imageSizePx,config.pixelSizeNm,config.psfSigmaPx);
fprintf(fid,'Signal %g ADU/emitter; background %g ADU/pixel; offset %g ADU.\n\n', ...
    config.signalADU,config.backgroundADU,config.offsetADU);
fprintf(fid,'EM gain %g; conversion %g e-/ADU; output read noise %g e-.\n\n', ...
    config.emGain,config.electronsPerADU,config.readNoiseOutputElectron);
fprintf(fid,'Known trajectory: true_drift.csv. Per-emitter coordinates: ground_truth.csv.\n\n');
fprintf(fid,'Camera noise seeds: %d + frame index.\n\n',config.noiseSeedBase);
fprintf(fid,'Generation checks passed: %d/%d. Saturated pixels: %d.\n', ...
    nnz(checks.passed),height(checks),sum(summary.saturated_pixels));
end
