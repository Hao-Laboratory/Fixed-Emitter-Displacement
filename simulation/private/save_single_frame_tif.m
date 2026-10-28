function save_single_frame_tif(I, outPath)
%SAVE_SINGLE_FRAME_TIF Save one frame as a TIFF image.

assert(~isempty(I), 'Input frame must not be empty.');
assert(ischar(outPath) || isstring(outPath), 'outPath must be a path string.');

outPath = char(outPath);
outDir = fileparts(outPath);
if ~isempty(outDir) && exist(outDir, 'dir') ~= 7
    mkdir(outDir);
end

imwrite(I, outPath, 'tif');
end
