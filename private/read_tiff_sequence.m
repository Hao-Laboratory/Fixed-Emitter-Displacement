function [filePaths, fileNames] = read_tiff_sequence(dataDir, ext)
files = dir(fullfile(dataDir, ext));

[~, idx] = sort({files.name});
files = files(idx);

fileNames = {files.name};
fileNames = fileNames(:);
filePaths = cell(numel(files),1);
for i = 1:numel(files)
    filePaths{i} = fullfile(files(i).folder, files(i).name);
end
end
