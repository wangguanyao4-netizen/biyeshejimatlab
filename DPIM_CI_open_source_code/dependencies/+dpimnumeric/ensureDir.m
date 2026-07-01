function ensureDir(dirPath)
%ensureDir Create a directory if it does not already exist.

if ~exist(dirPath, "dir")
    mkdir(dirPath);
end
end
