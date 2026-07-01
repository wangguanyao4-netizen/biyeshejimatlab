function writeText(filePath, textValue)
%writeText Write text to a UTF-8 file.

[parentDir, ~, ~] = fileparts(filePath);
if strlength(string(parentDir)) > 0
    dpimnumeric.ensureDir(parentDir);
end

fid = fopen(filePath, "w", "n", "UTF-8");
if fid < 0
    error("Could not open %s for writing.", filePath);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, "%s", char(string(textValue)));
end
