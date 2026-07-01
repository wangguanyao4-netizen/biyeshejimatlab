function writeJson(filePath, value)
%writeJson Write a JSON file, using pretty printing when supported.

[parentDir, ~, ~] = fileparts(filePath);
if strlength(string(parentDir)) > 0
    dpimnumeric.ensureDir(parentDir);
end

try
    jsonText = jsonencode(value, PrettyPrint=true);
catch
    jsonText = jsonencode(value);
end

dpimnumeric.writeText(filePath, sprintf("%s\n", char(jsonText)));
end
