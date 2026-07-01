function dpim_write_table_safe(tbl, path)
%DPIM_WRITE_TABLE_SAFE Create parent directory and write table.
folder = fileparts(path);
if ~isempty(folder) && ~exist(folder, 'dir')
    mkdir(folder);
end
writetable(tbl, path);
end
