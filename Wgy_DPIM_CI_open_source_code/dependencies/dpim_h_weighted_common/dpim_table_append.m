function tbl = dpim_table_append(tbl, row)
%DPIM_TABLE_APPEND Append one table row to a possibly empty table.
% MATLAB does not allow vertical concatenation between table() with zero
% variables and a nonempty table in many releases.  This helper makes all
% level scripts version-robust.
if isempty(tbl) || width(tbl) == 0
    tbl = row;
else
    tbl = [tbl; row]; %#ok<AGROW>
end
end
