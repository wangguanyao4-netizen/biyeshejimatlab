function row = dpim_weight_summary_row(curve, problemName, methodName)
%DPIM_WEIGHT_SUMMARY_ROW Return a numeric/string row as a struct.
w = curve.weights(:);
wd = curve.weightData;
row = struct();
row.problem = string(problemName);
row.method = string(methodName);
row.curve_index = curve.curve_index;
row.point_actual_method = string(getfield_or(curve, 'point_actual_method', 'unknown'));
row.point_fallback_used = logical(getfield_or(curve, 'point_fallback_used', false));
row.weight_source = string(curve.weight_source);
row.n = numel(w);
row.sum_weights = sum(w);
row.min_w = min(w);
row.max_w = max(w);
row.mean_w = mean(w);
row.std_w = std(w);
row.cv_w = std(w) / mean(w);
s2 = sum(w.^2);
s3 = sum(w.^3);
s4 = sum(w.^4);
row.s2_w = s2;
row.s3_w = s3;
row.s4_w = s4;
row.rho3_w = s3 / max(s2, realmin)^(3/2);
row.rho4_w = s4 / max(s2, realmin)^2;
row.n2_eff_w = 1 / max(s2, realmin);
row.n3_eff_w = max(s2, realmin)^3 / max(s3, realmin)^2;
row.n4_eff_w = max(s2, realmin)^2 / max(s4, realmin);
row.ess_w = row.n2_eff_w;
row.ess_ratio = row.ess_w / numel(w);
row.l1_from_equal = sum(abs(w - 1 / numel(w)));
row.max_over_equal = max(w) / (1 / numel(w));
row.empty_cell_count = getfield_or(wd, 'empty_cell_count', NaN); %#ok<GFLD>
row.auxiliary_sample_count = getfield_or(wd, 'auxiliary_sample_count', NaN); %#ok<GFLD>
end

function v = getfield_or(s, f, defaultVal)
if isstruct(s) && isfield(s, f)
    v = s.(f);
else
    v = defaultVal;
end
end
