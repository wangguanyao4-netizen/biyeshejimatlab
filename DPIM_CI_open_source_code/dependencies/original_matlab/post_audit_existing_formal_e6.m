function post_audit_existing_formal_e6(resultsRoot)
%POST_AUDIT_EXISTING_FORMAL_E6 Downgrade weak E6 effective-order overclaim.
%
% This is a lightweight re-audit for an already completed formal campaign.
% It does not rerun numerical experiments; it only updates audit artifacts
% when the scrambled RQMC effective-order slope is too weak for a positive
% main-text efficiency claim.

if nargin < 1 || strlength(string(resultsRoot)) == 0
    error("resultsRoot is required.");
end
resultsRoot = string(resultsRoot);
methodPath = fullfile(resultsRoot, "method_audit.csv");
claimPath = fullfile(resultsRoot, "claim_audit.csv");
orderPath = fullfile(resultsRoot, "formal_experiments", ...
    "E6_rqmc_effective_order", "effective_order_summary.csv");

methodAudit = readtable(methodPath, TextType="string", Delimiter=",");
claimAudit = readtable(claimPath, TextType="string", Delimiter=",");
orderSummary = readtable(orderPath, TextType="string", Delimiter=",");
rqmc = orderSummary(orderSummary.method == "sobol_scrambled", :);
if isempty(rqmc)
    error("No sobol_scrambled row in %s.", orderPath);
end

ord = rqmc.effective_order_variance(1);
ratio = rqmc.variance_ratio_at_max_n_vs_mc(1);
oldMethod = methodAudit.method_status(methodAudit.paper_experiment == "E6" ...
    & methodAudit.method == "sobol_scrambled");
oldClaim = claimAudit.claim_status(claimAudit.paper_experiment == "E6");

methodMask = methodAudit.paper_experiment == "E6" & methodAudit.method == "sobol_scrambled";
methodAudit.method_status(methodMask) = "rqmc_constructed_weak_order";
methodAudit.note(methodMask) = sprintf(['n=%g..%g; variance order=%.4g; ' ...
    'variance ratio at max n=%.4g; RQMC construction and weights are valid, ' ...
    'but effective-order gain is weak; use as supplement/diagnostic, not a positive speedup claim.'], ...
    rqmc.n_min(1), rqmc.n_max(1), ord, ratio);

claimMask = claimAudit.paper_experiment == "E6";
claimAudit.claim_status(claimMask) = "supplement_only";
claimAudit.main_text_decision(claimMask) = "rqmc_constructed_but_effective_order_weak";
claimAudit.note(claimMask) = sprintf(['Scrambled RQMC point pools and Voronoi weights are constructed ' ...
    'without fallback, but formal effective-order evidence is weak: variance order=%.4g, ' ...
    'variance ratio at max n=%.4g. Keep this as construction/weight-moment diagnostic or supplement, ' ...
    'not as a main positive efficiency claim.'], ord, ratio);

writetable(methodAudit, methodPath);
writetable(claimAudit, claimPath);
statusSummary = groupsummary(claimAudit, "claim_status");
writetable(statusSummary, fullfile(resultsRoot, "_formal_report", "claim_status_summary.csv"));

correction = table("E6", oldClaim(1), "supplement_only", oldMethod(1), ...
    "rqmc_constructed_weak_order", ord, ratio, ...
    "Post-audit correction: weak formal RQMC effective-order gain should not be overclaimed.", ...
    'VariableNames', {'paper_experiment','old_claim_status','new_claim_status', ...
    'old_method_status','new_method_status','effective_order_variance', ...
    'variance_ratio_at_max_n_vs_mc','note'});
writetable(correction, fullfile(resultsRoot, "post_audit_e6_correction.csv"));

claimTexPath = fullfile(resultsRoot, "_formal_report", "tables", "formal_claim_audit_table.tex");
if isfile(claimTexPath)
    tex = string(fileread(claimTexPath));
    tex = replace(tex, ...
        "E6 & rqmc\_effective\_order\_diagnostic & main\_text\_ready & rqmc\_effective\_order\_documented \\", ...
        "E6 & rqmc\_effective\_order\_diagnostic & supplement\_only & rqmc\_constructed\_but\_effective\_order\_weak \\");
    dpimnumeric.writeText(claimTexPath, tex);
end

post = sprintf(['# Post-audit E6 Correction\n\n' ...
    '- Applied to: `%s`\n' ...
    '- Reason: formal scrambled RQMC variance order is only `%.6g` and the variance ratio at max n is `%.6g`; this documents construction but not a convincing effective-order gain.\n' ...
    '- Updated: `claim_audit.csv`, `method_audit.csv`, `_formal_report/claim_status_summary.csv`, and `formal_claim_audit_table.tex`.\n' ...
    '- Paper implication: E6 can support the statement that sample pools use randomized RQMC and probability weights; do not claim a strong RQMC acceleration/effective-order theorem from this run.\n'], ...
    resultsRoot, ord, ratio);
postPath = fullfile(resultsRoot, "post_audit_e6_correction.md");
dpimnumeric.writeText(postPath, post);

appendTargets = [fullfile(resultsRoot, "paper_update_notes.md"), ...
    fullfile(resultsRoot, "_formal_report", "weighted_paper_formal_campaign_report.md")];
for i = 1:numel(appendTargets)
    p = appendTargets(i);
    if isfile(p)
        txt = string(fileread(p));
        dpimnumeric.writeText(p, txt + newline + newline + post);
    end
end

htmlPath = fullfile(resultsRoot, "_formal_report", "weighted_paper_formal_campaign_report.html");
if isfile(htmlPath)
    html = string(fileread(htmlPath));
    banner = sprintf(['<h2>Post-audit E6 Correction</h2><p>Scrambled RQMC construction passed, ' ...
        'but effective-order gain is weak: variance order %.6g, max-n variance ratio %.6g. ' ...
        'Treat E6 as supplement/diagnostic, not a main speedup claim.</p>'], ord, ratio);
    html = replace(html, "</body></html>", banner + "</body></html>");
    dpimnumeric.writeText(htmlPath, html);
end

fprintf("post_audit_done order=%.6g ratio=%.6g\n", ord, ratio);
end
