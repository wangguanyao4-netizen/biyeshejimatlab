function figureRoot = plot_original_plate_results(resultRoot)
%plot_original_plate_results Plot high-resolution original plate diagnostics.

if nargin < 1 || strlength(string(resultRoot)) == 0
    resultRoot = localLatestOriginalPlateRoot();
end
resultRoot = char(string(resultRoot));
summaryRoot = fullfile(resultRoot, "_summary");
if exist(fullfile(summaryRoot, "active_grid_summary.csv"), "file") ~= 2
    summarize_original_plate_results(resultRoot);
end

figureRoot = fullfile(summaryRoot, "figures");
dataRoot = fullfile(figureRoot, "figure_data");
dpimnumeric.ensureDir(figureRoot);
dpimnumeric.ensureDir(dataRoot);

E7 = readtable(fullfile(resultRoot, "E7_plate_SFEM", "summary.csv"));
E8 = readtable(fullfile(resultRoot, "E8_simultaneous_band", "summary.csv"));
active = readtable(fullfile(summaryRoot, "active_grid_summary.csv"));

localPlotAllGridCoverage(E7, E8, figureRoot, dataRoot);
localPlotFallbackInfRates(E7, E8, figureRoot, dataRoot);
localPlotActiveCoverage(active, figureRoot, dataRoot);
localPlotActiveCounts(active, figureRoot, dataRoot);

inventory = table( ...
    ["fig01_original_plate_allgrid_coverage_by_h"; ...
     "fig02_original_plate_fallback_inf_rates"; ...
     "fig03_original_plate_active_legacy_coverage"; ...
     "fig04_original_plate_active_grid_counts"], ...
    ["All-grid coverage by h"; ...
     "Fallback and Inf rates by h"; ...
     "Legacy active-grid coverage by h"; ...
     "Active grid counts under reference floors"], ...
    'VariableNames', {'figure_id', 'description'});
writetable(inventory, fullfile(figureRoot, "figure_inventory.csv"));
fprintf("Original plate figures written: %s\n", figureRoot);
end

function resultRoot = localLatestOriginalPlateRoot()
files = dir(fullfile(pwd, "results", "original_plate_full_*"));
files = files([files.isdir]);
if isempty(files)
    error("No results/original_plate_full_* directory found.");
end
[~, idx] = max([files.datenum]);
resultRoot = fullfile(files(idx).folder, files(idx).name);
end

function localPlotAllGridCoverage(E7, E8, figureRoot, dataRoot)
fig = localFigure();
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
localCoveragePanel(nexttile, E7, "E7 pointwise all-grid");
localCoveragePanel(nexttile, E8, "E8 pointwise and band all-grid");
localSaveFigure(fig, fullfile(figureRoot, "fig01_original_plate_allgrid_coverage_by_h"));
writetable([localFigureRows("E7", E7); localFigureRows("E8", E8)], ...
    fullfile(dataRoot, "fig01_original_plate_allgrid_coverage_by_h.csv"));
end

function localCoveragePanel(ax, tbl, titleText)
hold(ax, "on");
methods = ["Bootstrap-t", "Bootstrap-t fallback-rule", "Bootstrap-t simultaneous-band"];
colors = [0.121 0.466 0.705; 0.850 0.325 0.098; 0.000 0.500 0.300];
for i = 1:numel(methods)
    mask = string(tbl.method) == methods(i);
    if any(mask)
        plot(ax, tbl.h(mask), tbl.coverage(mask), "-o", "LineWidth", 1.6, ...
            "MarkerSize", 4, "Color", colors(i, :), "DisplayName", methods(i));
    end
end
yline(ax, 0.95, "k--", "Nominal 0.95", "LabelHorizontalAlignment", "left", ...
    "HandleVisibility", "off");
set(ax, "XScale", "log", "Box", "on", "FontName", "Times New Roman", "FontSize", 10);
xlabel(ax, "Smoothing coefficient h");
ylabel(ax, "Coverage");
title(ax, titleText, "FontWeight", "normal");
ylim(ax, [0, 1.05]);
grid(ax, "on");
legend(ax, "Location", "best", "Box", "off");
end

function localPlotFallbackInfRates(E7, E8, figureRoot, dataRoot)
fig = localFigure();
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
localRatePanel(nexttile, E7, "E7 fallback and Inf rates");
localRatePanel(nexttile, E8, "E8 fallback and Inf rates");
localSaveFigure(fig, fullfile(figureRoot, "fig02_original_plate_fallback_inf_rates"));
writetable([localFigureRows("E7", E7); localFigureRows("E8", E8)], ...
    fullfile(dataRoot, "fig02_original_plate_fallback_inf_rates.csv"));
end

function localRatePanel(ax, tbl, titleText)
hold(ax, "on");
fb = string(tbl.method) == "Bootstrap-t fallback-rule";
bt = string(tbl.method) == "Bootstrap-t";
plot(ax, tbl.h(fb), tbl.fallback_rate(fb), "-o", "LineWidth", 1.6, ...
    "MarkerSize", 4, "Color", [0.850 0.325 0.098], "DisplayName", "Fallback rate");
plot(ax, tbl.h(bt), tbl.bootstrap_t_inf_rate(bt), "-s", "LineWidth", 1.6, ...
    "MarkerSize", 4, "Color", [0.121 0.466 0.705], "DisplayName", "Raw BT Inf rate");
set(ax, "XScale", "log", "Box", "on", "FontName", "Times New Roman", "FontSize", 10);
xlabel(ax, "Smoothing coefficient h");
ylabel(ax, "Rate");
title(ax, titleText, "FontWeight", "normal");
ylim(ax, [0, 1.05]);
grid(ax, "on");
legend(ax, "Location", "best", "Box", "off");
end

function localPlotActiveCoverage(active, figureRoot, dataRoot)
legacy = active(string(active.active_definition) == "active_legacy", :);
fig = localFigure();
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
localActiveCoveragePanel(nexttile, legacy(string(legacy.experiment) == "E7_plate_SFEM", :), false, "E7 active-grid coverage");
localActiveCoveragePanel(nexttile, legacy(string(legacy.experiment) == "E8_simultaneous_band", :), true, "E8 active-grid coverage");
localSaveFigure(fig, fullfile(figureRoot, "fig03_original_plate_active_legacy_coverage"));
writetable(legacy, fullfile(dataRoot, "fig03_original_plate_active_legacy_coverage.csv"));
end

function localActiveCoveragePanel(ax, tbl, hasBand, titleText)
hold(ax, "on");
plot(ax, tbl.h, tbl.bootstrap_t_coverage, "-o", "LineWidth", 1.6, ...
    "MarkerSize", 4, "Color", [0.121 0.466 0.705], "DisplayName", "Raw BT");
plot(ax, tbl.h, tbl.bootstrap_t_fallback_rule_coverage, "-o", "LineWidth", 1.6, ...
    "MarkerSize", 4, "Color", [0.850 0.325 0.098], "DisplayName", "BT fallback");
if hasBand
    plot(ax, tbl.h, tbl.simultaneous_band_pointwise_coverage, "-^", "LineWidth", 1.6, ...
        "MarkerSize", 4, "Color", [0.000 0.500 0.300], "DisplayName", "Band pointwise");
end
yline(ax, 0.95, "k--", "Nominal 0.95", "LabelHorizontalAlignment", "left", ...
    "HandleVisibility", "off");
set(ax, "XScale", "log", "Box", "on", "FontName", "Times New Roman", "FontSize", 10);
xlabel(ax, "Smoothing coefficient h");
ylabel(ax, "Coverage on active grid");
title(ax, titleText, "FontWeight", "normal");
ylim(ax, [0, 1.05]);
grid(ax, "on");
legend(ax, "Location", "best", "Box", "off");
end

function localPlotActiveCounts(active, figureRoot, dataRoot)
fig = localFigure();
tiledlayout(1, 2, "TileSpacing", "compact", "Padding", "compact");
localActiveCountPanel(nexttile, active(string(active.experiment) == "E7_plate_SFEM", :), "E7 active grid count");
localActiveCountPanel(nexttile, active(string(active.experiment) == "E8_simultaneous_band", :), "E8 active grid count");
localSaveFigure(fig, fullfile(figureRoot, "fig04_original_plate_active_grid_counts"));
writetable(active, fullfile(dataRoot, "fig04_original_plate_active_grid_counts.csv"));
end

function localActiveCountPanel(ax, tbl, titleText)
hold(ax, "on");
defs = unique(string(tbl.active_definition), "stable");
colors = lines(numel(defs));
for iDef = 1:numel(defs)
    sub = tbl(string(tbl.active_definition) == defs(iDef), :);
    plot(ax, sub.h, sub.active_grid_count, "-o", "LineWidth", 1.4, ...
        "MarkerSize", 3.5, "Color", colors(iDef, :), "DisplayName", defs(iDef));
end
set(ax, "XScale", "log", "Box", "on", "FontName", "Times New Roman", "FontSize", 10);
xlabel(ax, "Smoothing coefficient h");
ylabel(ax, "Active grid count");
title(ax, titleText, "FontWeight", "normal");
grid(ax, "on");
legend(ax, "Location", "northwest", "Box", "off", "Interpreter", "none");
end

function rows = localFigureRows(experiment, tbl)
rows = tbl;
rows.experiment = repmat(string(experiment), height(rows), 1);
rows = movevars(rows, "experiment", "Before", 1);
end

function fig = localFigure()
fig = figure("Color", "w", "Units", "centimeters", "Position", [2, 2, 18, 8]);
end

function localSaveFigure(fig, basePath)
basePath = string(basePath);
exportgraphics(fig, basePath + ".png", "Resolution", 600);
exportgraphics(fig, basePath + ".tif", "Resolution", 600);
exportgraphics(fig, basePath + ".pdf", "ContentType", "vector");
close(fig);
end
