function dpim_save_current_figure(fig, pathBase)
%DPIM_SAVE_CURRENT_FIGURE Save current figure as PNG and PDF when possible.
[folder,~,~] = fileparts(pathBase);
if ~exist(folder, 'dir'); mkdir(folder); end
try
    exportgraphics(fig, [pathBase '.png'], 'Resolution', 600);
catch
    saveas(fig, [pathBase '.png']);
end
try
    exportgraphics(fig, [pathBase '.pdf'], 'ContentType', 'vector');
catch
    saveas(fig, [pathBase '.pdf']);
end
end
