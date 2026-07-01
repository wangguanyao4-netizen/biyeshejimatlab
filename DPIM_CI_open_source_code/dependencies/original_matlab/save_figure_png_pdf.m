function save_figure_png_pdf(fig, filePathNoExt, resolution)
%save_figure_png_pdf  Save figure as PNG (and PDF if supported)
if nargin < 3 || isempty(resolution); resolution = 180; end
pngPath = char(filePathNoExt);
exportgraphics(fig, [pngPath, '.png'], 'Resolution', resolution);
try
    exportgraphics(fig, [pngPath, '.pdf'], 'ContentType', 'vector');
catch
end
end
