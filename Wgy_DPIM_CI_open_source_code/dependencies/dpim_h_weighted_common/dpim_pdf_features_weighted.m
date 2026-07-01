function feat = dpim_pdf_features_weighted(yGrid, pdf)
%DPIM_PDF_FEATURES_WEIGHTED Basic density-shape diagnostics.
yGrid = yGrid(:); pdf = pdf(:);
pdf = max(pdf, 0);
area = trapz(yGrid, pdf);
if area > 0; pdf = pdf / area; end
[peakVal, peakIdx] = max(pdf);
half = peakVal/2;
leftIdx = find(pdf(1:peakIdx) <= half, 1, 'last');
rightIdx = find(pdf(peakIdx:end) <= half, 1, 'first');
if isempty(leftIdx); leftIdx = 1; end
if isempty(rightIdx); rightIdx = numel(pdf); else; rightIdx = rightIdx + peakIdx - 1; end
d = diff(pdf);
nModes = sum(d(1:end-1) > 0 & d(2:end) < 0);
if peakVal > 0; nModes = max(nModes, 1); end
feat = struct();
feat.area = area;
feat.peak_y = yGrid(peakIdx);
feat.peak_val = peakVal;
feat.FWHM = yGrid(rightIdx) - yGrid(leftIdx);
feat.n_modes = nModes;
end
