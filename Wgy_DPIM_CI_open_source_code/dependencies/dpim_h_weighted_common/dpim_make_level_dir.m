function levelRoot = dpim_make_level_dir(cfg, levelName)
levelRoot = fullfile(cfg.results_root, char(levelName));
if ~exist(levelRoot, 'dir'); mkdir(levelRoot); end
if ~exist(fullfile(levelRoot, 'data'), 'dir'); mkdir(fullfile(levelRoot, 'data')); end
if ~exist(fullfile(levelRoot, 'figures'), 'dir'); mkdir(fullfile(levelRoot, 'figures')); end
end
