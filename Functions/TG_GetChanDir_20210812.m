function chanDir = TG_GetChanDir_20210812(paths, subject, session)
% TG_GetChanDir_20210812 sorts channel directory

chanDir     = dir(fullfile(paths, subject, session, 'chan*'));
tmp         = split(transpose({chanDir.name}), 'n');
[~, I]      = sort(cellfun(@str2num, tmp(:, 2)));
chanDir     = chanDir(I);
clear tmp;

end

