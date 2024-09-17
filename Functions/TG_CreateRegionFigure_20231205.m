function fig = TG_CreateRegionFigure_20231205(regions, condBl, condEnc, condRec, incIdx, colorData)
% TG_CreateRegionFigure_20231205 Creates figures for the script
% TG_FigurePhaseLockingBrainRegion_20231205

% create figure
fig             = figure;

% results
res             = [condBl(incIdx, 1)'; condEnc(incIdx, 1)'; condRec(incIdx, 1)'];

% error bars
if size(condBl, 2) == 2 % SEM
    resErr      = [condBl(incIdx, 2)'; condEnc(incIdx, 2)'; condRec(incIdx, 2)'];
end

% bar plot
b               = bar(regions, res, 'FaceColor', 'flat');
hold on;

% calculate the number of conditions and regions
[nCond, nReg]   = size(res);

% loop through conditions
xBar            = nan(nCond, nReg);
for iCond = 1:nCond

    % get the x coordinate of the bars
    xBar(iCond, :) = b(iCond).XEndPoints;
    
    % adjust colors
    b(iCond).CData = colorData(iCond, :);
end

% plot error bars
if size(condBl, 2) == 2 % SEM
    errorbar(xBar, res, resErr, 'LineStyle', 'none', 'Color', 'k');
end

% settings
set(gca, 'tickDir', 'out');
box off;

end

