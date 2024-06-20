function [lineOut, fillOut] = TG_ShadeSEM_20210714(x, y, color, transparency)
% Plots mean as line and SEM as shaded area

meanData = mean(y, 1, 'omitnan'); % get mean over first dimension
stdData  = std(y, [], 1, 'omitnan'); % get standard deviation
semData  = stdData / sqrt(size(y, 1));
fillOut  = fill([x, fliplr(x)], [meanData + semData, fliplr(meanData - semData)], color, 'FaceAlpha', transparency, 'LineStyle', 'none');
hold on;
lineOut  = plot(x, meanData, 'color', color, 'LineWidth', 1, 'MarkerFaceColor', color, 'MarkerSize', 4);
% lineOut  = plot(x, meanData, 'color', color, 'LineWidth', 1);
end