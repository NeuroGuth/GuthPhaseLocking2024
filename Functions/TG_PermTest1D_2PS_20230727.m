function [rank, t, tSur] = TG_PermTest1D_2PS_20230727(x, y, nSur, randSeedNum)

% TG_PermTest1D_2PS_20230727 performs a permutation test for two paired
% samples using t-values as the test statistic. The observations are
% randomly flipped or kept between the two different samples.
%
% Input:
% x             first data vector
% y             second data vector
% nSur          number of permutations
% randSeedNum   random seed numbers
%
% Output:
% rank     rank of empirical t-value in surrogate t-value distribution
% t        empirical t-value
% tSur     distribution of surrogate t-values
%
% Tim Guth, 2023

% transpose x, if necessary
if size(x, 2) > size(x, 1)
    x = x';
end

% transpose y, if necessary
if size(y, 2) > size(y, 1)
    y = y';
end

% concatenate x and y
xy = [x, y];

% calculate empirical t-value
[~, ~, ~, stats] = ttest(xy(:, 1), xy(:, 2));
t                = stats.tstat;

% prealloacte
tSur = nan(nSur, 1);

% loop through surrogates
parfor iSur = 1:nSur
    
    % display progress
    disp(iSur);

    % set random seed
    rng(randSeedNum(iSur, 1), 'twister');

    % create random logical for label swapping
    bRand               = logical(randi([0, 1], size(xy, 1), 1));
    
    % randomly swap labels
    surXy               = xy;
    surXy(bRand, :)     = flip(surXy(bRand, :), 2);
    
    % calculate surrogate t-value
    [~, ~, ~, surStats] = ttest(surXy(:, 1), surXy(:, 2));
    tSur(iSur, 1)       = surStats.tstat;
end

% calculate rank
rank = sum(t > tSur) / nSur;

end
