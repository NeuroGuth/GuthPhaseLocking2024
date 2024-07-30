function [xgp, wt] = TG_generalized_phase_vector(x, Fs, lp)
%
% TG_generalized_phase_vector calculates the generalized phase of an input vector
%
% INPUT:
%   x       - data column vector (t, 1)
%   Fs      - sampling rate (Hz)
%   lp      - low-frequency data cutoff (Hz)
%
% OUTPUT:
%   xgp     - output datacube
%   wt      - instantaneous frequency estimate
%
% Adapted from Davis et al., Nature, 2020.
% https://github.com/mullerlab/generalized-phase
%
% Tim Guth, 2023

% parameter
nwin            = 3; % extension factor for identified negative frequency epochs

% handle input
assert(iscolumn(x), 'Column vector input required');

% number of samples
npts            = length(x);

% time difference between samples
dt              = 1 / Fs;

% function for rewrapping phase
rewrap          = @(xp) (xp - 2 * pi * floor((xp - pi) / (2 * pi)) - 2 * pi);

% function for piecewise cubic Hermite interpolation
naninterp       = @(xp) interp1(find(~isnan(xp)), xp(~isnan(xp)), find(isnan(xp)), 'pchip');

% analytic signal representation (single-sided Fourier approach, cf. Marple 1999)
xo              = fft(x, npts, 1);
h               = zeros(npts, 1);
if (npts > 0) && (mod(npts, 2) == 0)
    h([1, npts / 2 + 1])    = 1;
    h(2:npts / 2)           = 2;
else
    h(1)                    = 1;
    h(2:(npts + 1)/2)       = 2;
end
xo              = ifft(xo .* h(:, ones(1, size(x, 2))), [], 1);
ph              = angle(xo);
md              = abs(xo);

% calculate instantaneous frequency
wt              = zeros(size(xo));
wt(1:end - 1)   = angle(xo(2:end) .* conj(xo(1:end - 1))) ./ (2 * pi * dt);

% account for sign of instantaneous frequency
sign_if         = sign(mean(wt(:), 'omitnan'));
if sign_if == -1
    modulus         = abs(xo);
    ang             = sign_if .* angle(xo); % rectify rotation
    xo              = modulus .* exp(1i .* ang);
    ph              = angle(xo);
    md              = abs(xo);
    wt(1:end - 1)   = angle(xo(2:end) .* conj(xo(1:end - 1))) ./ (2 * pi * dt); % re-calculate IF
end

% check if all phase values are NaN
if all(isnan(ph))
    xgp             = nan(size(xo));
    return;
end

% find negative frequency epochs (i.e. less than LP cutoff)
idx             = (wt < lp);
idx(1)          = 0;
[L, G]          = bwlabel(idx);
for kk = 1:G
    idxs            = find(L == kk);
    upper_limit     = idxs(1) + ((idxs(end) - idxs(1)) * nwin);
    if upper_limit > npts
        idx(idxs(1):npts)           = true;
    else
        idx(idxs(1):upper_limit)    = true;
    end
end

% create a temporary variable 'p' to handle phase adjustments
p               = ph;

% unwrap the phase to ensure continuity
p               = unwrap(p);

% set phases corresponding to negative frequency epochs to NaN
p(idx)          = NaN;

% check if all phases are NaN after adjustment
if all(isnan(p))
    xgp             = nan(size(xo));
    return;
end

% interpolate NaN values in the phase
p(isnan(p))     = naninterp(p);
p               = rewrap(p);
ph              = p;

% output
xgp             = md .* exp(1i .* ph);
