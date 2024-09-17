function data = TG_taper(data, ratio)
%
% Input:
%   data    - signal to be tapered with welch window
%   ratio   - ratio of data's length to be tapered
%
% Output:
%   data    - tapered signal

% row vector:
if nargin ~= 2
    error('insert data and the ratio of data to be tapered');
end
if length(ratio) ~= 1 || ratio < 0 || ratio > 1  
    error('ratio has to be a scalar with value [0, 1]');
end
if size(data, 1) > size(data, 2)
    data = data';
end

% length of data
noOfSpData      = length(data);

% no of samples matching ratio
noOfSpToTaper   = round(noOfSpData * ratio);
N               = noOfSpToTaper - 1;

% welch - window
welch = 1 - (((0:N) - (N / 2)) / (N / 2)) .^ 2;

% mid-sample of taper range
m                           = floor(noOfSpToTaper / 2); % use floor to make sure it's an integer afterwards (to use as index)
data(1:m)                   = data(1:m) .* welch(1:m);
data((end - m + 1):end)     = data((end - m + 1):end).*welch((end - m + 1):end);