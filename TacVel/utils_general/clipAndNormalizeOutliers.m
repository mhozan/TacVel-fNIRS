function timeseries = clipAndNormalizeOutliers(timeseries, varargin)
% clipAndNormalizeOutliers - Clips extreme outliers and normalizes the data
%
% This function clips values in the input timeseries data that fall outside the
% specified percentile range, setting them to the nearest boundary. Then, it
% normalizes the data to a range of [0, 1], where 0 corresponds to the clipped
% lower boundary and 1 corresponds to the maximum value after clipping.
%
% Usage:
%   timeseries = clipAndNormalizeOutliers(timeseries, 'lowerprctile', 90, 'multiplyfactor', 10);
%
% INPUT:
%   timeseries (n x 1 vector or 1 x n row vector) - A vector containing the time-series data.
%   lowerprctile (optional)    - The lower percentile to use for clipping (default: 90).
%   multiplyfactor (optional)   - The factor by which the lower percentile is multiplied
%                                 to get the upper threshold (default: 10).
%
% OUTPUT:
%   timeseries (n x 1 vector)  - The input time-series data with extreme outliers clipped and normalized.
%
% Communication Neuroscience Laboratories,
% Center for Brain, Biology & Behavior
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2025

% Parse input arguments
p = inputParser;
addRequired(p, 'timeseries');
addParameter(p, 'lowerprctile', 97, @(x) isnumeric(x) && x >= 0 && x <= 100);
addParameter(p, 'multiplyfactor', 2, @(x) isnumeric(x) && x > 0);
parse(p, timeseries, varargin{:});

timeseries = timeseries(:);  % Convert to column vector if it's a row

% Get parsed values
lowerprctile = p.Results.lowerprctile;
multiplyfactor = p.Results.multiplyfactor;

% Calculate thresholds
threshlow = prctile(timeseries, lowerprctile);
threshhi = threshlow * multiplyfactor;

% Clip values outside the threshold
timeseries(timeseries > threshhi) = threshhi;
timeseries(timeseries < threshlow) = threshlow;

% Normalize and handle NaNs
timeseries = timeseries - threshlow;
timeseries = timeseries ./ max(timeseries);
timeseries(isnan(timeseries)) = 1; % Handle NaNs by setting them to 1

