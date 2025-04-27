function data = applyBandPassFilter(data, varargin)
% applyBandPassFilter - FIR band/low/high-pass filter for time-series (cols)
%
%% Usage:
%   data = applyBandPassFilter(data, 'lowpass', 30, 'highpass', 1);
%
% INPUT:
%   data: matrix (time x channels)
%% OPTIONS:
%   'Fs'       - Sampling rate (default: 1000)
%   'lowpass'  - Upper cutoff frequency (Hz)
%   'highpass' - Lower cutoff frequency (Hz)
%   'filter_order' - FIR filter order (default: 500)
%% EXAMPLES;
% Bandpass between 1 and 30 Hz
% filtered = applyBandPassFilter(rawData, 'lowpass', 30, 'highpass', 1);
% Just lowpass at 10 Hz
% filtered = applyBandPassFilter(rawData, 'lowpass', 10);
% See also nirs.cnl.BandPassFilter
% Communication Neuroscience Laboratories,
% Center for Brain, Biology & Behavior
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2025

p = inputParser;
addParameter(p, 'Fs', 1000);
addParameter(p, 'lowpass', []);
addParameter(p, 'highpass', []);
addParameter(p, 'filter_order', 1023);
parse(p, varargin{:});
Fs = p.Results.Fs;
lp = p.Results.lowpass;
hp = p.Results.highpass;
N = p.Results.filter_order;

if isempty(lp) && isempty(hp)
    disp('No action performed. Provide lowpass or highpass (or both) cutoff(s).')
    return;
end

% Normalized cutoff
Wn = [];
if ~isempty(hp)
    Wn(1) = hp / (Fs/2);
end
if ~isempty(lp)
    Wn(end+1) = lp / (Fs/2);
end

% Filter design
if length(Wn) == 1
    mode = 'low';
    if ~isempty(hp)
        mode = 'high';
    end
else
    mode = 'bandpass';
end
b = fir1(N, Wn, mode);
% figur(mfilename),clf, freqz(b,1,N)
% Apply filter
data = filtfilt(b, 1, data);


