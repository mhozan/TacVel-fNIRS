function [rate_movavg_smooth, varargout] = peaks2rate(timeseries, varargin)
%PEAKS2RATE Calculate the rate of occurrences of peaks in a given time series
%
% USAGE:
%   [rate_movavg_smooth, tvec_out, peak_indices] = peaks2rate(timeseries, 'Param1', val1, 'Param2', val2, ...)
%
% INPUTS:
%   timeseries      -   A vector containing the time series data. Required.
%
%   'Fs'              -   The sampling frequency of the time series data.
%                       Optional. Default value is 1000.
%   f_lowpass           -   The cut-off frequency for the lowpass filter. 
%                       Optional. Default is no filtering. 5 Hz is a good
%                       input
%   'MinPeakHeight'     -   The minimum height of the peaks to be counted.
%                           Optional. Default value is 0.5.
%
%   'MinPeakProminence' -   The minimum prominence of the peaks to be counted.
%                           Optional. Default value is 2.
%
%   'MinPeakWidth'     -   The minimum width of the peaks to be counted.
%                           Optional. Default value is 1.
%
%   'MinPeakDistance'  -   The minimum distance between peaks to be counted.
%                           Optional. Default value is 1.
%   'MovingAverageDuration'          -   The duration of the moving average window in seconds.
%                           Optional. Default value is 5.
% OUTPUTS:
%   rate_movavg_smooth  -   The rate of occurrences of peaks with moving average smoothing applied.
%
%   tvec_out varargout{1} -   The time vector corresponding to the time series data.
%                           Optional. Only output if requested using nargout.
%
%   peak_indices varargout{2} -   The indices where the peaks occur.
%                                 Optional. Only output if requested using nargout.
%
%
% This function calculates the rate of occurrences of peaks in a given time series.
%
% The function first finds the peaks in the time series using the `findpeaks` function.
% The function then calculates the rate of occurrences of the peaks using the `mean` or `movmean` function.
%
% The function can be used to calculate the rate of peaks in a variety of physiological signals, such as ECG, EEG, EMG, respiration, and blood pressure.
%
% The following are some examples of physiological signals that have peaks and their rates matter:
%
% * **Electrocardiogram (ECG)**: The ECG measures the electrical activity of the heart. The peaks in an ECG correspond to the different phases of the heart's electrical cycle. The heart rate is the number of times the heart beats per minute. A normal resting heart rate for adults is between 60 and 100 beats per minute.
% * **Electroencephalogram (EEG)**: The EEG measures the electrical activity of the brain. The peaks in an EEG correspond to different brain waves, which are associated with different levels of consciousness, alertness, and sleep.
% * **Electromyogram (EMG)**: The EMG measures the electrical activity of muscles. The peaks in an EMG correspond to muscle contractions. The muscle contraction rate is the number of times a muscle contracts per minute. A normal muscle contraction rate for adults is between 5 and 10 contractions per minute.
% * **Respiration rate (RR)**: The RR measures the number of breaths per minute. A normal resting RR for adults is between 12 and 16 breaths per minute.
% * **Blood pressure (BP)**: The BP measures the pressure of the blood flowing through the arteries. The peaks in a BP waveform correspond to the systolic and diastolic pressures. The systolic pressure is the highest pressure in the arteries, and the diastolic pressure is the lowest pressure in the arteries. A normal resting BP for adults is less than 120/80 mmHg.
%
% Communication Neuroscience Laboratories,
% Center for Brain, Biology & Behavior
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2023

% Create input parser
p = inputParser;
addRequired(p, 'timeseries', @isvector);
addOptional(p, 'Fs', 1000, @(x) isscalar(x) && isnumeric(x));
addOptional(p, 'MovingAverageDuration', 5, @(x) isscalar(x) && isnumeric(x));
addOptional(p,'f_lowpass',5,@isscalar);
addParameter(p, 'MinPeakHeight', 0.5, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'MinPeakProminence', 2, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'MinPeakWidth', .001, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'MinPeakDistance', .005, @(x) isscalar(x) && isnumeric(x));
parse(p, timeseries, varargin{:});

% Get the parsed values
timeseries = p.Results.timeseries(:);
Fs = p.Results.Fs; %Hz
f_lowpass = p.Results.f_lowpass;
mov_avg_dur = p.Results.MovingAverageDuration; %second
minpeak_height = p.Results.MinPeakHeight;
minpeak_prominence = p.Results.MinPeakProminence;
minpeak_width = p.Results.MinPeakWidth;
minpeak_distance = p.Results.MinPeakDistance;

%lowpass filter
if f_lowpass > 0
    timeseries = lowpass(timeseries(:),f_lowpass,Fs,'ImpulseResponse','fir');

% % Define filter parameters
% order = 100;  % Butterworth filter order (adjust as needed)
% fc = f_lowpass / (Fs / 2);  % Normalize cutoff frequency
% 
% % Design Butterworth lowpass filter
% [b, a] = butter(order, fc, 'low');
% 
% % Apply zero-phase filtering using filtfilt
% timeseries = filtfilt(b, a, timeseries(:));


end


% Find indices and time instances of peaks

[~,peak_indices] = findpeaks(timeseries, ...
    'MinPeakHeight', minpeak_height, ...
    'MinPeakProminence',minpeak_prominence, ...
    'MinPeakWidth', minpeak_width, ...
    'MinPeakDistance', minpeak_distance);
[~,instances] = findpeaks(timeseries,Fs, ...
    'MinPeakHeight', minpeak_height, ...
    'MinPeakProminence',minpeak_prominence,    ...
    'MinPeakWidth', minpeak_width, ...
    'MinPeakDistance', minpeak_distance);


% Calculate the inter-peak interval in seconds
ipi = diff(instances);

% Calculate the rate per second
rate_avg = 1 / mean(ipi);

% Calculate the rate using a moving average
mov_avg_samples = floor(mov_avg_dur * Fs);
peak_tvec = zeros(size(timeseries));
peak_tvec(peak_indices) = 1;
rate_continuous = movsum(peak_tvec, mov_avg_samples) / mov_avg_dur;
rate_movavg_smooth = movmean(rate_continuous, mov_avg_samples);

% Create output time vector if requested
if nargout >= 2
    varargout{1} = (0:length(timeseries)-1)'/Fs;
end

% Output peak indices if requested
if nargout == 3
    varargout{2} = peak_indices;
end
