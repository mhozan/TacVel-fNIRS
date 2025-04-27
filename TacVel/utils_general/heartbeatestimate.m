function [ hrate_bpm ,varargout] = heartbeatestimate(ECG_vec,varargin)
% HEARTBEATESTIMATE estimates the heart rate in beats per minute and
% returns  it in a timeseries of the same size as the input.
% It's basically a wrapper for peaks2rate function
% with specialized parameters and constraints to detect the heart rate
%
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2023
% See also
% nirs.muscle.estimate_heartrate, peaks2rate

p = inputParser;
addRequired(p, 'ECG_vec', @isvector);
addOptional(p, 'Fs', 1000, @(x) isscalar(x) && isnumeric(x));
addOptional(p, 'MovingAverageDuration', 5, @(x) isscalar(x) && isnumeric(x));
addOptional(p,'f_lowpass',5,@isscalar);
addParameter(p, 'MinPeakHeight', .1, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'MinPeakProminence', .3, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'MinPeakWidth', .02, @(x) isscalar(x) && isnumeric(x));
addParameter(p, 'MinPeakDistance', .43, @(x) isscalar(x) && isnumeric(x));
parse(p, ECG_vec, varargin{:});

% Get the parsed values
Heart_vec = p.Results.ECG_vec(:);
Fs = p.Results.Fs; %Hz
f_lowpass = p.Results.f_lowpass;
mov_avg_dur = p.Results.MovingAverageDuration; %second
minpeak_height = p.Results.MinPeakHeight;
minpeak_prominence = p.Results.MinPeakProminence;
minpeak_width = Fs*p.Results.MinPeakWidth;
minpeak_distance = Fs*p.Results.MinPeakDistance;

% Highpass to get rid of DC drifts.
% figure(23123), clf, plot(Heart_vec), hold on
f_hipass = 2; %Hz
% Heart_vec = highpass(Heart_vec,f_hipass,Fs);
Heart_vec = (highpass(Heart_vec,f_hipass,Fs,'ImpulseResponse','fir'));
% plot(eart_vec)


[Heart_rate,tvec,heartind] = peaks2rate(Heart_vec,...
    'Fs',Fs,...
    'f_lowpass',f_lowpass, ...
    'MinPeakProminence',minpeak_prominence, ...
    'MinPeakDistance',minpeak_distance, ...
    'MinPeakHeight',minpeak_height, ...
    'MinPeakWidth',minpeak_width, ...
    'MovingAverageDuration',mov_avg_dur);


hrate_bpm = 60*Heart_rate;
upperbound=140;
lowerbound=20;
hrate_bpm(hrate_bpm>upperbound) = upperbound;
hrate_bpm(hrate_bpm<lowerbound) = lowerbound;

% figure(11), clf; plot(obj.time,obj.data), hold on, scatter(y2(y3),ones(size(y3))),plot(obj.time,ata/60)

% Create output time vector if requested
if nargout >= 2
    varargout{1} = tvec;
end

% Output peak indices if requested
if nargout == 3
    varargout{2} =heartind;
end