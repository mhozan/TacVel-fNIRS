function [ mvgavgabs ] = movingavg_abs_emg(timeseries,fs,varargin)
%MOVINGAVG_ABS_EMG Returns a timeseries, which is the moving average of the
%absolute value of given timeseries over a 200ms windows, as is typical for
%EMG analysis.
%
%
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2021

abs_sig = abs(timeseries);

windur_s = 0.2;  %200 ms window duration
if nargin==3
    windur_s = varargin{1};
end

winLen = floor(windur_s*fs);
mvgavgabs = movmean(abs_sig,winLen);
