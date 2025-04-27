function [ hpftimeseries ,varargout] = highpass_emg(timeseries,fs,varargin)
%HIGHPASS_EMG Returns a high-passed filtered version of the given timeseries
% optimized for EMG analysis.
%
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2021
% See also
% fir_hpf_equiripple, firpm

f_hipass = 50;
if nargin==3
    f_hipass = varargin{1};
end

b = fir_hpf_equiripple(fs,f_hipass); %Z:\students\mhozan2\Matlab Libraries\utils_general\fir_hpf_equiripple.m
% figure(444), clf, freqz(b,1,length(b)-1,fs)

hpftimeseries = filtfilt(b,1,timeseries);

if nargout ==2
    varargout{1} = b;
end
    
% hpftimeseries = highpass(timeseries,f_hipass,fs,'ImpulseResponse','fir');

% hpftimeseries = movmean(abs_sig,winLen);
