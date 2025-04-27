function [ lpftimeseries ,varargout] = lowpass_accelerometer(timeseries,Fs,varargin)
% LOWPASS_ACCELEROMETER Returns a low-pass filtered version of the given 
% timeseries optimized for accelerometer analysis of body movements.
%
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2021
% See also
% fir_hpf_equiripple, firpm

f_lowpass = 20;
%In walking at natural velocity the bulk of acceleration power in the upper body ranges from 0.8–5 Hz, whereas the most abrupt accelerations occur at the foot in vertical direction during heel strike and sometimes amount up to 60 Hz [26]. By measuring these “worst case” accelerations at heel strike with a force platform Antonnson and Mann [27] demonstrated that 99% of the acceleration power during walking with bare feet is concentrated below 15 Hz. Higherfrequencies are caused by the impact between foot and walking surface and do not directly result from voluntary muscular work.
% C. V. C. Bouten, K. T. M. Koekkoek, M. Verduin, R. Kodde and J. D. Janssen, "A triaxial accelerometer and portable data processing unit for the assessment of daily physical activity," in IEEE Transactions on Biomedical Engineering, vol. 44, no. 3, pp. 136-147, March 1997, doi: 10.1109/10.554760.

if nargin==3
    f_lowpass = varargin{1};
end

b = fir_lpf_equiripple(Fs,f_lowpass); %Z:\students\mhozan2\Matlab Libraries\utils_general\fir_hpf_equiripple.m
% figure(444), clf, freqz(b,1,length(b)-1,Fs)

lpftimeseries = filtfilt(b,1,timeseries);

if nargout ==2
    varargout{1} = b;
end
    

