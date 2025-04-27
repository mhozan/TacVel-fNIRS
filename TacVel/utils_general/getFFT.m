function [fft_abs,fvec] = getFFT(timedomain_sig,fs,varargin)
%getDFT
%   Detailed explanation goes here
%
%   Input Variables:
%   timesomain_sig: input 1-D or 2-D timeseries in columnized format. 
%   fs: sampling frequency
%   f_bounds (3rd (optional) input): two-element numeric array: 
    %lower and upper bounds of the frequency band of interest.
%   default f_bounds is the full spectrum, i.e. [0,fs/2) Hz

%   Output Variables:
%   fft_abs : absolute value of the fft of the signal; mirrored values
%   across negative frequencies are exluded.
%   fvec : frequency vector corresponding to fft_abs in Hz based on f_bounds. 
%
%   Example(s):
%
%   see also
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan mhozan2@unl.edu
% Date 2018

narginchk(2,3)
if isvector(timedomain_sig)
    timedomain_sig = timedomain_sig(:); %columnize 1D inputs
end
f_bounds = [0 fs/2];
if nargin==3
f_bounds = varargin{1};
end

validateattributes(timedomain_sig,	{'numeric'},    {'2d'})
validateattributes(fs,              {'numeric'},    {'scalar','positive'})
validateattributes(f_bounds,       {'numeric'},    {'numel',2,'nonnegative','increasing','<=',fs/2}) 



nwind = size(timedomain_sig,1);
nfft = 2^(nextpow2(nwind)+2);
nfft2= nfft/2;
dim = 1; %fft across columns
fft_abs = abs(fft(timedomain_sig,nfft,dim))/nfft;
fft_abs = fft_abs(1:nfft2,:);
fvec_full = (fs/2) * (0:nfft2-1).'/(nfft2);

[~,fcut_indices]= min(abs(fvec_full-f_bounds));
selectedindices = fcut_indices(1):fcut_indices(2);

fvec = fvec_full(selectedindices);
fft_abs = fft_abs(selectedindices,:);


