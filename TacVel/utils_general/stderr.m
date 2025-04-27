function [ standard_error ] = stderr(data)
%STDERR calcultes the standard error
%given in varargin.
%   Detailed explanation goes here
%
%   Input Variables:
%
%   Output Variables:
%
%   Example(s):
%
%   see also
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan mhozan2@unl.edu
% Date 2017
data = data(:);
data(isnan(data)) = []; %ignores nan values.
L = length(data);
standard_error = std(data)/sqrt(L);    
