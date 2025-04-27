function [] = clfall()
%CLFALL clears all the open figures.
%   Detailed explanation goes here
%   
%   Input Variables:
%
%   Output Variables:
%
%   Example(s):
%   
%   see also 
% FIGURE, CLF

% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2017

fhz = get(0,'Children');
try
    for i=1:length(fhz)
        clf(fhz(i))
    end
catch
    
end
end


