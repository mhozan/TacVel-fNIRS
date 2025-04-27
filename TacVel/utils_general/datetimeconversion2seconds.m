function [ convmatrix6by1 ] = datetimeconversion2seconds(  )
%datetimeconversion2seconds can be multiplied by the output of datevec() to
%   
% Communication Neuroscience Laboratories,
% Center for Brain, Biology & Behavior
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2020
% Example:
%   dvec_seconds = datevec(now)*datetimeconversion2seconds;
%
%	see also 
%  datevec()
y = 365.24; %days in a year
convmatrix6by1 = [...
    y*24*60*60;...A year is 365.24 days long — that's why we have to skip a leap day every 100 years.
    (y/12)*24*60*60;...month
    24*60*60;...day
    60*60;...hour
    60;...minute
    1];

