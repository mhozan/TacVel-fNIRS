function flag=connected2web(varargin)
%CONNECTED2WEB checks if MATLAB has access to internet or not.
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2017
url2Bchecked = 'docs.google.com';
[~,stat]=dos(['ping -n 1 ',url2Bchecked]);
flag = contains(stat, 'Lost = 0');
