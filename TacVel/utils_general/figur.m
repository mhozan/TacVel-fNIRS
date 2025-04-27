function [ varargout ] = figur( varargin )
%FIGUR creates a figure with set properties and optionally returns the handle.
%   Detailed explanation goes here
%   
%   Input Variables:
%
%   Output Variables:
%
%   Example(s):
%   
%   see also 
% FIGURE

% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2017

nargoutchk(0,1)
narginchk(0,1)

if nargin==1
    figname = varargin{1};
else
    figname = '-';
end
fignum = 1 + sum(figname - 0);

fh=figure(fignum); %clf; 
set(fh,'Name',figname,'numbertitle','off','menubar','figure')
ax=gca;
ax.NextPlot='add';
if nargout==1
    varargout{1} = fh;
end

