function varargout=colorplot(varargin)
% COLORPLOT will function simlar to Matlab's built-in "plot" function, except uses a default
% colormap, in order to produce a gradient-colored 2D plot. see the link
% below for more information:
% https://stackoverflow.com/questions/42174917/2-d-line-gradient-color-in-matlab
% 
% Example(s):
% colorplot(1:1000,sin(.01:0.01:10))
% 
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan mhozan2@unl.edu
% Date Feb 2018
if isnumeric(varargin{1})
    lngth = length(varargin{1});
else  %in case the first input is the 'axes' for the plot.
    lngth = length(varargin{2});
end
    
clrmap = jet(lngth)*255;
% clrmap = hsv(lngth)*255;
% clrmap = parula(lngth)*255;



% cmap = colormap(clrmap)*255;
cmap = uint8([clrmap,   ones(lngth,1)].');
p = plot(varargin{1:end});
drawnow
set(p.Edge, 'ColorBinding','interpolated', 'ColorData',cmap)
drawnow

if nargout==1
varargout{1} = p;
end