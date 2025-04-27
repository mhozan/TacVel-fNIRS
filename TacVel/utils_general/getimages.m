function [fullpathnamez]=getimages(varargin)
if nargin == 0
    startdir = '\\secd.unl.edu\barlow\students\mhozan2\Matlab Scripts\plot3by1_tiffz\';
    [FileName,PathName,FilterIndex] = uigetfile('*.tiff','Select the images with different colormaps','MultiSelect', 'on',startdir);
elseif nargin==1
   PathName = varargin{1}; 
   FileName = dir([PathName,filesep,'*.tiff']);
   FileName = {FileName.name};
   excludeFINALtaggedTIFF = contains(FileName,'final','IgnoreCase',true);
   FileName = FileName(~excludeFINALtaggedTIFF);
end

if isequal(FileName,0)
   disp('User selected Cancel')
   return
else
   fullpathnamez = fullfile(PathName, FileName);
end