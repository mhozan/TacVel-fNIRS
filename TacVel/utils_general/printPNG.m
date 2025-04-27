function varargout = printPNG(fh,varargin)
%PrintALL saves a PNG version of the given figure.
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
% Mohsen Hozan hozan@huskers.unl.edu
% Date 2017
%% saving PNG files
%the first input is always the handle to the figure that will be saved in the current directory.
%the second input is optional, and is the path to where the figure will be saved. printPNG(fh,'F:\')
%the third output is optional, and is the name of the png file.: printPNG(fh,pathstr,'pngfilename(.png)')
nargoutchk(0,1)
narginchk(1,3)
pathstr = cd;
filename = ['ppo-',fh.Name,'-Figure ',num2str(fh.Number),datestr(now,'-yyyymmdd_HHMMSS')];
% filename = ['printpngoutput-',fh.Name,'-Figure ',num2str(fh.Number),datestr(now,'-yyyymmdd_HHMMSS')];

ext = '.png';

if nargin == 2
    pathstr = varargin{1};
    %     [pathstr,filename,ext] = fileparts(pathstring);
    if isempty(pathstr)
        pathstr = cd;
    elseif ~isdir(pathstr)
        mkdir(pathstr)
    end
elseif nargin == 3
    pathstr = varargin{1};
    if isempty(pathstr)
        pathstr = cd;
    elseif ~isdir(pathstr)
        mkdir(pathstr)
    end
    [~,filename,ext] = fileparts(varargin{2});
    if ~strcmpi(ext,'.png')
        filename = [filename,ext];
        ext = '.png';
    end
end

fullfilename = fullfile(pathstr,[filename,ext]);


% png_resolution  =	[3840 2160];
png_resolution  =	[3840 2160]./2;
% png_resolution  =	[2880 2160]./2; %4:3 Aspect Ratio
% png_resolution  =	[2000 2000]./2;

% png_resolution  =	[4000 3000]./2;
dpi = 100;
set(fh,'PaperUnits','inches','PaperPosition',[0 0 png_resolution(1)/dpi png_resolution(2)/dpi])

AdjustFontSize
% set(fh1,...
%     'PaperUnits','points',...
%     'PaperSize',[1800 1200])
% saveloc = 'F0comparisons\';
% saveloc = '';
% filename1=[saveloc,'test',num2str(obj.Dose),'-',datestr(now,'yymmdd')];
% filename1=['test-',datestr(now,'yymmdd')];

% print([filename1,'-fulllength.png'],'-dpng')
% print([filename1,'.svg'],'-dsvg','-cmyk')
print(fh,[fullfilename],['-r',num2str(dpi)],'-dpng','-cmyk')
disp(['PNG succesfully saved as :        ',fullfilename])


if nargout > 0
    disp(pathstr)
    varargout{1} = pathstr;    
end