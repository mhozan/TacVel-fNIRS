function printSVG(fh,varargin)
narginchk(1,2)
if nargin==2
    filestring=varargin{1};
    [pathstr,filename,ext]=fileparts(filestring);
    if ~isdir(pathstr)
        mkdir(pathstr)
    end
    filenamefinal = fullfile(pathstr,filename);
    
    
else
    filenamefinal = ['test-',datestr(now,'yymmdd')];
end

png_resolution  =	[3840 2160]./2;
png_resolution  =	[2880 2160]./2; %4:3 Aspect Ratio

fntsize = 10;
set(fh,'PaperUnits','inches','PaperPosition',[0 0 png_resolution(1)/150 png_resolution(2)/150])
fontsizeGroup=findall(gcf,'-property','FontSize');
for ii=1:length(fontsizeGroup)
    if fontsizeGroup(ii).FontSize<fntsize
    fontsizeGroup(ii).FontSize=fntsize;
    end
end
set(findall(gcf,'-property','FontName'),'FontName','Arial')
% set(fh1,...
%     'PaperUnits','points',...
%     'PaperSize',[1800 1200])
% saveloc = 'F0comparisons\';
% saveloc = '';
% filename1=[saveloc,'test',num2str(obj.Dose),'-',datestr(now,'yymmdd')];
% filename1=['test-',datestr(now,'yymmdd')];

% print([filename1,'-fulllength.png'],'-dpng')
% print([filename1,'.svg'],'-dsvg','-cmyk')
print([filenamefinal,'.svg'],'-r150','-dsvg','-cmyk')


end