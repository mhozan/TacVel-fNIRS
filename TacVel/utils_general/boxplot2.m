function [varargout]= boxplot2(caov, varargin)
% input data in form of cell array of r x c dimensions where r = number of
% groups (such as young/old mice - grouped along the x axis) and c is the 
% secondary group condition (such as doses). Each cell sould contain a 
% vector of N sample observations
% Optional inputs include:
%   Figure Title:  string
%   Legend labels:  cell of strings
%   Legend Title:   string
%   Y axis Labels:  cell of strings
%% Varargin Parsing
%   
p = inputParser;
paramnames = {...
   'title';...
   'legendlabels';...
   'legendtitle';...
   'ylabel';...
   'xlabel';...
   'facealpha';...
   'alpha'};
defaultz = {...
    ' ';...
    {'50','100','200','400'};...
    {'Dose (ug/kg)'};...
    {'Time (min)'};...
    {'Young','Old','All'};...
    .25;...
    0.05};
for pn = 1:numel(paramnames)
    addParameter(p,paramnames{pn},defaultz{pn})
end    
parse(p,varargin{1:end})

figtitle = p.Results.title;
lgdlabels = p.Results.legendlabels;
ldgtitle = p.Results.legendtitle;
ylabelstr = p.Results.ylabel;
xlabelstr = p.Results.xlabel;
facealpha = p.Results.facealpha;
alpha = p.Results.alpha;
%%
loopcounter = reshape(1:size(caov,1)*size(caov,2),size(caov,2),size(caov,1))';
colorz = distinguishable_colors(size(caov,2)+1);
% colorz =    [1 .9 .9;
%             1 .6 .6;
%             1 .3 .3;
%             1 .1 .1];

barW = 1;
endbarW = 0.75;
barB = 0.25;

segment = (barW+barB)*(1:size(caov,2));
scaler  = (2+size(caov,2))*[0:size(caov,1)-1]';
posmat = repmat(segment,3,1)+scaler;
hfig=figure;clf;
ax=gca;
for r=1:size(caov,1)
  for c=1:size(caov,2)
      lc = loopcounter(r,c);
      v=[caov{r,c}];      
%       v=[caov{r,c}]./3600;      
      mx=max(v);
      mn=min(v);
      ci = bootci2(v,alpha);
      xcent = posmat(r,c);
      x = [xcent-(barW/2),xcent-(barW/2),xcent+(barW/2),xcent+(barW/2)];
      y = [ci(1),ci(end),ci(end),ci(1)];
      mbarx = [xcent-(barW/2),xcent+(barW/2)];
      mbary = [ci(2),ci(2)];
      
      hs=scatter(ax,ones(length(v),1)*xcent,v);
      set(hs,'MarkerFaceColor', colorz(c+1,:),'MarkerEdgeColor','black');
      hold on
      hb(lc)=fill(ax,x,y,'k');
      set(hb(lc),'FaceColor',colorz(c+1,:),'FaceAlpha',facealpha,...
          'EdgeColor',colorz(c+1,:),'EdgeAlpha',1);
     
      hl=line(ax,mbarx,mbary,'color','r');
      hw(1) = line(ax,[xcent,xcent],[ci(end),mx]);
      hw(2) = line(ax,[xcent,xcent],[ci(1),mn]);
      set(hw,'LineStyle','--','Color','black');
      h_endbar(1) = line(ax,[xcent-(endbarW/2),xcent+(endbarW/2)],[mx,mx]);
      h_endbar(2) = line(ax,[xcent-(endbarW/2),xcent+(endbarW/2)],[mn,mn]);
      set(h_endbar,'color','black')
      
  end
end
ax.XTick = 3:6:18;
ax.XTickLabel = xlabelstr;
ax.FontSize = 12;
ax.FontWeight = 'bold';
ylabel(ax,ylabelstr)
uistack(hb,'top')
lgd=legend(hb(1:4),lgdlabels);
title(lgd,ldgtitle)
legend boxoff
title(ax,figtitle);
varargout = {ax};
