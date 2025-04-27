function [CInMEAN] = bootci2(inputmatrix,varargin)
% see also BOOTCI
% mohsen 10/28/2016

if nargin >= 2
    alpha = varargin{1};      
    if isempty(alpha)
        alpha = 0.05;
    end
else
    alpha = 0.05;
end
if nargin >=3   
    bootfun = varargin{2};
   	if isempty(bootfun)
        bootfun = 'mean';
    end
else
    bootfun = 'mean'; %'median';
end

nboot = 10000;
UseParallel     =   false;
Options = statset('UseParallel',UseParallel);


prctilez = [(alpha/2)*100,50,(1-alpha/2)*100];

[bootstat,~] =   bootstrp(nboot, bootfun,  inputmatrix,	'Options',Options);
CInMEAN      =   prctile(bootstat,prctilez);