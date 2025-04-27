function [bootstat,varargout]=bootstrp2(bootstrapmatrix,varargin)
%BOOTSTRP2 uses built-in BOOTSTRP function to calculate the overal
%empirical distribution confidence intervals based on the mean-value
% See also
% BOOTSTRP
%mohsen hozan 10/22/2016
narginchk(1,2)
nargoutchk(1,2)
validateattributes(bootstrapmatrix,{'numeric'},{'nonempty'})

nboot = 1000;
UseParallel     =   false;
Options = statset('UseParallel',UseParallel);
if nargin==2
    validateattributes(varargin{1},{'numeric'},{'nonempty','scalar','integer','>=',100,'<=',1e8})
    nboot = varargin{1};
end
[bootstat,~]         	=   bootstrp(nboot, @mean,  bootstrapmatrix,	'Options',Options);

% prctilez = [0.5,50,99.5];    %	99% CI and mean
prctilez = [2.5,50,97.5];               %   95% CI and Mean

if nargout==2
    varargout{1} = prctile(bootstat,prctilez); %95 percent confidence intervals (two sided) and the mean
end
end