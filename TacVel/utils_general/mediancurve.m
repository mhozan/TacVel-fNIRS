function [y,varargout] = mediancurve(X,varargin)

%MEDIANCURVE finds the one vector among the given vectors, which is defined
%as the least deviated curve from the mean of all the input curves.
%
% Syntax
% y = mediancurve(X)
% y = mediancurve(X,dim)
% [y,ind] = mediancurve(X,dim);


%   Input Variables:
%   X: An [MxN] matrix, where N is the number of curves (vectors) and M is the
%   length of each curve.
%   dim: dimension for calculation. 1: across columns, 2[default]: across rows(curves).
%
%   Output Variables:
%   y: the one curve selected out of N given curves, whose least deviated
%   from the mean of all N curves by a least-square criteria.
%   ind: index of y, a number in the range [1,N]
%
%
% see also
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan mhozan2@unl.edu
% Date 2019
narginchk(1,2)
nargoutchk(1,2)

dim = 2;
if nargin >= 2
	dim = varargin{1};
    validateattributes(dim,{'numeric'},{'scalar','nonempty','integer','>=',1,'<=',2})    
end


meanvec = mean(X,dim);
otherdim = ~(dim-1)+1;
LS_dist = sum((X-meanvec).^2,otherdim).';
[~,I] = sort(LS_dist);
ind = I(1);
if dim == 2
    y = X(:,ind);
elseif dim==1
    y=X(ind,:);
end
if nargout>1
    varargout{1}=ind;
end


%display
%figure, plot(X), hold on, plot(y,'.--g','LineWidth',2), plot(meanvec,'.--r','LineWidth',2)




