function [LS_dist,varargout] = curvedistance(X,varargin)

%CURVEDISTANCE finds the least-square distance between an array of vectors, 
% and a reference vector of the same length. 
%
% Syntax
% LS_dist = curvedistance(X)
% LS_dist = curvedistance(X,refvec)
% [LS_dist,indsorted] = curvedistance(X,refvec);


%   Input Variables:
%   X: An [MxN] matrix, where N is the number of curves (vectors) and M is the
%   length of each curve.
%   refvec : An [Mx1] reference vector to calculate the distance from. If empty, refvec will be
%   the mean of all vectors in X.
%
%   Output Variables:
%   LS_dist: a [1xN] vector consisting of all the least square distances.
%   indsorted: a [1xN] vector consisting of the sorted indices of vectors, from
%   the nearest to the refvec to the most distant.
%
%
% see also
% Communication Neuroscience Laboratories,
% University of Nebraska-Lincoln,
% Mohsen Hozan mhozan2@unl.edu
% Date 2019
narginchk(1,2)
nargoutchk(1,2)
[M,N] = size(X);
refvec = mean(X,2);
if nargin >= 2
	refvec = varargin{1};
    validateattributes(refvec,{'numeric'},{'vector','nonempty','size',[M,1]})
end



LS_dist = sum((X-refvec).^2,1).';

if nargout>1
    [~,indsorted] = sort(LS_dist);
    varargout{1}=indsorted;
end






