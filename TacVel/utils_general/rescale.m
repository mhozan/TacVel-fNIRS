function [scaledoutput]=rescale(input,varargin)
% [scaledoutput]=rescale(input) rescales the range of the input  to [0,1] range
% [scaledoutput]=rescale(input,[a b]) rescales the range of the input  to
% [a,b]; b must be greater than a;
% mohsen hozan@mit.edu 10/14/2016
rng=[0 1]; %rescaled range
validateattributes(input,{'numeric'},{'nonempty'})

if nargin==2
    validateattributes(varargin{1},{'numeric'},{'vector','numel',2,'increasing'})
    rng = varargin{1};    
end

% scaledoutput = (input-min(input(:))) ./ (max(input(:)-min(input(:))));
scaledoutput = (input-min(input(:))) ./ range(input(:));
scaledoutput = rng(1) + range(rng)*scaledoutput;