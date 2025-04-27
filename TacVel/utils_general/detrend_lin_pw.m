function y = detrend_lin_pw (x,varargin)
%linear piecwise detrending using built-in 'detrend'
% y=detrend_lin_pw(x), linearly detrends the signal x at breakpoints that
% are 1000 samples apart.
% % y=detrend_lin_pw(x,bp_dist), linearly detrends the signal x at breakpoints that
% are 'bp_dist' samples apart.
% see also DETREND
%mohsen hozan 10/16/2016
validateattributes(x,{'numeric'},{'nonempty'})
bp_dist = 2300;
if nargin==2
    validateattributes(varargin{1},{'numeric'},{'nonempty','scalar','integer','>',1})
    bp_dist = varargin{1};
end
bp = 1:bp_dist:length(x);
y = detrend(x,'linear',bp);



