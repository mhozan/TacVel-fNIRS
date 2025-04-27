function binnedsignal=hannbinning(inputsignal,WinLength,varargin)
% binnedsignal=HANNBINNING(inputsignal,WinLength) bins the inputsignal to %50
% overlapping bins/segments (via a hann window). WinLength is the length of the bins in samples.
% binnedsignal=OVERLAPBINNING(inputsignal,WinLength,prcntoverlap) bins the
% signal with a prcntoverlap given by the user, a value in [0 1] range.
% Default prcntoverlap is 0.5;
% The end columns of the output matrix might contain NaN.
% see also buffer
% Mohsen Hozan 11/2018

narginchk(2,3)
if nargin==3
    prcntoverlap = varargin{1};
    validateattributes(prcntoverlap,    {'numeric'},    {'scalar','>=',0,'<=',1})
    
else %only two inputs
    prcntoverlap = 0.5; %default overlap is 50%
end
wintype = 'hann';

validateattributes(inputsignal,     {'numeric'},    {'vector'})
validateattributes(WinLength,       {'numeric'},    {'scalar','integer','>',1,'<',length(inputsignal)}) %bin length must be at least 2 to make sense, right?

rect_binned = overlapbinning(inputsignal,WinLength,prcntoverlap);
windough = window(wintype,WinLength);
windoughs = repmat(windough,1,size(rect_binned,2));
binnedsignal = rect_binned.*windoughs;





