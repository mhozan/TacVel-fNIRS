function binnedsignal=overlapbinning(inputsignal,WinLength,varargin)
% binnedsignal=OVERLAPBINNING(inputsignal,WinLength) bins the inputsignal to %50
% overlapping bins/segments (via a rectangular window). WinLength is the length of the bins in samples.
% binnedsignal=OVERLAPBINNING(inputsignal,WinLength,prcntoverlap) bins the
% signal with a prcntoverlap given by the user, a value in [0 1] range.
% Default prcntoverlap is 0.5;
% The end columns of the output matrix might contain NaN.
% see also buffer
% Mohsen Hozan 9/23/2016
narginchk(2,3)
if nargin==3
    prcntoverlap = varargin{1};
    validateattributes(prcntoverlap,    {'numeric'},    {'scalar','>=',0,'<=',1})
    
else %only two inputs
    prcntoverlap = 0.5; %default overlap is 50%
end

validateattributes(inputsignal,     {'numeric'},    {'vector'})
validateattributes(WinLength,       {'numeric'},    {'scalar','integer','>',1,'<',length(inputsignal)}) %bin length must be at least 2 to make sense, right?
inputsignal = [inputsignal(:);nan(WinLength,1)]; %zero padding the signal at the end


overlap=round(prcntoverlap*WinLength);
if overlap == WinLength
    overlap = WinLength-1;
end


fprintf([mfilename,':\n\t\t\tWindow Length\t\t\t=\t%3.0f\n\t\t\tOverlap(#samples)\t\t=\t%3.1d\n'],WinLength,overlap)
stepforward = WinLength-overlap;
L = length(inputsignal);
nrow = WinLength;
ncol = ceil(L./stepforward);
% ncol = 1+ceil(L./stepforward)-overlap;
% ncol = ceil(L./WinLength);
% binnedsignal = zeros(nrow,ncol);
binnedsignal = zeros(nrow,ncol);
col_startindex = 1;
col_endindex = WinLength;
for col=1:ncol
    binnedsignal(:,col) = inputsignal(col_startindex:col_endindex);
    col_startindex = col_startindex+stepforward;
    col_endindex = col_startindex + WinLength - 1;
    if col_endindex>L;
        binnedsignal(:,col+1:end)=[];
        break;
    end
end



%
%

% binnedsignal =  zeros(nrow,ncol);
% for col=1:ncol
% binnedsignal(:,col) = inputsignal(1+(col-1)*stepforward:1+(col-1)*stepforward+WinLength);
% end



% prepadding = 1i*ones(WinLength-overlap,1);
% prepadding = [];
% binnedsignal = buffer([prepadding;inputsignal(:)],WinLength,overlap); %buffer function zero pads the inputsignal at the beginning; I do not want that.

% binnedsignal(:,1:ceil(WinLength./(WinLength-overlap))) = [];






