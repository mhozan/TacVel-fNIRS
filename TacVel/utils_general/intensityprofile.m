function [intensityprofile_curve,varargout] = intensityprofile(spctrgrm,varargin)
% [intensityprofile_curve] = INTENSITYPROFILE(spctrgrm); returns the
% intensity curve of the spctrgram along the frequency axis (dim 2)
% [intensityprofile_curve] = INTENSITYPROFILE(spctrgrm,1); returns the
% intensity curve of the spctrgram along the time axis (dim 1)
% [~,intensityprofile_matrix] = INTENSITYPROFILE(spctrgrm); returns the
% intensity matrix which is the same size as spctrgrm and is consisted of copies of the intensityprofile_curve in each column(or row if dim ==1). 
%mohsen hozan Oct 2016
narginchk(1,2)
nargoutchk(1,2)
dim = 2; %default dimension : 2. the frequency axis of the spctrgrm (rows)
validateattributes(spctrgrm,    {'numeric'},    {'2d','nonempty'});
if nargin==2
    validateattributes(varargin{1},    {'numeric'},    {'scalar','nonempty','integer','>=',1,'<=',2});
    dim = varargin{1};
end
%%

intensityprofile_curve = sum(spctrgrm,dim);
intensityprofile_curve = intensityprofile_curve./max(intensityprofile_curve);

if dim == 1
    intensityprofile_matrix = repmat(intensityprofile_curve,size(spctrgrm,dim),1);
else %dim == 2
    intensityprofile_matrix = repmat(intensityprofile_curve,1,size(spctrgrm,dim));
end
if nargout == 2
    varargout{1} = intensityprofile_matrix;
end

