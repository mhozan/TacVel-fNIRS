function [intensityprofile_curve,varargout] = intensityprofile2(f_vec,varargin)
% similar to intensityprofile but creates the profile based on an analytic
% formula.
% [intensityprofile_curve] = INTENSITYPROFILE2(f_vec); returns the
% intensity curve of the spctrgram along the frequency axis (dim 2)
% [intensityprofile_curve] = INTENSITYPROFILE(spctrgrm,1); returns the
% intensity curve of the spctrgram along the time axis (dim 1)
% [~,intensityprofile_matrix] = INTENSITYPROFILE(spctrgrm); returns the
% intensity matrix which is the same size as spctrgrm and is consisted of copies of the intensityprofile_curve in each column(or row if dim ==1). 
% See also 
% INTENSITYPROFILE
%mohsen hozan Oct 2016
narginchk(1,2)
nargoutchk(1,2)
validateattributes(f_vec,    {'numeric'},    {'vector','nonempty','nonnegative'});

%%
alfa=-4.0;
% zeroshift=+0.5;
zeroshift=-.0;
f_vec=f_vec+zeroshift;
intensityprofile_curve = f_vec.^8.*exp(alfa*f_vec);
intensityprofile_curve = rescale2(intensityprofile_curve,[0.2 1]);


% if spectrogram == 1
%     intensityprofile_matrix = repmat(intensityprofile_curve,size(f_vec,spectrogram),1);
% else %dim == 2
%     intensityprofile_matrix = repmat(intensityprofile_curve,1,size(f_vec,spectrogram));
% end
% if nargout == 2
%     varargout{1} = intensityprofile_matrix;
% end

