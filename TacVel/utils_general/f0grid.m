function [indices] = f0grid(f_vec,f0,varargin)
%[indices] = F0GRID(f_vec,f0) takes the frequency vector and a fundamental frequency, returns the
%indices of the F0 and all its harmonics(integer multiples) in the f_vec.
%[indices] = F0GRID(f_vec,f0,width) also controls the number of neighbor
%elements to F0 or its harmonics. default width is 1, i.e. only one index
%matches with each frequency.
% mohsen hozan 10/13/2016
validateattributes(f_vec,{'numeric'},{'nonempty','vector','real'})
validateattributes(f0,{'numeric'},{'nonempty','scalar','real','<=',max(f_vec),'>=',min(f_vec)})

width = 1;
if nargin==3
    width = varargin{1};
    validateattributes(width,{'numeric'},{'nonempty','scalar','real','integer','>=',1,'<=',length(f_vec)})
end
i=1;
while i*f0 <= max(f_vec)
        [~,indicesF(:,i)]=sort(abs(f_vec-i*f0));
        i=i+1;
end

indices = indicesF(1:width,:);
indices = sort(indices(:));



