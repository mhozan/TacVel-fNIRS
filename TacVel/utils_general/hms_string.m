function [timestring] = hms_string(SEC)
%[timestring] = HMS_STRING(SEC) returns an string output in the format MM:SS or HH:MM:SS
%
%   Example 1:
%      SEC = 3607;
%      [timestring] = hms_string(SEC);
%       timestring =
%                        01:00:07
%   Example 2:
%      SEC = 1223;
%      [timestring] = hms_string(SEC);
%       timestring =
%                        20:23
%   Example 3:
%      SEC = [1223,-3607];
%      [timestring] = hms_string(SEC);
%       timestring =
%                        00:20:23
%                       -01:00:07
%   Mohsen hozan@mit.edu 4/21/2016

validateattributes(SEC, {'numeric'},{'vector'})
SEC = SEC(:);
sgnstr = num2str(sign(SEC)); %to handle negative input times too
if ~isvector(sgnstr)
    sgnstr = sgnstr(:,1); %take the sign only, throw away the number
else %this happens when SEC is all positive
    sgnstr = [];
end
hh = fix(SEC / 3600);  %For positive SEC, the behavior of fix is the same as floor. For negative SEC, the behavior of fix is the same as ceil.
ss = SEC - 3600*hh;
mm = fix(ss / 60);
ss = round(ss - 60*mm);
if any(hh~=0) %
    timestring = [sgnstr,   num2str(abs([hh,mm,ss])   ,'%02.0f:%02.0f:%02.0f')];
else %exclude the hour string
    timestring = [sgnstr,   num2str(abs([mm,ss])      ,'%02.0f:%02.0f')];
end



