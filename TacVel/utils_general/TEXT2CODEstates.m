function statescode = TEXT2CODEstates(stateschar)
%hozan@mit.edu
%2/22/2016
if ~ischar(stateschar)  %must be char AND a vector, disp error otherwise
    statescode=[];
    disp('unknown input type: input must be "char".')
    return
end
if ~isvector(stateschar)  %must be char AND a vector, disp error otherwise
    statescode=[];
    disp('dimension error: input must be a vector of characters.')
    return
end
stateschar                  = stateschar(:);
stateschar                  = stateschar.';
statescode                  = 8*ones(length(stateschar),1); %8: No Stage: According to Auto Sleep Scoring Script for SPIKE2
WAKEindices                 = regexpi(stateschar,'W');
statescode(WAKEindices)     = 1; %1:W
NREMindices                 = regexpi(stateschar,'N');
statescode(NREMindices)     = 2; %2:N
REMindices                  = regexpi(stateschar,'R');
statescode(REMindices)      = 3; %3:R
DOUBTindices                = regexpi(stateschar,'D');
statescode(DOUBTindices)	= 4; %4:D

if ~isempty(DOUBTindices)
    fprintf('WARNING!!!: %3d DOUBT epochs detected. Consider rescoring via Spike2.\n',length(DOUBTindices))
end


