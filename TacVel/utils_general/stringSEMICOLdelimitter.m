function STRZcell= stringSEMICOLdelimitter(strIN)
%uses the delimitter in strings to look for filenames. Output is always a
%cell array of strings. Input is always a string.
% mohsen 2/17/2016
%sometimes the input is a cell array containing strings
if iscellstr(strIN)
    strIN = [strIN{:}];
end
if isnan(strIN)
    STRZcell = [];
    return
end

% validateattributes(strIN,{'char'},{'nonempty'})

validateattributes(strIN,{'char'},{})

% getting rid of "char(10) aka'\n' aka 'return key' aka 'ctrl/alt Enter' " in string names:
strIN = strrep(strIN,char(10),'');
% strIN=strrep(strIN,char(32),''); %removing spaces %better not to use. uncommenting this line makes NLX3 unable to process independent filepath-names outside the directory of primaryCSC.


if ~isempty(strfind(strIN,';')) %semicolon delimitter detected, split the string
    STRZcell = strsplit(strIN,';');
    STRZcell = STRZcell(~cellfun('isempty',STRZcell)); %Remove empty cells (This happens sometimes.e.g. 'TT1.ntt; ;TT2.ntt;')
else %no semicolon delimitter is used in the string--> return the input string 
    STRZcell = {strIN};
end


