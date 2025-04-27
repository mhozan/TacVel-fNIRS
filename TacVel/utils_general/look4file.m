function outSTR=look4file(strIN,pathname)
%returns inSTR into outSTR if inSTR refers to an existing file.
%returns given pathname+inSTR into outSTR if pathname\inSTR refers to a file.
%returns '' otherwise.
%mohsen Feb 2016
if iscellstr(strIN)
    strIN = [strIN{:}];
end


if isempty(strfind(strIN,';'))
%     assignin('base','pathname',pathname)
%     assignin('base','strIN',strIN)
    if isnan(strIN)
        %skip
        outSTR    	=   '';
    elseif doesitexist(strIN)
        outSTR    	=   strIN;
    elseif doesitexist(fullfile(pathname,strIN))
        outSTR    	=   fullfile(pathname,strIN);
    else
        if ~isempty(strIN)
            fprintf(['could not locate the file: "', strIN, '".\n'])
        end
        outSTR    	=   '';
    end
else %multiple files in the string separated via semicolon delimitter
    nonsplitted_filenames       =   stringSEMICOLdelimitter(strIN);
%     assignin('base','nonsplitted_filenames',nonsplitted_filenames)
    c1=1;
    for i=1:size(nonsplitted_filenames,2)
        filename= nonsplitted_filenames{1,i};
        if doesitexist(filename)
            outSTR(c1).filename    	=   filename;
            c1=c1+1;
        elseif doesitexist(fullfile(pathname,filename))
            outSTR(c1).filename    	=   fullfile(pathname,filename);                       
            c1=c1+1;
        else
            disp('could not locate the file:')
            disp(filename)
            outSTR(c1).filename    	=   [];
        end
%         
%         
%         if ~isempty(look4file(nonsplitted_filenames{1,i},pathname))
%             c1=c1+1;
%             outSTR(c1).filename  =   look4file(nonsplitted_filenames{1,i},pathname);
%         end
    end
end
end