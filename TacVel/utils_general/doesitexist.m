function [logicoutput]=doesitexist(inputstringz)
%checks the existence of input files for NLX3
%hozan@mit.edu Feb 2016
if iscellstr(inputstringz)
    N=length(inputstringz);
    testall=true(N,1);
    for n=1:N
        if exist(inputstringz{n},'file')~=2
            %             error('Could not locate the file: \n %s',primaryinput)
            testall(n)=false;
            if ~isempty(inputstringz{n})
                warning('Could not locate the file #%3d:\n\t',n);
                disp(inputstringz{n})
            else %empty primaryCSC input string; or empty row(s) in primaryCSC column of the master excel file.
                warning('\nFile #%3d is an empty string!\t A primaryCSC file is a mandatory input.',n)
%                 warning(wrnstr)
            end
        end
    end
    if all(testall)
        logicoutput=true;
    else
        logicoutput=false;
    end
elseif isstruct(inputstringz)
    N=length(inputstringz);
    feeldz =fieldnames(inputstringz);
    testall=true(N,1);
    for n=1:N
        for fld=1:length(feeldz)
            %         filename=getfield(inputstringz(n),feeldz{fld});
            filename=inputstringz(n).(feeldz{fld});
            if ~ischar(filename)
                testall(n)=false;
                disp('doesitexist cannot tell you whether a number exist or not; that''s too philosophical, doesitexist''s afraid.')
            elseif exist(filename,'file')~=2
                %             error('Could not locate the file: \n %s',primaryinput)
                testall(n)=false;
                fprintf('Could not locate the file #%3d: ',n)
                disp(filename)
                %             if ~isempty(inputstringz{n})
                %                 warning(['Could not locate the file #%3d:\n\t',inputstringz{n}],n)
                %             else %empty primaryCSC input string; or empty row(s) in primaryCSC column of the master excel file.
                %                 warning('\nFile #%3d is an empty string!\t A primaryCSC file is a mandatory input.',n)
                %             end
            end
        end
    end
    if all(testall)
        logicoutput=true;
    else
        logicoutput=false;
    end
elseif ischar(inputstringz)
    logicoutput=true;
    if exist(inputstringz,'file')~=2
        %             error('Could not locate the file: \n %s',primaryinput)
        logicoutput=false;
        %       	warning(['Could not locate the file: ',inputstringz])
        %      	disp(inputstringz);
    end
else %other types of inputs, e.g. double
    logicoutput = false;
    disp('doesitexist: the only supported types of input are cell,struct and char. have a nice day.')
end


end
