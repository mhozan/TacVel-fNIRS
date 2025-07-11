function data = loadDirectory2( rootFolder, folderHierarchy, loadFunc, fileExt )
%modified to automatically load adicht files as auxillary data via loadNIRX2
if nargin < 4
%     fileExt = {'.nirs','.oxy3','.wl1','Probe*.csv','_fnirs.csv','nir5','TXT'};
    fileExt = {'.wl1'}; 
end
if nargin < 3 || isempty(loadFunc)
    loadFunc = {@(file)nirs.cnl.loadNIRx2(file,false)};
%     loadFunc = {@nirs.io.loadDotNirs,@nirs.io.loadOxy3,@(file)nirs.cnl.loadNIRx2(file,false),@nirs.io.loadHitachi,@nirs.io.loadHitachiV2,@nirs.io.loadNIR5,@nirs.io.loadShimadzu};
end

if(~iscell(fileExt)); fileExt={fileExt}; end
if(~iscell(loadFunc)); loadFunc={loadFunc}; end

% remove trailing file separator
if rootFolder(end) == filesep
    rootFolder = rootFolder(1:end-1);
end

% default folder structure
if nargin < 2
    folderHierarchy = {};
end

if(~iscellstr(folderHierarchy))
    folderHierarchy={folderHierarchy};
end


% all files in subdirectory with correct extension
data = nirs.core.Data.empty;
for i=1:length(fileExt)
    if(contains(rootFolder,'*'))
        files = rdir(fullfile(rootFolder,'*',['*' fileExt{i}]));
    else
        files = rdir(fullfile(rootFolder,'**',['*' fileExt{i}]));
    end
    
    for iFile = 1:length( files )
        
        % load using load function
%         try
            %            disp(['loading: ' files(iFile).name]);
            tmp = loadFunc{i}( files(iFile).name );
%         catch err
%             if(~strcmp(fileExt{i},'TXT'))
%                 warning(['error reading file: ' files(iFile).name]);
%                 disp(err)
%             end
%             continue;
%         end
        % NIRx data uses folders instead of files... back up one
        if(contains(fileExt{i},'.wl'))
            files(iFile).name=[fileparts(files(iFile).name) filesep];
        end
        if ~isempty(tmp)
            
            if(contains(func2str(loadFunc{i}),'nirs.cnl.loadNIRx2') && ...
                    strcmp(func2str(loadFunc{i}),'@(file)nirs.cnl.loadNIRx2(file,false)') && isempty(data))
                disp('Loading NIRx file geometry from:')
                disp(['     ' files(iFile).name]);
                disp('      Note: This registration will be used for all subjects');
                disp('      To load all use "loadDirectory(<>,<>,@(file)nirs.cnl.loadNIRx2(file,true))"');
                tmp = nirs.cnl.loadNIRx2(files(iFile).name,true);
                probe=tmp.probe;
                loadFunc{i} = @(file)nirs.cnl.loadNIRx2(file,false);
                
            end
            
            if exist('probe','var')
                tmp.probe=probe;
            end
            
            data(end+1) = tmp;
            
            % moved to Nirs.cnl.loadNIRX2
            % % % split filename on separators
            % % fsplit = strsplit( files(iFile).name, filesep );
            % % rsplit = strsplit( rootFolder, filesep );
            % % 
            % % % add demographics variables based on folder names
            % % demo = fsplit(length(rsplit)+1:end-1);
            % % data(end).description=files(iFile).name;
            % % for iDemo = 1:min(length(folderHierarchy),length(demo))
            % %     data(end).demographics(folderHierarchy{iDemo}) = demo{iDemo};
            % % end
            % % data(end).demographics({'isFirstSession'}) = contains(data(end).demographics({'Session'}),'1');
            % % 
        end
    end
    
    
end
data = data';
end

