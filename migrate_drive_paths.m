%%% Google drive changed the way it syncs to the local PC, so all my old
%%% file paths are defunct. This function searches all files in the current
%%% folder for the old file path string, and replaces it with the new one

%% Setup
path_working = pwd; % current directory
filesAndFolders = dir(path_working); % find all files and folders in the directory
filesInDir = filesAndFolders(~([filesAndFolders.isdir]));  % Returns only the files in the directory
stringToBeFound = 'C:\Users\esmyt\Google Drive\'; % search for this within each file
replacementString = 'G:\My Drive\'; % replace with this
numOfFiles = length(filesInDir); % number of candidate files

%% Loop through files
file_counter = 0;
str_counter = 0;
for ii = 1:numOfFiles
    filename = filesInDir(ii).name;
    if filename(end-1:end) == '.m'
        % This is a matlab script
        file_counter = file_counter+1;
        
        %%% See if script has the string
        fid = fopen(filename);
        while(~feof(fid))                                           % Execute till EOF has been reached
            contentOfFile = fgetl(fid);                             % Read the file line-by-line and store the content
            found = strfind(contentOfFile,stringToBeFound);         % Search for the stringToBeFound in contentOfFile
            if ~isempty(found)
                str_counter = str_counter+1;
            end
        end
        fclose(fid); 
        
        %%% Read script
        fid = fopen(filename,'r');
        f=fread(fid,'*char')';
        fclose(fid);
        
        %%% Make edits
        f = strrep(f,stringToBeFound,replacementString);
        
        %%% Delete old file
        delete(filename)
        
        %%% Write new script
        fid  = fopen(filename,'w');
        fprintf(fid,'%s',f);
        fclose(fid);
        
    end
    
end

%% Display
disp(['Found and replaced ' num2str(str_counter) ' instances of the old GDrive file path in ' ...
    num2str(file_counter) ' files'])

