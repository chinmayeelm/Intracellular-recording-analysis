function [status, msg] = organiseFiles(rootFolderPath, folderName,fileName)
%ORGANISEFILES This function 
%                   1. creates nested raw and MothID folders
%                   2. Moves files from source folder to M1 or M2 folders
%                   based on moth ID

sourceFolderPath = join([rootFolderPath '\' folderName '\'],'');
cd (sourceFolderPath);
fileNameparts = split(fileName, '_');

mothFolderID = fileNameparts(1);
if ~isfolder("raw") 
    mkdir("raw");
    addpath("raw");   
end 

cd raw

if ~isfolder(mothFolderID)
    mkdir(mothFolderID);
    addpath(mothFolderID);
end

cd (sourceFolderPath);
destinationFolderPath = join([sourceFolderPath 'raw\' mothFolderID], '');
[status, msg] = movefile(join([mothFolderID "_*"],''), destinationFolderPath);

end

