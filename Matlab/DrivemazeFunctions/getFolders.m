function folders = getFolders(topLevelFolder)
%UNTITLED12 getFolders
%   gives table of imaging folder names 
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = ({subFolders(3:end).name})';
folders = cell2table(subFolderNames);
end