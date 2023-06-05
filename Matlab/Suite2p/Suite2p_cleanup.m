% suite2p cleanup 
% saves F.csv and Fneu.csv in My_V4_Miniscope folder of correct time folder
%16.05.2023
% extract ids of actual cells
%CAREFUL no cell 0 in Matlab!!!!!!!!!!!!!!!!
mergedCells = [79:81];
removeMerged = [5,4,61,12,39,29,65,69,26];
iscellMerged = iscell;
for i= 1:size(removeMerged,2)
    iscellMerged((removeMerged(i)+1),1) = 0; %mark individual cells of a merge as non cells   
end    
iscellMerged(:,3) = [1:size(iscellMerged,1)];%give cells a label according to suite2p+1 (0=1, 1=2 ect.)
onlyCells =iscellMerged;
nonCells = find(iscellMerged(:,1)==0);Cells =find(iscellMerged(:,1)==1);
onlyCells(nonCells,:) =[];
%get Frames and neuropil signals from onlyCells
F_Cells = F(Cells,:);
Fneu_Cells = Fneu(Cells,:);
%split signal and neuropil signal into blocks
%iterate through first cell and find first black frame of new block
BlackFrames = find(F_Cells(1,:) < 20);
for i = 1:(size(BlackFrames,2))
    if BlackFrames(i+1) == BlackFrames(i)+1 
       BlackFrames(i+1) = [];
    end 
    if BlackFrames(i+1) == BlackFrames(i)+2 
       BlackFrames(i+1) = [];
    end 
end 
BlackFrames(length(BlackFrames)+1) = (size(F_Cells,2)+1);% make the end of the last block
%cut signal (F_Cells and Fneu_Cells into blocks according to black frames
for i = 1:(size(BlackFrames,2)-1)
    Sessions(i).F = F_Cells(:,[BlackFrames(i):BlackFrames(i+1)-1]);
    Sessions(i).Fneu = Fneu_Cells(:,[BlackFrames(i):BlackFrames(i+1)-1]);
end    


%%save into original folders as results
topLevelFolder = uigetdir();
% Get a list of all files and folders in this folder.
folders = getFolders(topLevelFolder);
for i = 1:size(Sessions,2)
    ImgS = char(folders{i,1});
    savename = [topLevelFolder '\' ImgS '\My_V4_Miniscope\'];
    writematrix(Sessions(i).F, [savename 'F.csv']);
    writematrix(Sessions(i).Fneu, [savename 'Fneu.csv']);
end    








