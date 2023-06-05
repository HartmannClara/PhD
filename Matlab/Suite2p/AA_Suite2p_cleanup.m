% suite2p cleanup 
% saves F.csv and Fneu.csv in My_V4_Miniscope folder of correct time folder
%16.05.2023
clearvars -except ci_data;
% extract ids of actual cells
%CAREFUL no cell 0 in Matlab!!!!!!!!!!!!!!!!
mergedCells = [57:66];
removeMerged = [1 13 4 26 7 10 30 53 15 19 38 41 3 5 36 20 49 8 32 17 48 6 24];%actual names in suite2p check if x - 1 is correct in ROI image
iscellMerged = iscell;
for i= 1:size(removeMerged,2)
    iscellMerged((removeMerged(i)+1),1) = 0; %mark individual cells of a merge as non cells   
end    

iscellMerged(:,3) = [1:size(iscellMerged,1)];%give cells a label according to suite2p+1 (0=1, 1=2 ect.)
onlyCells =iscellMerged;
nonCells = find(iscellMerged(:,1)==0);Cells =find(iscellMerged(:,1)==1);
onlyCells(nonCells,:) =[]; S2pCells = onlyCells(:,3)-onlyCells(:,1);
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
%plot(F_Cells(1,:))
BlackFrames(length(BlackFrames)+1) = (size(F_Cells,2)+1);% make the end of the last block
%cut signal (F_Cells and Fneu_Cells into blocks according to black frames
for i = 1:(size(BlackFrames,2)-1)
    Sessions(i).F = F_Cells(:,[BlackFrames(i):BlackFrames(i+1)-1]);
    Sessions(i).Fneu = Fneu_Cells(:,[BlackFrames(i):BlackFrames(i+1)-1]);
end    


%%save into original folders as results
topLevelFolder = uigetdir();
% Get a list of all files and folders in this folder.
folders = getFolders(topLevelFolder);% if there are bad recordings, exclude them here!
folders = folders(2:6,1);
before = 0;
for i = 1:5%size(Sessions,2)%%%%split here if days were analyzed together
    ImgS = char(folders{(i-before),1});
    savename = [topLevelFolder '\' ImgS '\My_V4_Miniscope\'];
    writematrix(Sessions(i).F, [savename 'F.csv']);
    writematrix(Sessions(i).Fneu, [savename 'Fneu.csv']);
end    








