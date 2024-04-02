% suite2p cleanup 
% saves F.csv and Fneu.csv in My_V4_Miniscope folder of correct time folder
%16.05.2023
clearvars -except ci_data;
% extract ids of actual cells
%CAREFUL no cell 0 in Matlab!!!!!!!!!!!!!!!!
merged = true;
iscellMerged = iscell;
if merged == true
    mergedCells = [33];
    %removeMerged = [1 13 4 26 7 10 30 53 15 19 38 41 3 5 36 20 49 8 32 17 48 6 24];%actual names in suite2p check if x - 1 is correct in ROI image
    %removeMerged = [0 13]; %g2 14. [3 5 10] 25. [3 8 4]
    removeMerged = [6 2 5 8 105 54 19 120 11 9 37 22 35 84 30 13];
    for i= 1:size(removeMerged,2)
        iscellMerged((removeMerged(i)+1),1) = 0; %mark individual cells of a merge as non cells   
    end    
end
iscellMerged(:,3) = [1:size(iscellMerged,1)];%give cells a label according to suite2p+1 (0=1, 1=2 ect.)
onlyCells =iscellMerged;
nonCells = find(iscellMerged(:,1)==0);Cells =find(iscellMerged(:,1)==1);
onlyCells(nonCells,:) =[]; S2pCells = onlyCells(:,3)-onlyCells(:,1);
%get Frames and neuropil signals from onlyCells
F_Cells = F(Cells,:);
Fneu_Cells = Fneu(Cells,:);
%% sort cells here if they have to be matched with previous session S2pCells is the order they are in
matching = false;
%matchedOrder = [2 0 7 13 6 19 25 26 17 43 45]; % order of cells in S2P!!! so it corresponds to 14.07 [0 1 2 4 6 7 11 15 17 18 31]%g2 25.07 
matchedOrder = [33 5 1 3 6 7 11 2 29 4 8 9 24 12 14 18 16 32 10 17]; % g4 0209 matched to 0109 [0 1 2 3 4 5 6 7 8 9 11 12 13 14 16 17 19 20 22 24]
%Cells is how they are ordered now, matchedOrderList is how they need to be 
%S2p cells is the label !!
if matching == true
    matchedOrderlist = (matchedOrder + 1)' ;
    F_Cells_unmatched = F_Cells; Fneu_Cells_unmatched = Fneu_Cells;
    F_Cells = []; Fneu_Cells = [];
    for i = 1:size(matchedOrderlist,1)
        F_Cells(i,:) = F_Cells_unmatched(find(Cells == matchedOrderlist(i)),:);
        Fneu_Cells(i,:) = Fneu_Cells_unmatched(find(Cells == matchedOrderlist(i)),:);
    end
end
%%
%split signal and neuropil signal into blocks
%iterate through first cell and find first black frame of new block
BlackFrames = find(F_Cells(1,:) < 20);%play with this value until all starts are found
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
%% Import Log file and check for matching frame numbers
%LogFrames = Log(2:2:end,1); FoundFrames = 
%FrameDiff = LogFrames - size(Sessions.F,2)


%%save into original folders as results
topLevelFolder = uigetdir();
% Get a list of all files and folders in this folder.
folders = getFolders(topLevelFolder);% if there are bad recordings, exclude them here!
%folders = folders(2:6,1);
before = 0;
for i = 1:size(Sessions,2)%%%%split here if days were analyzed together
    ImgS = char(folders{(i-before),1});
    savename = [topLevelFolder '\' ImgS '\My_V4_Miniscope\'];
    writematrix(Sessions(i).F, [savename 'F.csv']);
    writematrix(Sessions(i).Fneu, [savename 'Fneu.csv']);
end    








