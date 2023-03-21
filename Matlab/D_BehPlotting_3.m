% plotting behavior data
clearvars -except entriesBehavior; close all;
%animals 5,2,4
%% transitions
animal = "g5_2";
%animal = "all";
%mat = entriesBehavior.Swap.matrix_averg;%# A 5-by-5 matrix of random values from 0 to 1
mat = entriesBehavior.bsl.matrix(:,:,6);
clims = [0 0.72];
imagesc(mat,clims);            %# Create a colored plot of the matrix values
%colormap(flipud(gray));  %# Change the colormap to gray (so higher values are
                       %#   black and lower values are white)
textStrings = num2str(mat(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:6);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');

midValue = mean(get(gca,'CLim'));  %# Get the middle value of the color range
textColors = repmat(mat(:) < midValue,1,3);  %# Choose white or black for the
                                             %#   text color of the strings so
                                             %#   they can be easily seen over
                                             %#   the background color
set(hStrings,{'Color'},num2cell(textColors,2));  %# Change the text colors
set(gca, 'XAxisLocation','top')
xticklabels({'food','run','expl','soc','drink', 'home'});xlabel('Trial + 1');
yticklabels({'food','run','expl','soc','drink', 'home'});ylabel('Trial');
title(animal);

%% cosine distance

all = reshape(entriesBehavior.bsl.matrix_averg, [1, 36]);
g12 = reshape(entriesBehavior.bsl.matrix(:,:,1), [1, 36]);
g10 = reshape(entriesBehavior.bsl.matrix(:,:,2), [1, 36]);
g5 = reshape(entriesBehavior.bsl.matrix(:,:,3), [1, 36]);
g2 = reshape(entriesBehavior.bsl.matrix(:,:,4), [1, 36]);
g4 = reshape(entriesBehavior.bsl.matrix(:,:,5), [1, 36]);


csg12= pdist2(all,g12,'cosine');
csg10= pdist2(all,g10,'cosine');
csg5= pdist2(all,g5,'cosine');
csg2= pdist2(all,g2,'cosine');
csg4= pdist2(all,g4,'cosine');

y=[csg12,csg10, csg5, csg2, csg4];
x=[1:5]
b = bar(x,y); ylim([0 0.35])
ylabel('Cosine distance from mean');
xticklabels({'g12', 'g10', 'g5', 'g2','g4'});

%%
g5_1 = reshape(entriesBehavior.bsl.matrix(:,:,3), [1, 36]);
g5_2 = reshape(entriesBehavior.bsl.matrix(:,:,6), [1, 36]);
g5swap = reshape(entriesBehavior.Swap.matrix(:,:,1), [1, 36]);
g5FEDrem = reshape(entriesBehavior.FEDrem.matrix(:,:,1), [1, 36]);

csg5_2 =pdist2(g5_1,g5_2,'cosine');
csg5Swap =pdist2(g5_1,g5swap,'cosine');
csg5FEDrem =pdist2(g5_1,g5FEDrem,'cosine');

y=[csg5_2,csg5Swap,csg5FEDrem];
x=[1:3]
b = bar(x,y)
ylabel('Cosine distance from baseline (g5_1)');
xticklabels({'bsl_2','g5 Swap', 'g5 FEDrem'});

%%
g12Swap = reshape(entriesBehavior.Swap.matrix(:,:,3), [1, 36]);

csg12Swap =pdist2(g12,g12Swap,'cosine');

y=[csg12Swap];
x=1
b = bar(x,y)
ylabel('Cosine distance from baseline');
xticklabels({'g12 Swap'});


%% plotting nr of singles vs nr of runs
% make a mean and std of each type of run in each type of entry
% fn = fieldnames(SR_classif.all);
%     for i =1:numel(fn)
%         for r = 1:size(SR_classif.all.ef,1)
%             SR_classif.all.(fn{i})(:,4) = mean(SR_classif.all.(fn{i}),2);
%             SR_classif.all.(fn{i})(:,5) = std(SR_classif.all.(fn{i}),0,2);
%         end
%     end    

fn = fieldnames(entriesBehavior.bsl.SR_classif.all);
    for i =1:numel(fn)
        entriesBehavior.bsl.SR_classif.all.(fn{i})(:,6) = sum(entriesBehavior.bsl.SR_classif.all.(fn{i}),2);
    end   
%% plotting counts singles vs runs g5
f1 = figure; lim = [0 150];
figure(f1); 
subplot(1,5,1);hold on; y = entriesBehavior.FEDrem.SR_classif.g5.ef; b= bar(y); b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];  %errorbar(y,(SR_classif.all.ef(:,5)/sqrt(3)), '.k', LineWidth=1); 
title(['\fontsize{14}food']);ylabel('\fontsize{18}Counts');xlabel('\fontsize{18}Consecutive entries in a run');
xticks([1 2 3 4 5]);xticklabels({'S','r2','r3','r4','r5+'});ylim(lim);
subplot(1,5,2);hold on; y = entriesBehavior.FEDrem.SR_classif.g5.er; b= bar(y); b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5]; %errorbar(y,(SR_classif.all.er(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{14}run']);xticks([1 2 3 4 5]);;ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,3);hold on; y = entriesBehavior.FEDrem.SR_classif.g5.ee; b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.ee(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{14}expl']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,4);hold on; y = entriesBehavior.FEDrem.SR_classif.g5.es; b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.es(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{14}soc']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,5);hold on; y = entriesBehavior.FEDrem.SR_classif.g5.ed; b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.ed(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{14}drink']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});


%% plotting counts singles vs runs
f1 = figure; lim = [0 550];
figure(f1); 
subplot(1,5,1);hold on; y = entriesBehavior.bsl.SR_classif.all.ef(:,6); b= bar(y); b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];  %errorbar(y,(SR_classif.all.ef(:,5)/sqrt(3)), '.k', LineWidth=1); 
title(['\fontsize{14}food']);ylabel('\fontsize{18}Counts');xlabel('\fontsize{18}Consecutive entries in a run');
xticks([1 2 3 4 5]);xticklabels({'S','r2','r3','r4','r5+'});ylim(lim);
subplot(1,5,2);hold on; y = entriesBehavior.bsl.SR_classif.all.er(:,6); b= bar(y); b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5]; %errorbar(y,(SR_classif.all.er(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{14}run']);xticks([1 2 3 4 5]);;ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,3);hold on; y = entriesBehavior.bsl.SR_classif.all.ee(:,6); b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.ee(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{14}expl']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,4);hold on; y = entriesBehavior.bsl.SR_classif.all.es(:,6); b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.es(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{14}soc']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,5);hold on; y = entriesBehavior.bsl.SR_classif.all.ed(:,6); b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.ed(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{14}drink']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});


%% consuption counts comparison

%% plotting counts cons vs nc
fn = fieldnames(entriesBehavior.bsl.CNC_classif.all);
    for i =1:numel(fn)
        entriesBehavior.bsl.CNC_classif.all.(fn{i})(:,6) = sum(entriesBehavior.bsl.CNC_classif.all.(fn{i}),2);
    end   
%%

f1 = figure; lim = [0,100];
figure(f1); 
subplot(1,2,1);hold on; y = entriesBehavior.bsl.CNC_classif.g12.ef; b= bar(y);b.FaceColor = 'flat'; ylim(lim);  
b.CData(1,:) = [.5 0 .5];b.CData(2,:) = [.5 0 .5] ; 
title(['\fontsize{14}food']);ylabel('\fontsize{18}Counts'); xticks([1:8]); xticklabels({'s c', 's nc', '1r c','1r nc', 'r c', 'r nc', 'lr c', 'lr nc'});
subplot(1,2,2);hold on; y = entriesBehavior.bsl.CNC_classif.g12.ed; b= bar(y);ylim(lim);  
b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];b.CData(2,:) = [.5 0 .5] ;
title(['\fontsize{14}drink']);ylabel('\fontsize{18}Counts'); xticks([1:8]); xticklabels({'s c', 's nc', '1r c','1r nc', 'r c', 'r nc', 'lr c', 'lr nc'});
%% 
f1 = figure; lim = [0,1.2];
figure(f1); 
subplot(1,2,1);hold on; 
y = [entriesBehavior.Swap.CNC_classif.g12.ef(1,:), entriesBehavior.Swap.CNC_classif.g12.ef(2,:), (entriesBehavior.Swap.CNC_classif.g12.ef(3,:) + entriesBehavior.Swap.CNC_classif.g12.ef(5,:) + entriesBehavior.Swap.CNC_classif.g12.ef(7,:)), (entriesBehavior.Swap.CNC_classif.g12.ef(4,:) + entriesBehavior.Swap.CNC_classif.g12.ef(6,:) + entriesBehavior.Swap.CNC_classif.g12.ef(8,:))]; 
y= [y(1)/(y(1)+y(2)) y(2)/(y(1)+y(2)) y(3)/(y(3)+y(4)) y(4)/(y(3)+y(4)) ];
b= bar(y);b.FaceColor = 'flat'; ylim(lim);  
b.CData(1,:) = [.5 0 .5];b.CData(2,:) = [.5 0 .5] ; 
title(['\fontsize{14}food']);ylabel('\fontsize{18}Percentage'); xticks([1:4]); xticklabels({'s c', 's nc', 'r c','r nc'});
subplot(1,2,2);hold on; 
y = [entriesBehavior.Swap.CNC_classif.g12.ed(1,:), entriesBehavior.Swap.CNC_classif.g12.ed(2,:), (entriesBehavior.Swap.CNC_classif.g12.ed(3,:) + entriesBehavior.Swap.CNC_classif.g12.ed(5,:) + entriesBehavior.Swap.CNC_classif.g12.ed(7,:)), (entriesBehavior.Swap.CNC_classif.g12.ed(4,:) + entriesBehavior.Swap.CNC_classif.g12.ed(6,:) + entriesBehavior.Swap.CNC_classif.g12.ed(8,:))];  
y= [y(1)/(y(1)+y(2)) y(2)/(y(1)+y(2)) y(3)/(y(3)+y(4)) y(4)/(y(3)+y(4)) ];
b= bar(y);ylim(lim);  
b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];b.CData(2,:) = [.5 0 .5] ;
title(['\fontsize{14}drink']);ylabel('\fontsize{18}Percentage'); xticks([1:4]); xticklabels({'s c', 's nc', 'r c','r nc'});


%% consumption percentages 






