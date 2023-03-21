% plotting behavior data
clearvars -except entriesBehavior; close all;
%animals 5,2,4
%% transitions
animal = "g12";
%mat = entriesBehavior.bsl.matrix_averg;%# A 5-by-5 matrix of random values from 0 to 1
mat = entriesBehavior.bsl.matrix(:,:,1);
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
b = bar(x,y)
ylabel('Cosine distance from mean');
xticklabels({'g12', 'g10', 'g5', 'g2','g4'});



%% singles vs runs
% makes an average over all trials per cell
fn=fieldnames(events_classif);me = {'_m'};
%loop through the fields to make mean all bsl trials
for ii=1: (numel(fn)-1)
    fn1=fieldnames(events_classif.(fn{ii}));
        for j=6:numel(fn1)
            fn2=fieldnames(events_classif.(fn{ii}).(fn1{j}));
            for k =1:numel(fn2)
                events_classif.(fn{ii}).(fn1{j}).(char(strcat(fn2{k},me))) = nanmean(events_classif.(fn{ii}).(fn1{j}).(fn2{k}),3);
            end    
        end
end    

%% concatenating over averg from all cells (all cells together) separated into singles and runs for all entry types
fn1=fieldnames(events_classif.g5);
    for j=6:numel(fn1)
        fn2 = fieldnames(events_classif.g5.(fn1{j}));
        for jj =4:numel(fn2)
            me= fn2(jj);
            events_classif.all.(strcat((fn1{j}),"_",me)) = cat(2,events_classif.g5.(fn1{j}).(fn2{jj}),events_classif.g2.(fn1{j}).(fn2{jj}),events_classif.g4.(fn1{j}).(fn2{jj}));
        end
    end    




%% plotting traces singles vs first entry of a run

lim = [-3 3];limaverg2= [-1 1];win= 61;
f9=figure;% averages across all sessions
figure(f9);hold on;
subplot(2,4,1),imagesc([0:win],[1:size(events_classif.all.exf_averg_bsl_singles_m,2)],events_classif.all.exf_averg_bsl_singles_m');xline(31,'--b');
        title(['\fontsize{10}exit food singles',' trials: 75-143']);clim(lim);
        xticks([1 11 21 31 41 51 61]); xticklabels({'-3','-2','-1','0', '1', '2', '3'});;ylabel('\fontsize{15}Cells');xlabel('\fontsize{15}time [s]');
        sp7= subplot(2,4,5), stdshade(events_classif.all.exf_averg_bsl_singles_m',0.3,[0 0 1]);set(sp7, 'Position', [0.13,0.43,0.1566,0.15]);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0', '1', '2', '3'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');
        ylim(limaverg2);
subplot(2,4,2),imagesc([0:win],[1:size(events_classif.all.exf_averg_bsl_FR_m,2)],events_classif.all.exf_averg_bsl_FR_m');xline(31,'--b');
        title(['\fontsize{10}exit food 1st run',' trials: 5-10']);clim(lim);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0'});
        sp8= subplot(2,4,6), stdshade(events_classif.all.exf_averg_bsl_FR_m',0.3,[0 0 1]);set(sp8, 'Position', [0.3361,0.43,0.1566,0.15]);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0', '1', '2', '3'});
        ylim(limaverg2);
subplot(2,4,3),imagesc([0:win],[1:size(events_classif.all.exs_averg_bsl_singles_m,2)],events_classif.all.exs_averg_bsl_singles_m');xline(31,'--b');
        title(['\fontsize{10}exit social singles',' trials: 60-141']);clim(lim);
        xticks([1 11 21 31 41 51 61]); xticklabels({' ','1','2',' '});
        sp9= subplot(2,4,7), stdshade(events_classif.all.exs_averg_bsl_singles_m',0.3,[0 0 1]);set(sp9, 'Position', [0.5422,0.43,0.1566,0.15]);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0', '1', '2', '3'}); 
        ylim(limaverg2);
subplot(2,4,4),imagesc([0:win],[1:size(events_classif.all.exs_averg_bsl_FR_m,2)],events_classif.all.exs_averg_bsl_FR_m');xline(31,'--b');
        title(['\fontsize{10}exit social 1st run',' trials: 2-17']);clim(lim);
        xticks([1 11 21 31 41 51 61]);xticklabels({' ','1','2',' '});
        sp10= subplot(2,4,8), stdshade(events_classif.all.exs_averg_bsl_FR_m',0.3,[0 0 1]);set(sp10, 'Position', [0.7484,0.43,0.1566,0.15]);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0', '1', '2', '3'});
        ylim(limaverg2);
        h = axes(f9,'visible','off');
        c = colorbar(h,'Position',[0.93 0.45 0.02 0.45]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 2.5;
        caxis(h,lim);        



%% plotting nr of singles vs nr of runs
% make a mean and std of each type of run in each type of entry
% fn = fieldnames(SR_classif.all);
%     for i =1:numel(fn)
%         for r = 1:size(SR_classif.all.ef,1)
%             SR_classif.all.(fn{i})(:,4) = mean(SR_classif.all.(fn{i}),2);
%             SR_classif.all.(fn{i})(:,5) = std(SR_classif.all.(fn{i}),0,2);
%         end
%     end    

fn = fieldnames(SR_classif.all);
    for i =1:numel(fn)
        SR_classif.all.(fn{i})(:,4) = sum(SR_classif.all.(fn{i}),2);
    end   

%% plotting counts singles vs runs
f1 = figure; lim = [0 325];
figure(f1); 
subplot(1,5,1);hold on; y = SR_classif.all.ef(:,4); b= bar(y); b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];  %errorbar(y,(SR_classif.all.ef(:,5)/sqrt(3)), '.k', LineWidth=1); 
title(['\fontsize{12}food']);ylabel('\fontsize{14}Counts');xlabel('\fontsize{14}Consecutive entries in a run');
xticks([1 2 3 4 5]);xticklabels({'S','r2','r3','r4','r5+'});ylim(lim);
subplot(1,5,2);hold on; y = SR_classif.all.er(:,4); b= bar(y); b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5]; %errorbar(y,(SR_classif.all.er(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{12}run']);xticks([1 2 3 4 5]);;ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,3);hold on; y = SR_classif.all.ee(:,4); b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.ee(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{12}expl']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,4);hold on; y = SR_classif.all.es(:,4); b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.es(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{12}soc']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});
subplot(1,5,5);hold on; y = SR_classif.all.ed(:,4); b= bar(y);  b.FaceColor = 'flat'; b.CData(1,:) = [.5 0 .5];%errorbar(y,(SR_classif.all.ed(:,5)/sqrt(3)),  '.k', LineWidth=1); 
title(['\fontsize{12}drink']);xticks([1 2 3 4 5]);ylim(lim);xticklabels({'S','r2','r3','r4','r5+'});



%% all trials from one cell
fnEA = fieldnames(events_classif);
anim = 2;
fnEAanimal = fieldnames(events_classif.(fnEA{anim}));% baselined individual traces
Event = 4; CellID =9;%5%6
%
Cell1 = events_classif.(fnEA{anim}).(fnEAanimal{Event}).singles;Cell1 = Cell1(:,CellID,:);
C1 = reshape(Cell1 ,[],size(Cell1,3));
%entry
figure; plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');

Cell1 = events_classif.(fnEA{anim}).(fnEAanimal{Event}).FR;Cell1 = Cell1(:,CellID,:);
C1 = reshape(Cell1 ,[],size(Cell1,3)); 
%entry
figure; plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');



%% loop to plot all cells from g2 for es
fnEA = fieldnames(events_classif);
anim = 2;Event = 4;
fnEAanimal = fieldnames(events_classif.(fnEA{anim}));col = 6; st = 7;
for CellID = 8:13

    Cell1 = events_classif.(fnEA{anim}).(fnEAanimal{Event}).singles;Cell1 = Cell1(:,CellID,:);
    C1 = reshape(Cell1 ,[],size(Cell1,3));

    subplot(2,col,CellID-st);
    plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
    xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);ylim([-5 10]);
    ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');

    Cell1 = events_classif.(fnEA{anim}).(fnEAanimal{Event}).FR;Cell1 = Cell1(:,CellID,:);
    C1 = reshape(Cell1 ,[],size(Cell1,3));

    subplot(2,col,(CellID-st)+col);
    plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
    xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);ylim([-5 10]);
    ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)'); 
end


%% loop to plot all cells from g2 for exit social!!!!!!!!!!!!
fnEA = fieldnames(events_classif);
anim = 2;Event = 9;
fnEAanimal = fieldnames(events_classif.(fnEA{anim}));col = 7; st = 0;
for CellID = 1:7

    Cell1 = events_classif.(fnEA{anim}).(fnEAanimal{Event}).singles;Cell1 = Cell1(:,CellID,:);
    C1 = reshape(Cell1 ,[],size(Cell1,3));

    subplot(1,col,CellID-st);
    %plot(C1, LineWidth=1);hold on; 
    plot((mean(C1,2)),'b' ,LineWidth=3);hold on;
    Cell1 = events_classif.(fnEA{anim}).(fnEAanimal{Event}).FR;Cell1 = Cell1(:,CellID,:);
    C1 = reshape(Cell1 ,[],size(Cell1,3));
    plot((mean(C1,2)),'m' ,LineWidth=3);
    xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1', '0', '1', '2', '3'});xlabel('\fontsize{12}time from event[s]');xlim([1 61]);ylim([-2 2.5]);
    ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');
end




%% loop to plot all cells from g4 for ef 23cells
fnEA = fieldnames(events_classif);
anim = 3;Event = 1;
fnEAanimal = fieldnames(events_classif.(fnEA{anim}));col = 7; st = 17;
for CellID = 22:23%8:14

    Cell1 = events_classif.(fnEA{anim}).(fnEAanimal{Event}).singles;Cell1 = Cell1(:,CellID,:);
    C1 = reshape(Cell1 ,[],size(Cell1,3));

    subplot(2,col,CellID-st);
    plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
    xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);
    ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');

    Cell1 = events_classif.(fnEA{anim}).(fnEAanimal{Event}).FR;Cell1 = Cell1(:,CellID,:);
    C1 = reshape(Cell1 ,[],size(Cell1,3));

    subplot(2,col,(CellID-st)+col);
    plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
    xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);
    ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');
end









