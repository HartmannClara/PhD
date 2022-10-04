%plotting pellet experiment results and comparing activity before and after
clear all;close all;
% f1 = figure
% f2 = figure
% f3 = figure
% f4 = figure
% f5 = figure
%f6 = figure
set(0,'defaultAxesFontSize',8)
%% import and clean up of results files
%bgr is subtracted during motion correction preprocessing data is raw df/f
% 0601_FD, 1101_FD, 1401_FD, 1801_FD, 1701_nFD, 1801_nFD, 1901_nFD, 2001_nFD
%% g2
datasets= table(['2022_04_06'; '2022_04_07'; '2022_04_08']);
dates = table(['11_01_21'; '10_50_46'; '10_47_57']);
foodframes = [5000, 5000, 5000];
f1 = figure;

for s=1:3
    stack=cellstr(table2cell(datasets(s,1)));
    stack2=cellstr(table2cell(dates(s,1)));
    filename=['D:\Exp_1_PFC-Pellet\g2\' stack{1, 1} '\' stack2{1, 1} '\My_V4_Miniscope\Results.csv'];
    activity_raw = readtable(filename);
    activity_raw = activity_raw{:,3:4:end};
    activity_zs = zscore(activity_raw,1,1);

    f2 = figure
    figure(f2);
    stackedplot(activity_zs);
    xlabel('Frames (10fps)');
    title('Activity (z-scored)');


    food_frame = foodframes(s);
    activity_bfF=activity_zs(1:food_frame,:);
    activity_aF=activity_zs(food_frame:end,:);

    mean_bfF= mean(activity_bfF);
    mean_aF= mean(activity_aF);
    
    comp(:,1) = mean_bfF' ;
    comp(:,2) = mean_aF' ;
    
    clear ylabel
    figure(f1);hold on;
    subplot(2,3,s),coordLineStyle = 'k.';boxplot(comp(1:23,1:2), 'Symbol', coordLineStyle, 'labels', {'before food', 'after food'}); hold on;
    parallelcoords(comp(1:23,1:2), 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
      'Marker', '.', 'MarkerSize', 10);
    ylabel('Mean activity (z-scored)');

    [h1,p1]= ttest(comp(:,1),comp(:,2))% over all cells: overall activity bfF diff from aF?



end

%% g1
clearvars -except f1 
datasets= table(['2022_01_17'; '2022_01_18'; '2022_01_20']);
dates = table(['10_12_55'; '17_15_04'; '15_11_38']);
foodframes = [5443, 5191, 5275];
for s=1:3
    stack=cellstr(table2cell(datasets(s,1)));
    stack2=cellstr(table2cell(dates(s,1)));
    filename=['D:\Exp_1_PFC-Pellet\g1\01_Raw\no_FD\' stack{1, 1} '\' stack2{1, 1} '\My_V4_Miniscope\Results.csv'];
  
    activity_raw = readtable(filename);
    activity_raw = activity_raw{:,3:4:end};
    activity_zs = zscore(activity_raw,1,1);

    f2 = figure
    figure(f2);
    stackedplot(activity_zs);
    xlabel('Frames (10fps)');
    title('Activity (z-scored) g1');


    food_frame = foodframes(s);
    activity_bfF=activity_zs(1:food_frame,:);
    activity_aF=activity_zs(food_frame:end,:);

    mean_bfF= mean(activity_bfF);
    mean_aF= mean(activity_aF);
    
    comp(:,1) = mean_bfF' ;
    comp(:,2) = mean_aF' ;
    
    
    figure(f1);hold on;
    subplot(2,3,s+3),coordLineStyle = 'k.';boxplot(comp(1:23,1:2), 'Symbol', coordLineStyle, 'labels', {'before food', 'after food'}); hold on;
    parallelcoords(comp(1:23,1:2), 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
      'Marker', '.', 'MarkerSize', 10);
    ylabel('Mean activity (z-scored)');


    [h2,p2]= ttest(comp(:,1),comp(:,2))% over all cells: overall activity bfF diff from aF?
end



