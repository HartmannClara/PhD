clear all;close all;
f1 = figure
f2 = figure
f3 = figure
f4 = figure
f5 = figure
%f6 = figure
set(0,'defaultAxesFontSize',25)
%% import and clean up of results files
%ROIs 24-26 are background measures %bgr_mean = (bgr(:,1)+ bgr(:,2)+bgr(:,3))/3;
% just one bgr measurement
% 0601_FD, 1101_FD, 1401_FD, 1801_FD, 1701_nFD, 1801_nFD, 1901_nFD, 2001_nFD
datasets=table(['2022_04_06'; '2022_04_07'; '2022_04_08'])
stack=cellstr(table2cell(imaging_datasets(s,1)));
filename=['D:\Exp_1_PFC-Pellet\g2\' stack{1, 1} '\11_01_21\My_V4_Miniscope\Results.csv'];

activity_raw = readtable(filename);
%activity = activity_raw{:,3:4:end};
%bgr= activity(:,24);%
%activity_bgr = activity(:, 1:23)-bgr;%

figure(f1);
stackedplot(activity_bgr);
xlabel('Frames (10fps)')
title('Activity (bg subtracted)')
%figure(f2);
%stackedplot(activity); 

%% loop to make sliding median
%medians = [];
%for ROI = 1:size(activity_bgr,2)
%    cell_n = [];
%    cell_n = activity_bgr(:,ROI); 
%    medians_cell = movmedian(cell_n, 1200);%window size in frames %120s
%    medians = [medians, medians_cell];     
%end


%df_f_trace = [];
%for cell_median = 1: size(medians,2) % Loop df/f
%    v = [];
%    d = [];
%    v = medians(:,cell_median);%select column from medians
%    d = activity_bgr(:,cell_median);%select column from activity_bgr
%    for f = 1: size(d,1)
%        a = [];
%        a =  v(f) ;
%        b = (d(f) - 0.7*a)/0.7*a;
%        df_f_trace(f,cell_median) = b;
%    end
%end

%figure(f2);
%stackedplot(df_f_trace)
%xlabel('Frames (10fps)')
%title('Activity (dF/F)')


%% z scoring
bgr_zsc = zscore(activity_bgr,1,1);% z-scored using the SD of the population (,1,), mean and SD taken from column(,,1); How far away (in SD) is my value from the mean of the population

figure(f3);
sk = stackedplot(bgr_zsc)
xlabel('Frames (10fps)')
title('Activity (z-scored)')
sk.DisplayLabels = {'Cell 1','Cell 2','Cell 3','Cell 4','Cell 5','Cell 6','Cell 7','Cell 8','Cell 9','Cell 10','Cell 11','Cell 12','Cell 13','Cell 14','Cell 15','Cell 16','Cell 17','Cell 18','Cell 19','Cell 20','Cell 21','Cell 22','Cell 23'}


%% list of cells
cells=1:size(bgr_zsc,2);
pairwise_combinations=nchoosek(cells,2);

%% xcorr
for pair=1:size(pairwise_combinations,1)
    a=bgr_zsc(:,pairwise_combinations(pair,1));
    b=bgr_zsc(:,pairwise_combinations(pair,2));
    [c(:,pair),lags]=xcov(a,b,100,'coeff');%cross-covariance
    pearson(pair)=c(101,pair);
end
%figure(f4);
%plot(pearson);

corr = pairwise_combinations;
corr(:,3)= pearson;

%% loop through correlation results to make matrix for plotting

pearson_tbl = [];%create an empty matrix
for i = 1 : size(corr,1);
    pearson_tbl(corr(i,1),corr(i,2)) = corr(i,3);%coordinates of column 1 and 2 
end
pearson_tbl(pearson_tbl == 0)=NaN;
%M = M';

%figure(f2);
%imagesc(pearson_tbl');
%colorbar
figure(f4);
heatmap(pearson_tbl');
colormap cool;

%% find all highly correlating pairs 

%result=pairwise_combinations(find(pearson>0.5),:)
%result(:,3)=pearson(find(pearson>0.5))

%figure;
%plot(lags,c(:,(find(pearson>0.5))));

%result_tbl = [];
%for i = 1 : size(result,1)
 %   result_tbl(result(i,1),result(i,2)) = result(i,3); 
%end

%result_tbl(result_tbl == 0)=NaN

%figure(f3)
%heatmap(result_tbl')
%colormap cool

%% change analysis window, separate into before and after food pellet
% 0601_FD, 1101_FD, 1401_FD, 1801_FD, 1701_nFD, 1801_nFD, 1901_nFD, 2001_nFD
% 3113, 3123, 3096, 3259, 3245, 3094, 3257, 3143

food_frame = 3096 %%%% change
activity_bfF=bgr_zsc(1:food_frame,:);
activity_aF=bgr_zsc(food_frame:end,:);
clear a b c1 c2
for pair=1:size(pairwise_combinations,1)
    a=activity_bfF(:,pairwise_combinations(pair,1));
    b=activity_bfF(:,pairwise_combinations(pair,2));
    [c1(:,pair),lags]=xcov(a,b,100,'coeff');%cross-covariance
    pearson1(pair)=c1(101,pair);
end

for pair=1:size(pairwise_combinations,1)
    a=activity_aF(:,pairwise_combinations(pair,1));
    b=activity_aF(:,pairwise_combinations(pair,2));
    [c2(:,pair),lags]=xcov(a,b,100,'coeff');
    pearson2(pair)=c2(101,pair);
end

m1= mean(pearson1);
m2= mean(pearson2);

%f5 = figure;
%figure(f5) ;
%scatter(pearson1,pearson2);hold on; %shows relationship of correlation bf food and after (if high R before then also after)
%plot([min([pearson1 pearson2]) max([pearson1 pearson2])],[min([pearson1 pearson2]) max([pearson1 pearson2])]);
%xlabel('before food');
%ylabel('after food');
%scatter(m1,m2,'m', 'filled')

%boxplots comparing correlations before food and after food
pearson_bf_af(:,1) = pearson1;
pearson_bf_af(:,2) = pearson2;
clear ylabel
f5 = figure;
figure(f5) ;
coordLineStyle = 'k.';
boxplot(pearson_bf_af, 'Symbol', coordLineStyle, 'labels', {'before food', 'after food'});
ylabel('Correlation (r)')



[h1,p1] = ttest(pearson_bf_af(:,1),pearson_bf_af(:,2))
delta_correlations=mean(pearson_bf_af(:,2))-mean(pearson_bf_af(:,1));
%parallelcoords(pearson_bf_af, 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
 % 'Marker', '.', 'MarkerSize', 10);


%figure (f2);
%plot(pearson1);hold on;
%figure (f3);
%plot(pearson2);

%delta=pearson2-pearson1;
%pairwise_combinations(find(delta>0.5),:)

%% t-test
clear h p

activity_bfF=bgr_zsc(1:food_frame,:);
activity_aF=bgr_zsc(food_frame:end,:);
mean_bfF= mean(activity_bfF);
mean_aF= mean(activity_aF);

comp(:,1) = mean_bfF' ;
comp(:,2) = mean_aF' ;

clear ylabel
figure;
coordLineStyle = 'k.';
boxplot(comp(1:23,1:2), 'Symbol', coordLineStyle, 'labels', {'before food', 'after food'}); hold on;
parallelcoords(comp(1:23,1:2), 'Color', 0.7*[1 1 1], 'LineStyle', '-',...
  'Marker', '.', 'MarkerSize', 10);
ylabel('Mean activity (z-scored)')


[h2,p2]= ttest(comp(:,1),comp(:,2))% over all cells: overall activity bfF diff from aF?

%uneven samples?



