%%plotting of adlib data 
%5.11.23
%takes in Adlib_pellets and Adlib_drinks files
clearvars -except Adlib_pellets Adlib_drinks; close all;
load('ci_data_adlib_mf.mat');
%load('Adlib_pellets_bsl.mat');load('Adlib_drinks_bsl.mat');
%% remove drop stats from AdlibDri
for se = 1:size(Adlib_drinks,1)
    event = Adlib_drinks{se,2};
    ds = find(event == 'drop_stats');
    event(ds) = [];
    ts = Adlib_drinks{se,3};ts(ds)=[];
    Adlib_drinks{se,2} = event;
    Adlib_drinks{se,3} = ts;
end    
%% sort Adlib session into those with 3+ pellets
NrPel = cell2mat(Adlib_pellets(:,4));
adlib3 = NrPel >3;
AdlibPel = Adlib_pellets(adlib3,:); 

NrDrink = cell2mat(Adlib_drinks(:,4));
adlib3d = NrDrink >3;
AdlibDri = Adlib_drinks(adlib3d,:); 

%% separate out good cells
% collective neuropil cell is last cell in array and individual np is already subtracted

for i=1:numel(AdlibPel(:,1))
    if AdlibPel{i,5} =='g10'
        %good_cells=[1 2 4 5 8 10 11 13 14 15];% last is np
        good_cells=[1 2 5 8 15];
    elseif AdlibPel{i,5} =='g12'
        %good_cells=[1 2 3 4 5 6 7 8 9 14 15];
        good_cells=[2 5 15];
    elseif AdlibPel{i,5} =='g15' 
        %good_cells=[1 2 4 5 6 7 8 11 12 18 19];
        good_cells=[1 2 6 7 8 11 19];
    end
    snip = AdlibPel{i,1};snip=snip(:,good_cells);
    AdlibPel{i,1} = snip;
end


for i=1:numel(AdlibDri(:,1))
    if AdlibDri{i,5} =='g10'
        good_cells=[1 2 4 5 8 10 11 13 14];
    elseif AdlibDri{i,5} =='g12'
        good_cells=[1 2 3 4 5 6 7 8 9 14];
    elseif AdlibDri{i,5} =='g15'
        good_cells=[1 2 4 5 6 7 8 11 12 18];;    
    end
    snip = AdlibDri{i,1};snip=snip(:,good_cells);
    AdlibDri{i,1} = snip;
end
%% add new fluorescence across whole population as last "cell'
for ss = 1:size(AdlibPel,1)
    trace = AdlibPel{ss,1};
    trace_m = mean(trace(:,1:end-1),2);
    AdlibPel{ss,1}(:,end+1) = trace_m;
end

 %% try finding and subtracting baseline with matlab function msbackadj
% AdlibPel_raw = AdlibPel; %save original values, not raw bc df/f and zscored, but not baselined locally
% AdlibPel = AdlibPel_raw;
% for ps = 1:size(AdlibPel,1)
%     Meal = AdlibPel{ps,1};
%     X =[1:(size(Meal,1))]';
%     for ce = 1: size(Meal,2)
%         Cell = Meal(:,ce);
%         bsl = msbackadj(X,Cell,'WindowSize', 50,'StepSize', 4,'Quantile', 0,'ShowPlot',  'no');
%         Cell_bsl = Cell-(bsl*0.7);%subtract baseline by factor 0.7
%         Meal(:,ce) = Cell_bsl;
%     end
%     AdlibPel{ps,1} = Meal;
% end    

%% make simple plots of trace of each cell with marked pellet retrieval
figure;
session =16;
space=0;
for p= 1:size(AdlibPel{session,1},2)
   
    meal = AdlibPel{session,1};
    %meal(:,p) = meal(:,p) - meal(:,end);
    meal(:,p) = smooth(meal(:,p),6);% signal is smoothed!
    plot(meal(:,p)+ space, LineWidth=1.5);hold on
    space= space-4;
end
for line =1:size(AdlibPel{session,3},1)% plots lines for entry pellet and exit
    ev = AdlibPel{session,3};
    xline(ev(line), LineWidth=1 );
end    
%make a small axis label for traces
plot([100 50], [5 5],LineWidth= 3,Color='k');%[X X][Y Y]% 5s aka 50 frames
plot([50 50],[7 5],LineWidth= 3,Color= 'k'); %2s.d.
hold off;

%% traces of only one cell in one session

figure;
session =7;
space=0;
for p = 4%cell
   
    meal = AdlibPel{session,1};
    %meal(:,p) = normalize(meal(:,p));
    meal(:,p) = meal(:,p) - meal(:,end);
    meal(:,p) = smooth(meal(:,p),9);% signal is smoothed!
    
    plot(meal(:,p)+ space, LineWidth=1.5);hold on
    space= space-4;
end
for line =1:size(AdlibPel{session,3},1)% plots lines for entry pellet and exit
    ev = AdlibPel{session,3};
    xline(ev(line), LineWidth=1 );
end
%make a small axis label for traces
plot([100 50], [3 3],LineWidth= 3,Color='k');%[X X][Y Y]% 5s aka 50 frames
plot([50 50],[3 4],LineWidth= 3,Color= 'k'); %2s.d.
hold off;




%% make example for one cell response to first 2 pell midlle two and last two?

c22s16 = AdlibPel{16,1}(:,5);
c22s16BH = AdlibPel{16,3};
figure;
subplot(1,3,1); plot(smooth(c22s16(90:390),10),LineWidth=2); hold on;
    ylim([-1.5 3]);xlim([0 300]);
    ev= c22s16BH(3:7);%2, 3, 4, 5, 6 
    for p = 1:size(ev,1)
        xline(ev(p)-90, LineWidth=3, Color= 'r', Alpha=0.5 );
    end    

subplot(1,3,2); plot(smooth(c22s16(900:1200),10),LineWidth=2);
    ylim([-1.5 3]);xlim([0 300]);
    ev= c22s16BH(12:14);%pel 11 and 12, 13
    for p = 1:size(ev,1)
        xline(ev(p)-900, LineWidth=3, Color= 'r', Alpha=0.5 );
    end 

subplot(1,3,3); plot(smooth(c22s16(2000:2300),10),LineWidth=2);
    ylim([-1.5 3]);xlim([0 300]);
    ev= c22s16BH(21:22);%last 2 (20, 21)
    for p = 1:size(ev,1)
        xline(ev(p)-2000, LineWidth=3, Color= 'r', Alpha=0.5 );
    end 
hold off;
%% plot all in subplots
fA = figure;
figure(fA);
for session = 1: size(AdlibPel,1)-9
    space=0;
    for p= 1:(size(AdlibPel{session,1},2)-1)
        meal = AdlibPel{session,1};
        meal(:,p) = meal(:,p) - meal(:,end);%np subtracted!
        meal(:,p) = smooth(meal(:,p),5);% signal is smoothed! 
        h1 = subplot(5,2,session) ;plot(meal(:,p)+ space, LineWidth=1);hold on
        space= space-3;
    end
    for line =2:size(AdlibPel{session,3},1)-1% plots lines for entry pellet and exit
        ev = AdlibPel{session,3};
        subplot(5,2,session);xline(ev(line), LineWidth=1 );
    end   
end


fK = figure;
figure(fK);
for session = 10: size(AdlibPel,1)
    space=0;
    for p= 1:(size(AdlibPel{session,1},2)-1)
        meal = AdlibPel{session,1};
        meal(:,p) = meal(:,p) - meal(:,end);%np subtracted!
        meal(:,p) = smooth(meal(:,p),5);% signal is smoothed! 
        h1 = subplot(5,2,session-9) ;plot(meal(:,p)+ space, LineWidth=1);hold on
        space= space-3;
    end
    for line =2:size(AdlibPel{session,3},1)-1% plots lines for entry pellet and exit
        ev = AdlibPel{session,3};
        subplot(5,2,session-9);xline(ev(line), LineWidth=1 );
    end   
end








%% DRinks
% f2 = figure;
% figure(f2);
% session = 20;
% space=0;
% for p= 1:size(AdlibDri{session,1},2)
%     bout = AdlibDri{session,1};
% 
%     bout(:,p) = smooth(bout(:,p),5);% signal is smoothed!
%     plot(bout(:,p)+ space, LineWidth=1);hold on
%     space= space-3;
% end
% for line =1:size(AdlibDri{session,3},1)% plots lines for entry pellet and exit
%     ev = AdlibDri{session,3};
%     xline(ev(line), LineWidth=1 );
% end 

%% separate into separate animals
g10=[];g12=[];g15=[];
for el = 1:size(AdlibPel,1)
    label = AdlibPel{el,5};
    if label == "g10"
        a=AdlibPel(el,:);
        g10 = [g10; a];
    elseif label == "g12"
        a=AdlibPel(el,:);
        g12 = [g12; a];
    elseif label == "g15"
        a=AdlibPel(el,:);
        g15 = [g15; a];    
    end
end
Adlib.Pel.g10=g10;
Adlib.Pel.g12=g12;
Adlib.Pel.g15=g15;

%% snip baseline window and average across all trials to subtract from snip
% % later
% fn=fieldnames(Adlib);
% %loop through the fields
% for ii=1: numel(fn)%Pel or Drink
%     fn1=fieldnames(Adlib.(fn{ii}));
%     for j=1: numel(fn1)%animal
%         animal = Adlib.(fn{ii}).(fn1{j});
%         Bsnips = [];
%         for ses=1:size(animal,1)
%             win = 29; % 30 frames after exit
%             sn = 1;
%             snip = animal{ses,1}; snip= snip(end-win:end,:);
%             Bsnips(:,:,ses) = snip;
%         end
%         Bsnips_mean = mean(Bsnips,3);Bsnips_mean = mean(Bsnips_mean,1); %averg across session and across bsl window
%         BslSnips.(fn{ii}).(fn1{j}).bslMean = Bsnips_mean;
%         BslSnips.(fn{ii}).(fn1{j}).snips = Bsnips;
%     end
% end
% %concatenate and remove np and pop aver
% avergBslPel = cat(1,(BslSnips.Pel.g10.bslMean(1,1:end-2))',(BslSnips.Pel.g12.bslMean(1,1:end-2))',(BslSnips.Pel.g15.bslMean(1,1:end-2))');
% avergBslDri = cat(1,(BslSnips.Dri.g10.bslMean(1,1:end-2))',(BslSnips.Dri.g12.bslMean(1,1:end-2))',(BslSnips.Dri.g15.bslMean(1,1:end-2))');

%% collect snips around pellets and label which one in the meal they are
% determine snip size
clear s ts pe
TsDiff = [];
for s = 1:size(AdlibPel,1)
    ts = AdlibPel{s,3};ts = ts(2:end-1,:);
    for pe = 1:(size(ts,1)-1)
        diff = ts(pe+1)-ts(pe);
        TsDiff =[TsDiff,diff]; %collect all diff to sort later
    end    
end    
TS= sort(TsDiff);TS=TS(1,1); %this is the smallest frame difference between two pellet retrievals

fn=fieldnames(Adlib);
ColSnips ={};
ColInfo = {};
%loop through the fields
for ii=1: numel(fn)%Pel or Drink
    fn1=fieldnames(Adlib.(fn{ii}));
    for j=1: numel(fn1)%animal
        animal = Adlib.(fn{ii}).(fn1{j});
        for ses=1:size(animal,1)
            PelSnips = [];
            PelInfo = [];
            win = (TS-1)/2; %frames around pellet retrieval (picks the smallest frame difference between two retrievals) = 16s 8 on each side
            sn = 1;
            snip = animal{ses,1}; 
            snip = snip - snip(:,end); snip(:,end) = []; %subtract np and remove!!!
            animal{ses,1} = snip;
            for pel = 2:size(animal{ses,3},1)-1% starting from 2 bc 1 is entry and -1 bc last is exit
                evt= animal{ses,3};
                snip = animal{ses,1}; 
                snip= snip(evt(pel)-win:evt(pel)+win,:);
                PelSnips(:,:,sn) = snip;
                sn=sn+1;
                anim = animal{ses,5};
                info= {ses,pel-1,anim};
                PelInfo = vertcat(PelInfo, info);
            end
            ColSnips.(fn{ii}).(fn1{j}){ses} = PelSnips;
            ColInfo.(fn{ii}).(fn1{j}){ses} = PelInfo;
        end
    end
end
%% plot response amplitude vs pellet retrieval


fn=fieldnames(ColSnips.Pel);
%loop through the fields
    for j=1: numel(fn)%animal
        animal = ColSnips.Pel.(fn{j});
        for pel=1:size(animal,2)
            retr = animal{1,pel};
            m_pel_coll = [];
            for c = 1:size(retr,3)
                nr_pel = retr(:,:,c);
                m_pel = mean(nr_pel,1);
                m_pel_coll = vertcat(m_pel_coll, m_pel);
            end
            ColSnips.Pel.(fn{j}){2,pel}= m_pel_coll;
        end
    end

%% plot

fn=fieldnames(ColSnips.Pel);
for j = 1:size(fn,1)
    animal = ColSnips.Pel.(fn{j});
    for jj =1:size(animal,2)
        meal = ColSnips.Pel.(fn{j}){2,jj};
        for cell=1:size(meal,2)
            x = [1:size(meal(:,cell),1)];
            figure(s1);
            subplot(2,4,jj);plot(x,(meal(:,cell))',Color=colors(cell));hold on;      
        end
    end
end












s1=figure;
colors = ['r','g','b','c','m','y','k','#0072BD','#D95319','#EDB120','#7E2F8E',"#77AC30","#4DBEEE"];
%x=[1:30];
fn=fieldnames(ColSnips.Pel);
for j = 1:size(fn,1)
    animal = ColSnips.Pel.(fn{j});
    for jj =1:size(animal,2)
        meal = ColSnips.Pel.(fn{j}){2,jj};
        for cell=1:size(meal,2)
            x = [1:size(meal(:,cell),1)];
            figure(s1);
            subplot(2,4,jj);plot(x,(meal(:,cell))',Color=colors(cell));hold on;      
        end
    end
end













%%  PSTH

% concatenate over one animal
Pelg10 = [];
for s= 1:size(ColSnips.Pel.g10,2)
    sess = ColSnips.Pel.g10{1,s};
    sess = sess(:,1:end-1,:);%remove averg across population
    Pelg10  = cat(3,Pelg10, sess);
end
Pelg10m = mean(Pelg10,3);%average across all cells for one animal
%imagesc(Pelg10m');

%g12
Pelg12 = [];
for s= 1:size(ColSnips.Pel.g12,2)
    sess = ColSnips.Pel.g12{1,s};
    sess = sess(:,1:end-1,:);
    Pelg12  = cat(3,Pelg12, sess);
end
Pelg12m = mean(Pelg12,3);%average across all cells for one animal
%imagesc(Pelg12m');

%g15
Pelg15 = [];
for s= 1:size(ColSnips.Pel.g15,2)
    sess = ColSnips.Pel.g15{1,s};
    sess = sess(:,1:end-1,:);
    Pelg15  = cat(3,Pelg15, sess);
end

Pelg15m = mean(Pelg15,3);%average across all cells for one animal
%imagesc(Pelg15m');


Pel_averg = (cat(2,Pelg10m(:,1:end-1,:),Pelg12m(:,1:end-1,:),Pelg15m(:,1:end-1,:)))';%without np
figure;
imagesc(Pel_averg);hold on
xline(win+1,LineWidth=2);
xticks([4 9 14]); xticklabels({'-0.5','0' ,'0.5' });
xlabel('\fontsize{12}time from pellet retrieval [s]'); ylabel('Cells');

%sorted 
Pelaverg_mean = mean(Pel_averg(:,7:11),2);Pelaverg_mean(:,2)=[1:26];
Pelaverg_mean_sort = sortrows(Pelaverg_mean);

for r = 1:size(Pelaverg_mean_sort,1)
    k = Pelaverg_mean_sort(r,2);
    Pel_averg_sort(r,:) = Pel_averg(k,:) 
end
f2 = figure; 
figure(f2);
imagesc(Pel_averg_sort);hold on
xline(win+1,LineWidth=2);
xticks([4 9 14]); xticklabels({'-0.5','0' ,'0.5' });
xlabel('\fontsize{12}time from pellet retrieval [s]'); ylabel('Cells');
%% std shade [xLeft, yBottom, width, height]
clear subplot
f9=figure;
figure(f9);hold on;
sp1 = subplot(2,1,1);imagesc(Pel_averg_sort);set(sp1, 'Position', [0.13,0.58,0.3,0.37]);
    ylim([1 26]); xlim([1 17]);
    xline(win+1,LineWidth=1,Color='w',LineStyle='--');
    xticks([4 9 14]); xticklabels({'','' ,'' });
    ylabel('Cells');
    %yline(8.5,LineWidth=1,Color='k',LineStyle='-');
    %yline(17.5,LineWidth=1,Color='k',LineStyle='-');

sp2= subplot(2,1,2);stdshade(Pel_averg_sort,0.3,[0 0 1]);set(sp2, 'Position', [0.13,0.475,0.3,0.10]);
    xticks([4 9 14]); xticklabels({'-0.5','0' ,'0.5' });
    xlabel('\fontsize{12}time from pellet retrieval [s]');
    ylabel('Ca+-response (z-scored, s.d.)');
    ylim([-0.5 0.5]);xlim([1 17]);
    yticks([-0.3 0 0.3]);yticklabels({'-0.3','0' ,'0.3' });
    %xline(win+1,LineWidth=1,Color='w',LineStyle='--');
hold off;


%% find peaks and classify cells into preceding and following
clear cell C MaxI Cmax
for cell = 1:size(Pel_averg,1)
    C = Pel_averg(cell,:);
    [Cmax(1,cell),MaxI(1,cell)] = max(C);%get largest pos or neg value
    Cmax(1,cell) = C(MaxI(1,cell));%get the actual not absolute value
    if Cmax(cell)<0
        MaxI(2,cell) = 1;
    elseif Cmax(cell)>=0   
        MaxI(2,cell) = 0;
    end
end 
% classify peaks as positive and negative for plotting
% U(1,:) = 1:17;
% U(2,:) = zeros(1,17);%collect nr of negative
% U(3,:) = zeros(1,17);%collect nr of positive at that timepoint
% for un = 1:size(U,2)
%     f = find(MaxI(1,:) == U(1,un));
%     for t = 1:(size(f,2))
%         if MaxI(2,f(t)) ==1
%            U(2,un)= U(2,un)+1;
%         elseif MaxI(2,f(t)) == 0
%            U(3,un)= U(3,un)+1 ;
%         end
%     end
% end    
% %plot as stacked bars
% pea = U(2:3,:);
% figure; bar(pea', 'stacked');
% xticks([2 3 4 6 7 10 11 12 13 14 15 16 17 ]); xticklabels({'-.7','-.6','-.5','-.3','-.2','.1','.2','.3','.4','.5','.6','.7','.8'});
% xlabel('\fontsize{12}peak from pellet retrieval [s]'); ylabel('Nr of Cells');
% ylim([0 6]);


%plot histogram of cells
figure;
histogram(MaxI(1,:)); hold on;
ylim([0 9]);
xticks([1 2 3 4 5 6 7 10 11 13 14 15 17 ]); xticklabels({'-.8','-.7','-.6','-.5','-.4','-.3','-.2','.1','.2','.4','.5','.6','.8',});
xlabel('\fontsize{12}peak from pellet retrieval [s]'); ylabel('Nr of Cells');

figure;
histogram(MaxI(1,:),2);hold on;
ylim([0 16]);
xticks([4 14]); xticklabels({'peak pre pellet ','peak post pellet'});
ylabel('Nr of Cells');
% 
% nr_pre=find(MaxI < 9);
% nr_post= find (MaxI > 9);
%% classify cells as On or Off
% all cells marked with 1 in row 2 of Classif
% Classif = [1:size(Pel_averg,1)];
% clear cell
% for cell = 1:size(Pel_averg,1)
%     Cellmean(cell) = mean(Pel_averg(cell,:),2);
%     if Cmax(cell) > 0.4;%on cells
%        Classif(2,cell) = 1 ; 
%     elseif Cmax(cell) < -0.4% off cells ARBITRARY Threshold !!!!!!!!!!!!!!!!!!!!!
%        Classif(2,cell) = 2 ;   
%     else
%        Classif(2,cell) = 0 ;   
%     end
% end
%% plot On and off cells
% OnCells = find(Classif(2,:)==1);%index of on cells
% OffCells = find(Classif(2,:)==2);%index of on cells
% NoneCells = find(Classif(2,:)==0);
% Pel_On =Pel_averg([OnCells],:);
% Pel_Off =Pel_averg([OffCells],:);
% Pel_None =Pel_averg([NoneCells],:);
% 
% CellsClassif = cat(1,Pel_On,Pel_None,Pel_Off); 
% figure;
% p1=subplot(1,1,1); imagesc(CellsClassif);hold on;set(p1, 'Position', [0.13,0.58,0.3,0.37]);
%     ylim([1 26]); xlim([1 17]);
%     xline(win+1,LineWidth=1,Color='w',LineStyle='--');
%     xticks([4 9 14]); xticklabels({'-.5','0' ,'.5' });
%     xlabel('\fontsize{12}time from pellet retrieval [s]');
%     yline(size(Pel_On,1)+0.5,LineWidth=1,Color='w',LineStyle='-');
%     yline(size(Pel_None,1)+size(Pel_On,1)+0.5,LineWidth=1,Color='w',LineStyle='-');
%     ylabel('Cells');
% 
% 
% C = Pel_On(1,:);
% [Cmax,MaxI] = max(C);
% clear cell C
% for cell = 1:size(Pel_On,1)
%     C = Pel_On(cell,:);
%     [Cmax(cell),MaxI(cell)] = max(C);
% end  

%plot histogram of cells
% figure;
% histogram(MaxI); hold on;
% ylim([0 4]);
% xticks([ 2  4  6 7 10 11 14 15 17 ]); xticklabels({'-.8','-.5','-.3','-.2','.1','.2','.5','.6','.8',});
% yticks([1 2 3 4]);
% xlabel('\fontsize{12}peak from pellet retrieval [s]'); ylabel('Nr of Cells');
% 
% figure;
% histogram(MaxI,2);hold on;
% ylim([0 10]);
% xticks([4 14]); xticklabels({'peak pre pellet ','peak post pellet'});
% ylabel('Nr of Cells');
%% plot histogram/stacked bar of on and off cells
% clear cell C MaxI Cmax
% for cell = 1:size(Pel_averg,1)
%     C = Pel_averg(cell,:);
%     [Cmax(1,cell),MaxI(1,cell)] = max(abs(C));%get largest pos or neg value
%     Cmax(1,cell) = C(MaxI(1,cell));%get the actual not absolute value
%     if Cmax(cell)<0
%         MaxI(2,cell) = 1;%neagtive peak
%     elseif Cmax(cell)>=0   
%         MaxI(2,cell) = 0;%positive peak
%     end
% end 
% % remove NOne cells from MaxI
% MaxI(:,NoneCells) = [];
% % classify peaks as positive and negative for plotting
% clear U;
% U(1,:) = 1:17;
% U(2,:) = zeros(1,17);%collect nr of negative
% U(3,:) = zeros(1,17);%collect nr of positive at that timepoint
% for un = 1:size(U,2)
%     f = find(MaxI(1,:) == U(1,un));
%     for t = 1:(size(f,2))
%         if MaxI(2,f(t)) ==1
%            U(2,un)= U(2,un)+1;
%         elseif MaxI(2,f(t)) == 0
%            U(3,un)= U(3,un)+1 ;
%         end
%     end
% end    
% %plot as stacked bars
% pea = U(2:3,:);
% figure; bar(pea', 'stacked');
% xticks([3 6 10 11 12 13 14 16 17 ]); xticklabels({'-.6','-.3','.1','.2','.3','.4','.5','.7','.8'});
% xlabel('\fontsize{12}peak from pellet retrieval [s]'); ylabel('Nr of Cells');
% ylim([0 6]);



%% None cells could be onlz responding to first pellet or entry?



%% plot individual trials of ON cells?
%cell 5 (22) g15
hold off;
cell5g15 = Pelg15(:,5,:);
cell5= squeeze(cell5g15);
figure; imagesc(cell5');hold on;
xline(9,LineWidth=2);
title(['Cell 22']);
xticks([4 9 14]); xticklabels({'-.5','0' ,'.5' });
xlabel('\fontsize{12}time from pellet retrieval [s]');
ylabel('Trials');

%cell 4 g10
hold off;
cell4g10 = Pelg10(:,4,:);
cell4= squeeze(cell4g10);
figure; imagesc(cell4');hold on;
xline(9,LineWidth=2);
title(['Cell 4']);
xticks([4 9 14]); xticklabels({'-.5','0' ,'.5' });
xlabel('\fontsize{12}time from pellet retrieval [s]');
ylabel('Trials');

%cell 1 g10 OFF
hold off;
cell1g10 = Pelg10(:,1,:);
cell1= squeeze(cell1g10);
figure; imagesc(cell1');hold on;
xline(9,LineWidth=2);
title(['Cell 1']);
xticks([4 9 14]); xticklabels({'-.5','0' ,'.5' });
xlabel('\fontsize{12}time from pellet retrieval [s]');
ylabel('Trials');



%% cross correlation
% make binary behavior vector of 0/1 for pellet retrieval
% cross correlate with imaging signal with lag of 20 frames?

% make binary behavior array
PelCorr = AdlibPel;
for ses = 1:size(PelCorr,1)%loop through sessions
    Session = PelCorr{ses,1};
    binaryBh = zeros(size(Session,1),1);
    Peltimes = PelCorr{ses,3}(2:end-1,:); %fills pellet timestamps with 1
    binaryBh(Peltimes)=1;
    PelCorr{ses,6} = binaryBh;
end    

PelCorr(:,[2 3 4])=[];

%concatenate cells activity together
Catg10 = cat(1,PelCorr{1:8,1});Catg10(:,[end-1:end])=[] ;Bhg10 = cat(1,PelCorr{1:8,3});
Catg12 = cat(1,PelCorr{9:15,1});Catg12(:,[end-1:end])=[] ;Bhg12 =cat(1,PelCorr{9:15,3});
Catg15 = cat(1,PelCorr{16:19,1});Catg15(:,[end-1:end])=[] ;Bhg15 =cat(1,PelCorr{16:19,3});
CorrSess.g10{1,1} = Catg10; CorrSess.g10{1,2} = Bhg10;
CorrSess.g12{1,1} = Catg12; CorrSess.g12{1,2} = Bhg12;
CorrSess.g15{1,1} = Catg15; CorrSess.g15{1,2} = Bhg15;
% separate out only ON cells for xcorr
% OnCorrSess = CorrSess;
% OnCorrSess.g10{1,1} = Catg10(:,[1 2 4]);
% OnCorrSess.g12{1,1} = Catg12(:,[9 16-9 ]);
% OnCorrSess.g15{1,1} = Catg15(:,[18-17 18-17 21-17 22-17 24-17 25-17]);


%% loop though all cells to get xcorr and p values of r value
clear ce an behavior fluo cell
PValue=[];
RValue=[];
MaxLags=[];
AllCs=[];
lag = 15;
fn=fieldnames(CorrSess);%fn=fieldnames(OnCorrSess);
for an = 1:size(fn,1)%loop through animals
    behavior = CorrSess.(fn{an}){1,2};
    fluo = CorrSess.(fn{an}){1,1};
    for ce =1: size(fluo,2)%loop through cells of that animal
        cell = fluo(:,ce);
        [c, lags] = xcorr(cell,behavior,lag, 'normalized');%figure;plot(lags,c);
        AllCs=[AllCs, c];% collects cross-correlogram
        [~, ii] = max(abs(c));
        maxLag = lags(ii);
        MaxLags=[MaxLags, maxLag];
        % Shift one of the signals by the maxLag
        shiftedbh = circshift(behavior, maxLag);%shifts f data by maxshift
        [r, p] = corrcoef(cell, shiftedbh);
        PValue = [PValue,p(1,2)];
        RValue = [RValue,r(1,2)];
    end
end
%% plot individual correlograms
figure;
p1= plot(AllCs(:,4));hold on;
xticks([1 6 11 16 21 26 31]);xticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'});
xlim([1 31]);
xlabel('\fontsize{12}lag [s]');
ylabel('R-Value');
xline(16);yline(0);
p1.LineWidth = 2;
title('Cell 4');
hold off;

% plot averg correlogramm
meanxCor = mean(AllCs,2);
figure;
%plot(AllCs,'b');
stdshade(AllCs',0.2,[0.9290 0.6940 0.1250]); hold on;
%plot(meanxCor,LineWidth=3,Color='b');
xticks([1 6 11 16 21 26 31]);xticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'});
xlim([1 31]);ylim([-0.09 0.09]);
xlabel('\fontsize{12}lag [s]');
ylabel('XCorr');
xline(16);yline(0);hold off;




% separate cells according to PValue
sigC=find(PValue <= 0.05);
sigLags=MaxLags(sigC);%lags of sigC
Mprelags=mean(sigLags(find(sigLags < 0)));
Mpostlags=mean(sigLags(find(sigLags > 0)));

preCells =sigC(find(sigLags < 0)); PreC=AllCs(:,preCells);
postCells =sigC(find(sigLags > 0));PostC= AllCs(:,postCells);
%plot average xcorr of pre and post cells
figure; 
%plot(PreC,'b');hold on;
%plot(PostC,'m');
%plot(mean(PreC,2),'b','LineWidth',2);
%plot(mean(PostC,2),'m','LineWidth',2);
stdshade(PreC',0.2,[0.4940 0.1840 0.5560]);hold on;
stdshade(PostC',0.2,[0.8500 0.3250 0.0980]);
xticks([1 6 11 16 21 26 31]);xticklabels({'-1.5','-1','-0.5','0','0.5','1','1.5'});
xlim([1 31]);ylim([-0.09 0.09]);
xlabel('\fontsize{12}lag [s]');
ylabel('XCorr');
xline(16);yline(0);
% lag from 0 for averg sign corr cells?
%mPre=mean(PreC,2); mPreMax=find(mPre==(maxabs(mPre)));
%mPost=mean(PostC,2);mPostMax=find(mPost==(max(mPost)))

xline(11,'--g',LineWidth=2);xline(23,'--g',LineWidth=2);


hold off;

%% make imagesc for sign xcorr cells
signCells = cat(2,preCells,postCells);
NoSignCells= Pel_averg; NoSignCells(signCells,:) = [];
pre = Pel_averg(preCells,:);
post = Pel_averg(postCells,:);
Pel_averg_xcor= cat(1,NoSignCells,pre,post);

clear subplot
f9=figure;
figure(f9);hold on;
sp1 = subplot(2,1,1);imagesc(Pel_averg_xcor);set(sp1, 'Position', [0.13,0.58,0.3,0.37]);
    ylim([1 26]); xlim([1 17]);
    xline(win+1,LineWidth=1,Color='w',LineStyle='--');
    xticks([4 9 14]); xticklabels({'','' ,'' });
    ylabel('Cells');
    yline(6.5,LineWidth=1,Color='k',LineStyle='-');
    yline(11.5,LineWidth=1,Color='k',LineStyle='-');

sp2= subplot(2,1,2);set(sp2, 'Position', [0.13,0.475,0.3,0.10]);
    stdshade(pre,0.2,[0.4940 0.1840 0.5560]);hold on; 
    stdshade(post,0.2,[0.8500 0.3250 0.0980]);
    xticks([4 9 14]); xticklabels({'-0.5','0' ,'0.5' });
    xlabel('\fontsize{12}time from pellet retrieval [s]');
    ylabel('Ca+-response (z-scored, s.d.)');
    ylim([-0.5 0.5]);xlim([1 17]);
    yticks([-0.3 0 0.3]);yticklabels({'-0.3','0' ,'0.3' });
    %xline(win+1,LineWidth=1,Color='w',LineStyle='--');
hold off;



%% 

%plot all z in nice way?
%scatterplot
Y=1:(size(PValue,2));
corrsize = find(PValue <= 0.05);CorSize =zeros((size(PValue,2)),1);CorSize(corrsize)=30;CorSize=CorSize+30;

figure; scatter(MaxLags,Y,CorSize,RValue,'filled');hold on;%color to mark r value and shape for significance?
xlim([-16 16]);ylim([0 (size(PValue,2))+1]);
xticklabels({'-1.5','-1','-.5','0', '.5', '1','1.5'});xlabel('\fontsize{12}maximum lag [s]')
ylabel('\fontsize{12}Cells');colormap("copper");clim([-0.15 0.15]);
hold off;


%% enter food pod cells xcor
% make binary behavior vector of 0/1 for entry into food pod
% cross correlate with imaging signal with lag of 20 frames?

% make binary behavior array
PelCorr = AdlibPel;
for ses = 1:size(PelCorr,1)%loop through sessions
    Session = PelCorr{ses,1};
    binaryBh = zeros(size(Session,1),1);
    Peltimes = PelCorr{ses,3}(1,:); %fills ef timestamps with 1
    binaryBh(Peltimes)=1;
    PelCorr{ses,6} = binaryBh;
end    

PelCorr(:,[2 3 4])=[];

%concatenate cells activity together
Catg10 = cat(1,PelCorr{1:8,1});Catg10(:,[end-1:end])=[] ;Bhg10 = cat(1,PelCorr{1:8,3});
Catg12 = cat(1,PelCorr{9:15,1});Catg12(:,[end-1:end])=[] ;Bhg12 =cat(1,PelCorr{9:15,3});
Catg15 = cat(1,PelCorr{16:19,1});Catg15(:,[end-1:end])=[] ;Bhg15 =cat(1,PelCorr{16:19,3});
CorrSess.g10{1,1} = Catg10; CorrSess.g10{1,2} = Bhg10;
CorrSess.g12{1,1} = Catg12; CorrSess.g12{1,2} = Bhg12;
CorrSess.g15{1,1} = Catg15; CorrSess.g15{1,2} = Bhg15;
%% loop though all cells to get xcorr and p values of r value
clear ce an behavior fluo cell
PValue=[];
RValue=[];
MaxLags=[];
AllCs=[];
lag = 15;
fn=fieldnames(CorrSess);
for an = 1:size(fn,1)%loop through animals
    behavior = CorrSess.(fn{an}){1,2};
    fluo = CorrSess.(fn{an}){1,1};
    for ce =1: size(fluo,2)%loop through cells of that animal
        cell = fluo(:,ce);
        [c, lags] = xcorr(cell,behavior,lag,'normalized');%figure;plot(lags,c);
        AllCs=[AllCs, c];
        [~, ii] = max(abs(c));
        %maxIndex = find(c == max(abs(c)));
        maxLag = lags(ii);
        MaxLags=[MaxLags, maxLag];
        % Shift one of the signals by the maxLag
        shiftedbh = circshift(behavior, maxLag);%shifts f data by maxshift
        [r, p] = corrcoef(cell, shiftedbh);
        PValue = [PValue,p(1,2)];
        RValue = [RValue,r(1,2)];
    end
end

%plot in nice way?
%scatterplot
Y=1:(size(PValue,2));
corrsize = find(PValue <= 0.05);CorSize =zeros((size(PValue,2)),1);CorSize(corrsize)=30;CorSize=CorSize+30;
figure; scatter(MaxLags,Y,CorSize,RValue,'filled');hold on;%color to mark r value and shape for significance?
xlim([-16 16]);ylim([0 (size(PValue,2))+1]);
xticklabels({'-1.5','-1','-.5','0', '.5', '1','1.5'});xlabel('\fontsize{12}maximum lag [s]')
ylabel('\fontsize{12}Cells');hold off;

%% collect snips around entries 
clear s ts pe
TS = 41;
fn=fieldnames(Adlib);
%loop through the fields
for ii=1: numel(fn)%Pel or Drink
    fn1=fieldnames(Adlib.(fn{ii}));
    for j=1: numel(fn1)%animal
        animal = Adlib.(fn{ii}).(fn1{j});
        for ses=1:size(animal,1)
            PelSnips = [];
            PelInfo = [];
            win = (TS-1)/2; % 20 frames around pellet retrieval
            sn = 1;
            for pel = 1% starting from 2 bc 1 is entry and -1 bc last is exit
                evt= animal{ses,3};
                snip = animal{ses,1}; snip= snip(evt(pel)-win:evt(pel)+win,:);
                PelSnips(:,:,sn) = snip;
                sn=sn+1;
                anim = animal{ses,5};
                info= {ses,pel-1,anim};
                PelInfo = vertcat(PelInfo, info);
            end
            ColSnips.(fn{ii}).(fn1{j}){ses} = PelSnips;
            ColInfo.(fn{ii}).(fn1{j}){ses} = PelInfo;
        end
    end
end

%%  PSTH

% concatenate over one animal
Pelg10 = [];
for s= 1:size(ColSnips.Pel.g10,2)
    sess = ColSnips.Pel.g10{1,s};
    sess = sess(:,1:end-1,:);%remove averg across population
    Pelg10  = cat(3,Pelg10, sess);
end
Pelg10m = mean(Pelg10,3);%average across all cells for one animal
%imagesc(Pelg10m');

%g12
Pelg12 = [];
for s= 1:size(ColSnips.Pel.g12,2)
    sess = ColSnips.Pel.g12{1,s};
    sess = sess(:,1:end-1,:);
    Pelg12  = cat(3,Pelg12, sess);
end
Pelg12m = mean(Pelg12,3);%average across all cells for one animal
%imagesc(Pelg12m');

%g15
Pelg15 = [];
for s= 1:size(ColSnips.Pel.g15,2)
    sess = ColSnips.Pel.g15{1,s};
    sess = sess(:,1:end-1,:);
    Pelg15  = cat(3,Pelg15, sess);
end

Pelg15m = mean(Pelg15,3);%average across all cells for one animal
%imagesc(Pelg15m');


Pel_averg = (cat(2,Pelg10m(:,1:end-1,:),Pelg12m(:,1:end-1,:),Pelg15m(:,1:end-1,:)))';%without np
figure;
imagesc(Pel_averg);hold on
xline(win+1,LineWidth=2);
%xticks([4 9 14]); xticklabels({'-0.5','0' ,'0.5' });
xlabel('\fontsize{12}time from pellet retrieval [s]'); ylabel('Cells');

%sorted 
Pelaverg_mean = mean(Pel_averg(:,7:11),2);Pelaverg_mean(:,2)=[1:26];
Pelaverg_mean_sort = sortrows(Pelaverg_mean);

for r = 1:size(Pelaverg_mean_sort,1)
    k = Pelaverg_mean_sort(r,2);
    Pel_averg_sort(r,:) = Pel_averg(k,:) 
end
f2 = figure; 
figure(f2);
imagesc(Pel_averg_sort);hold on
xline(win+1,LineWidth=2);
xticks([4 9 14]); xticklabels({'-0.5','0' ,'0.5' });
xlabel('\fontsize{12}time from pellet retrieval [s]'); ylabel('Cells');

%plot(Pel_averg');
%plot(Pel_averg_bsl')
%% ethogram plots
%good_cells=[1 2 4 5 8 10 11 13 14];animal ='g10';
good_cells=[1 2 3 4 5 6 7 8 9 14]; animal ='g12';
%good_cells=[1 2 4 5 6 7 8 11 12 18];animal ='g15';

factor = 0.9;
fn=fieldnames(ci_data.adlib);
ii=2;%animal
fn1=fieldnames(ci_data.adlib.(fn{ii}));
j=1 ; s=2;%day and session
fn2=fieldnames(ci_data.adlib.(fn{ii}).(fn1{j}));
k= 1;
fn3= fieldnames(ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}));
                
 event_type = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(4);
 raw_f_res=ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
 Fneu = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2});
 subneuro_trace = subneuro(raw_f_res,Fneu,factor); % subtracts neuropil trace from raw trace by factor x subneuro(rawF, Fneuropil, factor)
 subneuro_trace(:,(size(subneuro_trace,2)+1)) = mean(Fneu,2);%last "cell" is mean neuropil trace
                
                
 df_f_trace = dFoF(subneuro_trace);%df/f
 raw_f= zscore(df_f_trace,0,1);% zscored using sample sd
 raw_f = raw_f(:,good_cells);% pick out good cells and np(last cell in array) % second to last is averg but across all not just good cells
 mean_raw_f = mean(raw_f(:,1:end-1),2);
 raw_f= [raw_f,mean_raw_f];%averg across selected pop is last one

x = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(5);%frames

y=NaN(size(x));
frameplot=raw_f;
%y(find(event_type=='block_start'))=2;
y(find(event_type=='imaging_start'))=2; 
y(find(event_type=='enter_feed'))=4; %food entry
y(find(event_type=='enter_run'))=3; %wheel entry
y(find(event_type=='enter_social'))=4.5; %
y(find(event_type=='enter_drink'))=2.5; %
y(find(event_type=='enter_explore'))=5; 
%consume
y(find(event_type=='retrieve_pellet'))=4;%pellet retrieval
y(find(event_type=='drink'))=2.5; %
y(find(event_type=='run'))=3; %wheel entry
pel=x(find(event_type=='retrieve_pellet'));
drink=x(find(event_type=='drink'));
run=x(find(event_type=='run'));
%exit
y(find(event_type=='block_end'))=2;
y(find(event_type=='imaging_stop'))=2;
y(find(event_type=='exit_feed'))=3.5;
y(find(event_type=='exit_drink'))=3.5;
y(find(event_type=='exit_run'))=3.5;
y(find(event_type=='exit_social'))=3.5;
y(find(event_type=='exit_explore'))=3.5;

yy=fillmissing(y,'linear');

figure;a=1;
subplot(1,1,a),plot([x(1) x(end)],[3.5 3.5],'b');hold on%decisionpoint
%subplot(1,1,a),plot([x(1) x(end)],[2 2],'b');%nest
subplot(1,1,a),plot(x,yy,'k','LineWidth',1);hold on
subplot(1,1,a),plot(x,y,'k.','MarkerSize', 20);hold on

title(animal);
yticks([2 2.5 3 3.5 4 4.5 5]);
yticklabels({'nest','drink','run','decision point','food','social','explore'});
ylim([-12 5]);xlim([-10 size(raw_f,1)+10]);
xticks([0 6000 12000 18000 24000 30000]);
xticklabels({'0','10','20','30','40', '50'});
xlabel("time [min]")

%plot consumption
for i=1:1:length(pel)
    subplot(1,1,a),plot([pel(i) pel(i)],[3.8 4.2],'y','LineWidth',3);
end
for i=1:1:length(drink)
    subplot(1,1,a),plot([drink(i) drink(i)],[2.3 2.7],'c','LineWidth',3);
end
for i=1:1:length(run)
    subplot(1,1,a),plot([run(i) run(i)],[2.8 3.2],'m','LineWidth',3);
end

j=0;
for i=1:size(raw_f,2)
   subplot(1,1,a),plot(smooth(raw_f(4:end,i),5)/6-j,'Color',[i/size(raw_f,2) 1-i/size(raw_f,2) i/size(raw_f,2)],'LineWidth',1);hold on;
   j=j+1;
end

%subplot(1,1,a),plot([300 600],[1 1],'k','LineWidth',2);%30s
%subplot(1,1,a),plot([300 300],[0.98 0.98+(2/6)],'k','LineWidth',2);%2Std
subplot(1,1,a),plot([3400 3700],[1 1],'k','LineWidth',2);%30s
subplot(1,1,a),plot([3400 3400],[0.98 0.98+(2/6)],'k','LineWidth',2);%2Std
ylabel("Ca-activity (z-scored)");




