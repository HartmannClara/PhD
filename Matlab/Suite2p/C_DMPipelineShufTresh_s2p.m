clearvars -except events_averg ci_data; close all;
set(0,'defaultAxesFontSize',10);
%% for shuffled data
fn=fieldnames(events_averg); Iter = 100;
for ii=1:size(fn,1)
    fn1=fieldnames(events_averg.(fn{ii})); 
        for j=20:numel(fn1) % only process baselined data
            % divide trials by then nr of iterations 
            start = 1;
            Eend = (size(events_averg.(fn{ii}).(fn1{j}),3)/Iter);
            for k =1:Iter
            events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr = events_averg.(fn{ii}).(fn1{j})(:,:,start:Eend);
            start = start + (size(events_averg.(fn{ii}).(fn1{j}),3)/Iter);
            Eend = Eend + (size(events_averg.(fn{ii}).(fn1{j}),3)/Iter);
            end
        end
end

% duplicate events which are analyzed in multiple ways
fn = fieldnames(events_averg.shuf);
for ii = 1:size(fn,1)
    events_averg.shuf.(fn{ii}).be_averg = events_averg.shuf.(fn{ii}).blockend_averg_bsl;
    events_averg.shuf.(fn{ii}).home_averg =events_averg.shuf.(fn{ii}).blockend_averg_bsl;
    events_averg.shuf.(fn{ii}).drink_antic = events_averg.shuf.(fn{ii}).drink_averg_bsl;
    events_averg.shuf.(fn{ii}).eat_antic = events_averg.shuf.(fn{ii}).eat_averg_bsl;
    events_averg.shuf.(fn{ii}) = rmfield(events_averg.shuf.(fn{ii}), 'blockend_averg_bsl');% remove blocked field
end    
events_averg1= events_averg;
%make averages
%% makes an average over all trials per cell
fn=fieldnames(events_averg.shuf);
%loop through the fields to make mean all bsl trials
for ii=1: numel(fn)
    fn1=fieldnames(events_averg.shuf.(fn{ii}));
        for j=1:numel(fn1)
            %different analysis windows for 
            arr = events_averg.shuf.(fn{ii}).(fn1{j}).Nr;empt = isempty(arr); % if array is empty move on
            if empt == true
                j=j+1;
            elseif j == 4 || j ==15%soc expl
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(10:25,:),1);% takes 1.5 s window 1s after beam break,most likelihood 
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            elseif j == 7 || j == 8 %drink eat
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(31:45,:),1);% 1.5s from beam break 
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            elseif   j ==20%home
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(30:89,:),1);% 6s window 3 s after exit beam
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            elseif  j == 21 || j == 22% antic
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(16:30,:),1);%1.5s right before beam break
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            else % all the other ones 
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(16:30,:),1);
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            end 
        end
end 
%manually delete run

fn=fieldnames(events_averg.shuf);
%combine shuffled means in one array
for ii=1:numel(fn1)
    fn1=fieldnames(events_averg.shuf.(fn{ii}));
        for j=1:numel(fn1)
            arr = events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr;empt = isempty(arr); % if array is empty move on
            if empt == true
                j=j+1;
            else    
                shufAnimal = [];
                for k = 1:Iter
                    shufAnimal(k,:) = events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean;
                end
                events_averg.shufCells.(fn{ii}).(fn1{j}) = shufAnimal;
            end
        end    
end

% combine all cells together
fn1=fieldnames(events_averg.shufCells.g5);
    for j=1:numel(fn1)
        events_averg.shufCells.all.(fn1{j}) = cat(2,events_averg.shufCells.g5.(fn1{j}),events_averg.shufCells.g4.(fn1{j}),events_averg.shufCells.g2.(fn1{j}),events_averg.shufCells.g12.(fn1{j}));
    end 
% get 5 and 95 percentile values
fn = fieldnames(events_averg.shufCells.all);
for ii = 1:numel(fn)
    for k = 1:115
       perc95(ii,k) = prctile(events_averg.shufCells.all.(fn{ii})(:,k), 95, 1);
       perc5(ii,k) =  prctile(events_averg.shufCells.all.(fn{ii})(:,k), 5, 1);
    end   
end    
cl =[1:115];
Perc5 = array2table(perc5,"RowNames",fn); %This only contains the event where there is trials from each animal (no run er or exr)
Perc95= array2table(perc95,"RowNames",fn);

% only animals with run cells
fn1=fieldnames(events_averg.shufCells.g5);
    for j=1:numel(fn1)
        events_averg.shufCells.all.(fn1{j}) = cat(2,events_averg.shufCells.g2.(fn1{j}),events_averg.shufCells.g12.(fn1{j}));
    end 
% get 5 and 95 percentile values 
fn = fieldnames(events_averg.shufCells.all);
for ii = 1:numel(fn)
    for k = 6
       perc95(ii,1) = prctile(events_averg.shufCells.all.(fn{ii})(:,k), 95, 1);
       perc5(ii,1) =  prctile(events_averg.shufCells.all.(fn{ii})(:,k), 5, 1);
    end
    for k = 6
       perc95(ii,1) = prctile(events_averg.shufCells.all.(fn{ii})(:,k), 95, 1);
       perc5(ii,1) =  prctile(events_averg.shufCells.all.(fn{ii})(:,k), 5, 1);
    end 
end    
cl =[1:66];


Perc5_2 = array2table(perc5,"RowNames",fn);
Perc95_2= array2table(perc95,"RowNames",fn);


Perc5 = array2table(perc5,"RowNames",fn);
Perc95= array2table(perc95,"RowNames",fn);


%% makes an average over all trials per cell
fn=fieldnames(events_averg);me = {'_m'};
%loop through the fields to make mean all bsl trials
for ii=1: numel(fn)
    fn1=fieldnames(events_averg.(fn{ii}));
        for j=17:numel(fn1)
        events_averg.(fn{ii}).(char(strcat(fn1{j},me))) = nanmean(events_averg.(fn{ii}).(fn1{j}),3);
        end
end    

%% concatenating over averg from all cells (all cells together)
fn1=fieldnames(events_averg.g5);
    for j=37:numel(fn1)
        events_averg.all.(strrep((fn1{j}), '_averg_bsl_m', '')) = cat(2,events_averg.g5.(fn1{j}),events_averg.g2.(fn1{j}),events_averg.g4.(fn1{j}));
    end    
events_averg.all.be = events_averg.all.blockend; events_averg.all = rmfield(events_averg.all, "blockend");


%% averages around events (gives one value for a specific window around the event to compare to bsl)
fn1=fieldnames(events_averg.all);mel = {'_averg'};
    for j=1:numel(fn1)
        meanV =mean(events_averg.all.(fn1{j})(15:30,:),1); %15:30 
        events_averg.all.(char(strcat(fn1{j},mel))) = meanV;
    end 
% averages with different event windows
events_averg.all.soc_averg = mean(events_averg.all.soc(10:25,:),1);%10:25
events_averg.all.drink_averg = mean(events_averg.all.drink(30:45,:),1);%30:45
events_averg.all.eat_averg = mean(events_averg.all.eat(30:45,:),1);%30:45
events_averg.all.home_averg = mean(events_averg.all.be(30:90,:),1);%30:90
events_averg.all.expl_averg = mean(events_averg.all.expl(10:25,:),1);%10:25
events_averg.all.drink_antic_averg = mean(events_averg.all.drink(14:29,:),1);%
events_averg.all.eat_antic_averg = mean(events_averg.all.eat(14:29,:),1);%

%% binary
% std of the population to compare to 
std_ed = std(events_averg.all.ed(1:10,:), 1);
std_ef = std(events_averg.all.ef(1:10,:), 1);
std_es = std(events_averg.all.es(1:10,:), 1);
std_ee = std(events_averg.all.ee(1:10,:), 1);
std_er = std(events_averg.all.er(1:10,:), 1);
std_be = std(events_averg.all.be(1:10,:), 1);

%figure; plot(events_averg.all.ed_averg); hold on;plot(std_ed);
%figure; plot(events_averg.all.ef_averg); hold on;plot(std_ef);
%figure; plot(events_averg.all.drink_averg); hold on;plot(std_ed);
%figure; plot(events_averg.all.eat_averg); hold on;plot(std_ef);
%figure; plot(events_averg.all.es_averg); hold on;plot(std_es);
%figure; plot(events_averg.all.soc_averg); hold on;plot(std_es);
%figure; plot(std_ed);hold on;plot(std_ef);plot(std_es);plot(std_es);plot(std_ee);plot(std_er);plot(std_be);
%cll =30;% plots the average over all cells
%figure; plot(events_averg.all.ed(:,cll), LineWidth=2);hold on;plot(events_averg.all.ef(:,cll), LineWidth=2);plot(events_averg.all.es(:,cll), LineWidth=2);plot(events_averg.all.ee(:,cll), LineWidth=2);plot(events_averg.all.er(:,cll), LineWidth=2);plot(events_averg.all.be(1:31,cll), LineWidth=2);

fact= 5; % deciding factor, X times the standard deviation of the population in the baseline window
events=[];
column= 0;
fn3 = fieldnames(events_averg.all);
for hh=21:43
    column= column+1;
    for kk=1:59
        if ((hh == 21) || (hh ==27) || (hh ==29) || (hh ==36)|| (hh ==37)|| (hh ==42)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ed(kk))% all events compared to the std of entry drink
            events(column,kk) = 1;
            %kk = kk + 1;
        elseif  ((hh == 21) || (hh ==27) || (hh ==29) || (hh ==36)|| (hh ==37)|| (hh ==42)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ed(kk)))
               events(column,kk) = 2;
  
        elseif  ((hh == 22)||(hh == 28)||(hh == 30)|| (hh ==38)|| (hh ==39)|| (hh ==43)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ef(kk))
               events(column,kk) = 1;
        elseif  ((hh == 22)||(hh == 28)||(hh == 30)|| (hh ==38)|| (hh ==39)|| (hh ==43)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ef(kk)))
               events(column,kk) = 2;    

        elseif ( (hh == 23)||(hh == 24)||(hh == 31)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_es(kk))
               events(column,kk) = 1;
        elseif ( (hh == 23)||(hh == 24)||(hh == 31)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_es(kk)))
               events(column,kk) = 2; 


        elseif  ((hh == 25)||(hh == 32)||(hh == 35))& (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ee(kk))
               events(column,kk) = 1;
        elseif  ((hh == 25)||(hh == 32)||(hh == 35)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ee(kk)))
               events(column,kk) = 2; 

        elseif ( (hh == 26)||(hh == 33)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_er(kk))
               events(column,kk) = 1;
        elseif ( (hh == 26)||(hh == 33) ) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_er(kk)))
               events(column,kk) = 2; 

        elseif ( (hh == 40)||(hh == 41)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_be(kk))
               events(column,kk) = 1;
        elseif ( (hh == 40)||(hh == 41))  & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_be(kk)))
               events(column,kk) = 2;        

        else    
            events(column,kk) = 0;
        end    
    end
    
end

[r,c]=find(events==1);%finds ones excited cells
[r2,c2]=find(events ==2);%finds 2 inhibited cells

fn3 = fieldnames(events_averg.all);
fn4 = fn3(21:43,1);
for ev = 1:numel(fn4)
    events_averg.cells.(char(strrep((fn4{ev}), '_averg', '_cells'))) = c(find(r==ev));
    events_averg.cells.(char(strrep((fn4{ev}), '_averg', '_cells_inh'))) = c2(find(r2==ev));
end
% antic cells are only cells which are not also entry cells
keep1=ismember(events_averg.cells.drink_antic_cells, events_averg.cells.ed_cells);keep2=ismember(events_averg.cells.eat_antic_cells, events_averg.cells.ef_cells);
events_averg.cells.drink_antic_cells = events_averg.cells.drink_antic_cells(find(keep1 ==0));events_averg.cells.eat_antic_cells = events_averg.cells.eat_antic_cells(find(keep2 ==0));%without entry cells
keep1=ismember(events_averg.cells.drink_antic_cells_inh, events_averg.cells.ed_cells_inh);keep2=ismember(events_averg.cells.eat_antic_cells_inh, events_averg.cells.ef_cells_inh);
events_averg.cells.drink_antic_cells_inh = events_averg.cells.drink_antic_cells_inh(find(keep1 ==0));events_averg.cells.eat_antic_cells_inh = events_averg.cells.eat_antic_cells_inh(find(keep2 ==0));%without entry cells

% as vars in workspace
ed_cells =c(find(r==1));ef_cells =c(find(r==2));es_cells =c(find(r==3));ee_cells =c(find(r==5));er_cells =c(find(r==6));
drink_antic_cells =c(find(r==22));drink_cells =c(find(r==7));eat_antic_cells =c(find(r==23));eat_cells =c(find(r==8));
be_cells =c(find(r==20));home_cells =c(find(r==21));soc_cells=c(find(r==4));run_cells=c(find(r==14));expl_cells=c(find(r==15));
exd_cells =c(find(r==9));exf_cells =c(find(r==10));exs_cells =c(find(r==11));exe_cells =c(find(r==12));exr_cells =c(find(r==13));
ed_c_cells =c(find(r==17));ef_c_cells =c(find(r==19));ed_nc_cells =c(find(r==16));ef_nc_cells =c(find(r==18));
keep1=ismember(drink_antic_cells, ed_cells);keep2=ismember(eat_antic_cells, ef_cells);
drink_antic_cells = drink_antic_cells(find(keep1 ==0));eat_antic_cells = eat_antic_cells(find(keep2 ==0));%without entry cells

ed_cells_inh =c2(find(r2==1));ef_cells_inh =c2(find(r2==2));es_cells_inh =c2(find(r2==3));ee_cells_inh =c2(find(r2==5));er_cells_inh =c2(find(r2==6));
drink_antic_cells_inh =c2(find(r2==22));drink_cells_inh =c2(find(r2==7));eat_antic_cells_inh =c2(find(r2==23));eat_cells_inh =c2(find(r2==8));
be_cells_inh =c2(find(r2==20));home_cells_inh =c2(find(r2==21));soc_cells_inh=c2(find(r2==4));run_cells_inh=c2(find(r2==14));expl_cells_inh=c(find(r2==15));
exd_cells_inh =c2(find(r2==9));exf_cells_inh =c2(find(r2==10));exs_cells_inh =c2(find(r2==11));exe_cells_inh =c2(find(r2==12));exr_cells_inh =c2(find(r2==13));
ed_c_cells_inh =c(find(r2==17));ef_c_cells_inh =c(find(r2==19));ed_nc_cells_inh =c(find(r2==16));ef_nc_cells_inh =c(find(r2==18));
keep1=ismember(drink_antic_cells_inh, ed_cells_inh);keep2=ismember(eat_antic_cells_inh, ef_cells_inh);
drink_antic_cells_inh = drink_antic_cells_inh(find(keep1 ==0));eat_antic_cells_inh = eat_antic_cells_inh(find(keep2 ==0));%without entry cells

%% percentages
all=59;
fn1=fieldnames(events_averg.cells);
for j=1:numel(fn1)
    events_averg.perc.(fn1{j}) = [(size(events_averg.cells.(fn1{j}),1)/all)*100 ; size(events_averg.cells.(fn1{j}),1)];
end
%% 
noEventCells = find((sum(events,1)) == 0); size(noEventCells,2)
specializedCells = find((sum(events,1)) == 1);%cells which only have one function
%% Specialization score
specScore = [];
for evnt = 1:size(events,2)
    specNr = 0;
    for cl = 1:size(events,1)
        if events(cl,evnt) ~= 0 
        specNr = specNr + 1;
        else
        specNr = specNr; 
        end
    end
    specScore(1,evnt) = specNr;
end    

dd = unique(specScore);NrSpecCells = [0:(size(dd,2) -1)];
for cNr = 1:size(dd,2)
    NrSpecCells(2,cNr) = size(find(specScore == dd(cNr)),2);
end    
% plot specialization score
figure; bar(NrSpecCells(2,:));hold on;
xticks([1:size(dd,2)]);xticklabels(NrSpecCells(1,:));xlabel('\fontsize{12}Nr of encoded events');
ylabel('\fontsize{12}Nr of cells');hold off;title('\fontsize{11}Specialization score')

%% piechart event cells, specialized and no event cells
figure;label={' \fontsize{12}no event cells', '\fontsize{12}specialized cells',' \fontsize{12}event cells' };pieCells = [size(noEventCells,2), size(specializedCells,2), 59-(size(noEventCells,2) + size(specializedCells,2))];
pie(pieCells, '%.1f%%');lgd = legend(label);

%% compare nr of cells before entry vs consum

ExConsC = unique(cat(1,c(r == 4),c(r == 7),c(r == 8), c(r == 21)));
ExEntryC = unique(cat(1,c(r == 1),c(r == 2),c(r == 3), c(r == 5) ,c(r == 6) ,c (r == 20)));
isboth = ismember(ExConsC,ExEntryC);ExBoth = ExConsC(isboth == 1);
isExConsConly = ismember(ExConsC, ExBoth); ExConsConly = ExConsC(isExConsConly == 0); 
isExEntryConly = ismember(ExEntryC, ExBoth); ExEntryConly = ExEntryC(isExEntryConly == 0);

InhConsC = unique(cat(1,c2(r2 == 4),c2(r2 == 7),c2(r2 == 8), c2(r2 == 21)));
InhEntryC = unique(cat(1,c2(r2 == 1),c2(r2 == 2),c2(r2 == 3), c2(r2 == 5) ,c2(r2 == 6) ,c2 (r2 == 20)));
isboth = ismember(InhConsC,InhEntryC);InhBoth = InhConsC(isboth == 1);
isInhConsConly = ismember(InhConsC, InhBoth); InhConsConly = InhConsC(isInhConsConly == 0); 
isInhEntryConly = ismember(InhEntryC, InhBoth); InhEntryConly = InhEntryC(isInhEntryConly == 0);

ystacked = [size(InhEntryConly,1) size(InhBoth,1) size(InhConsConly,1);size(ExEntryConly,1) size(ExBoth,1) size(ExConsConly,1)];
figure; bar(ystacked,'stacked'); hold on; xticklabels({'Inh', 'Exc'});clear ylabel;ylabel('\fontsize{12}Nr of cells');
label={' \fontsize{12}Entry', '\fontsize{12}Both',' \fontsize{12}Consumption'};lgd = legend(label);

%% all trials from one cell
fnEA = fieldnames(events_averg);
anim = 3;
fnEAanimal = fieldnames(events_averg.(fnEA{anim}));fnEAanimal = fnEAanimal(17:end,1);% baselined individual traces
Event = 2; CellID =7;
%
Cell1 = events_averg.(fnEA{anim}).(fnEAanimal{Event})(:,CellID,:);
C1 = reshape(Cell1 ,[],size(Cell1,3));
%entry
figure; plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');%title(['\fontsize{10}entry drink inhibited cell']);

%cons
%figure; plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
%xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1', '0', '1', '2', '3'});xlabel('\fontsize{12}time from event[s]');xlim([1 61]);
%ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');%title(['\fontsize{10}entry drink inhibited cell']);

%% averg traces from multiple cells
fnAll = fieldnames(events_averg.all);
Event = 1;
Cells = [14 16 5 57 6 36];

figure; plot(events_averg.all.(fnAll{Event}) (:, Cells), LineWidth=2);hold on; xlim([1 30]); 
xticks([1 5 10 15 20 25 30]); xticklabels({'-3','-2.5' ,'-2','-1.5' ,'-1' , '-0.5', '0'});xlabel('\fontsize{12}time from event [s]'); ylabel('\fontsize{12}averg Ca-activity (z-scored, s.d.)');
plot([1 10], [-0.9 -0.9],"k", LineWidth=3);plot([15 30], [-0.9 -0.9],"k", LineWidth=3);
text(2.5,-0.8,'baseline');text(16,-0.8,'event window');
legend({'\fontsize{10}c14(enter drink ex)','\fontsize{10}c16(enter drink ex)', '\fontsize{10}c5  (enter drink inh)','\fontsize{10}c43 (enter food ex)','\fontsize{10}c3 (enter food inh)','\fontsize{10}c36 (no event)'}, Location="northwest");

%% raw traces of cells





%% imgsc plots of individual cells
cellNr = 7; tit = ['\fontsize{14}Food anticipatory cell'];
cell = events_averg.g4.eat_averg(:,cellNr,:);
cNr = reshape(cell,[],size(cell,3));
cNr = cNr(:,25:100);

%lim = [-1 9];
fA = figure;
figure(fA);hold on;
        s1= subplot(2,1,1);imagesc([0:3],[1:size(cNr,2)],cNr');xlim([0 3]);
        title(tit);xticklabels({'','','','','','',''});ylabel('\fontsize{12}Trials');
        c = colorbar(s1,'Position',[0.93 0.583837209302326 0.022 0.341162790697675]);
        s2= subplot(2,1,2);stdshade(cNr',0.3,[0 0 1]);
        xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]')
        ylabel('\fontsize{12}Ca-activity (z-score)');
        set(s2, 'Position', [0.13,0.38,0.775,0.20]);


%% imgsc plots of averages
% entries
% trials 
fntrials = fieldnames(events_averg);Entries = [1 2 3 5 6 9];clear NrTrials;
for an = 1:3
    fntrials2 = fieldnames(events_averg.(fntrials{an}));
    for entr = 1:size(Entries,2)
        NrTrials(entr,an) = size(events_averg.(fntrials{an}).(fntrials2{Entries(entr)}),3);
    end
end
NrTrials = array2table(NrTrials,'VariableNames',{'g5','g2','g4'},'RowNames',{'ed','ef','es','ee', 'er', 'be'});


win= 30;
lim = [-2 2];
limaverg = [-0.3 0.5];
f7=figure;% averages across 2all sessions
figure(f7);hold on;
subplot(2,6,1),imagesc([-win:0],[1:size(events_averg.all.ed,2)],events_averg.all.ed');xline(-20,'--b');
        title(['\fontsize{10}enter drink',' trials: 62-88']);clim(lim);
        xticklabels({' ',' ',' ',''});ylabel('\fontsize{15}Cells');
        sp1= subplot(2,6,7); stdshade(events_averg.all.ed',0.3,[0 0 1]);set(sp1, 'Position', [0.13,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-scored)');
        ylim(limaverg);
subplot(2,6,2),imagesc([-win:0],[1:size(events_averg.all.ef,2)],events_averg.all.ef');xline(-20,'--b');
        title(['\fontsize{10}enter feed',' trials: 97-160']);clim(lim)
        xticklabels({' ',' ',' ',''});
        sp2= subplot(2,6,8), stdshade(events_averg.all.ef',0.3,[0 0 1]);set(sp2, 'Position', [0.2645,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
subplot(2,6,3),imagesc([-win:0],[1:size(events_averg.all.es,2)],events_averg.all.es');xline(-20,'--b');
        title(['\fontsize{10}enter social',' trials: 75-148']);clim(lim);
         xticklabels({' ',' ',' ',''});xlabel('\fontsize{15}time [s]');
        sp3= subplot(2,6,9); stdshade(events_averg.all.es',0.3,[0 0 1]);set(sp3, 'Position', [0.3991,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
subplot(2,6,4),imagesc([-win:0],[1:size(events_averg.all.ee,2)],events_averg.all.ee');xline(-20,'--b');
        title(['\fontsize{10}enter explore',' trials: 20-33']);clim(lim);
         xticklabels({' ',' ',' ',''});
        sp4= subplot(2,6,10); stdshade(events_averg.all.ee',0.3,[0 0 1]);set(sp4, 'Position', [0.5336,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
subplot(2,6,5),imagesc([-win:0],[1:size(events_averg.all.er,2)],events_averg.all.er');xline(-20,'--b');
        title(['\fontsize{10}enter run',' trials: 8-25']);clim(lim);
        xticklabels({'-3','-2','-1',' '});
        sp5= subplot(2,6,11); stdshade(events_averg.all.er',0.3,[0 0 1]);set(sp5, 'Position', [0.6682,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
subplot(2,6,6),imagesc([-win:0],[1:size(events_averg.all.be(1:31,:),2)],events_averg.all.be(1:31,:)');xline(-20,'--b');
        title(['\fontsize{10}block end',' trials: 10-71']);clim(lim);
        xticks([-30 -20 -10 0]); xticklabels({' ',' ',' ',' '});
        sp6= subplot(2,6,12); stdshade((events_averg.all.be(1:31,:))',0.3,[0 0 1]);set(sp6, 'Position', [0.8027,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
        h = axes(f7,'visible','off');
        c = colorbar(h,'Position',[0.93 0.425 0.022 0.5]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 3.2;
        caxis(h,lim);


%% consumption
fntrials = fieldnames(events_averg);Entries = [7 8 4 16 15 9];clear NrTrialsCons;
for an = 1:3
    fntrials2 = fieldnames(events_averg.(fntrials{an}));
    for entr = 1:size(Entries,2)
        NrTrialsCons(entr,an) = size(events_averg.(fntrials{an}).(fntrials2{Entries(entr)}),3);
    end
end
NrTrialsCons = array2table(NrTrialsCons,'VariableNames',{'g5','g2','g4'},'RowNames',{'drink','eat','social','explore', 'run ', 'home'});


lim = [-3 3];limaverg2= [-1 1];
f9=figure;% averages across all sessions
figure(f9);hold on;
subplot(2,4,1),imagesc([-win:win],[1:size(events_averg.all.drink,2)],events_averg.all.drink');xline(0,'--w');
        title(['\fontsize{10}drink',' trials: 46-70']);clim(lim)
        xticks([1 11 21 31 41 51 61]); xticklabels({'-3','-2','-1','0','1','2','3'});;ylabel('\fontsize{15}Cells');xlabel('\fontsize{15}time [s]');
        sp7= subplot(2,4,5), stdshade(events_averg.all.drink',0.3,[0 0 1]);set(sp7, 'Position', [0.13,0.43,0.1566,0.15]);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');
        ylim(limaverg2);
subplot(2,4,2),imagesc([-win:win],[1:size(events_averg.all.eat,2)],events_averg.all.eat');xline(0,'--w');
        title(['\fontsize{10}retrieve pellet',' trials: 35-97']);clim(lim)
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});
        sp8= subplot(2,4,6), stdshade(events_averg.all.eat',0.3,[0 0 1]);set(sp8, 'Position', [0.3361,0.43,0.1566,0.15]);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});
        ylim(limaverg2);
subplot(2,4,3),imagesc([0:win],[1:size(events_averg.all.soc,2)],events_averg.all.soc');
        title(['\fontsize{10}social interaction',' trials: 75-148']);clim(lim)
        xticks([1 11 21 31]); xticklabels({' ','1','2',' '});
        sp9= subplot(2,4,7), stdshade(events_averg.all.soc',0.3,[0 0 1]);set(sp9, 'Position', [0.5422,0.43,0.1566,0.15]);
        xticks([1 11 21 31]);xticklabels({'0','1','2','3'}); 
        ylim(limaverg2);
subplot(2,4,4),imagesc([0:100],[1:size(events_averg.all.be(31:131,:),2)],events_averg.all.be(31:131,:)');xline(0,'--w');
        title(['\fontsize{10}home',' trials: 10-71']);clim(lim)
        xticks([1 21 41 61 81 101]);xticklabels({' ','2','4','6','8',' '});
        sp10= subplot(2,4,8), stdshade((events_averg.all.be(31:131,:))',0.3,[0 0 1]);set(sp10, 'Position', [0.7484,0.43,0.1566,0.15]);
        xticks([1 21 41 61 81 101]);xticklabels({'0','2','4','6','8','10'});
        ylim(limaverg2);
        h = axes(f9,'visible','off');
        c = colorbar(h,'Position',[0.93 0.45 0.02 0.45]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 2.5;
        caxis(h,lim);        


%% entries consumption vs non consumption trials
fntrials = fieldnames(events_averg);Entries = [33:36];clear NrTrialsNonCons;
for an = 1:3
    fntrials2 = fieldnames(events_averg.(fntrials{an}));
    for entr = 1:size(Entries,2)
        NrTrialsNonCons(entr,an) = size(events_averg.(fntrials{an}).(fntrials2{Entries(entr)}),3);
    end
end
NrTrialsNonCons = array2table(NrTrialsNonCons,'VariableNames',{'g5','g2','g4'},'RowNames',{'ed_nc','ed_c','ef_nc','ef_c'});

lim = [-3 3];limaverg2= [-1 1];
f10=figure;
figure(f10);hold on;
subplot(2,4,1),imagesc([-win:0],[1:size(events_averg.all.ed_nc,2)],events_averg.all.ed_nc');
        title(['\fontsize{10}enter no drink',' trials: 42']);clim(lim);
        xticks([1 11 21 31 ]); xticklabels({'-3','-2','-1','0'});;ylabel('\fontsize{15}Cells');xlabel('\fontsize{15}time [s]');
        sp7= subplot(2,4,5); stdshade(events_averg.all.ed_nc',0.3,[0 0 1]);set(sp7, 'Position', [0.13,0.43,0.1566,0.15]);
        xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');
        ylim(limaverg2);
subplot(2,4,2),imagesc([-win:0],[1:size(events_averg.all.ed_c,2)],events_averg.all.ed_c');
        title(['\fontsize{10}enter drink',' trials: 46-82']);clim(lim);
        xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});
        sp8= subplot(2,4,6), stdshade(events_averg.all.ed_c',0.3,[0 0 1]);set(sp8, 'Position', [0.3361,0.43,0.1566,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg2);
subplot(2,4,3),imagesc([-win:0],[1:size(events_averg.all.ef_nc,2)],events_averg.all.ef_nc');
        title(['\fontsize{10}enter no feed',' trials: 9-73']);clim(lim);
        xticks([1 11 21 31]); xticklabels({' ','1','2',' '});
        sp9= subplot(2,4,7); stdshade(events_averg.all.ef_nc',0.3,[0 0 1]);set(sp9, 'Position', [0.5422,0.43,0.1566,0.15]);
        xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg2);
subplot(2,4,4),imagesc([-win:0],[1:size(events_averg.all.ef_c,2)],events_averg.all.ef_c');
        title(['\fontsize{10}enter feed',' trials: 34-97']);clim(lim);
        xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});
        sp10= subplot(2,4,8); stdshade(events_averg.all.ef_c',0.3,[0 0 1]);set(sp10, 'Position', [0.7484,0.43,0.1566,0.15]);
       xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg2);
        h = axes(f10,'visible','off');
        c = colorbar(h,'Position',[0.93 0.45 0.02 0.45]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 2.5;
        caxis(h,lim);        

%% ethogram plots
s=5;
event_type=table2array(ci_data.bsl.g5.d_2022_09_08(s).session.events(:,4));
raw_f_res =(ci_data.bsl.g5.d_2022_09_08(s).session.raw_f);
%df_f_trace = dFoF(raw_f_res);%df/f
%z_scored_raw= zscore(df_f_trace,0,1);
z_scored_raw= zscore(raw_f_res,0,1);

x = (table2array(ci_data.bsl.g5.d_2022_09_08(s).session.events(:,5)));
y=NaN(size(x));
frameplot=raw_f;


%y(find(event_type=='block_start'))=2;
y(find(event_type=='imaging_start'))=2.5; 
y(find(event_type=='imaging_stop'))=2;
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
y(find(event_type=='block_end'))=2.5;
y(find(event_type=='imaging_stop'))=2;
y(find(event_type=='exit_feed'))=3.5;
y(find(event_type=='exit_drink'))=3.5;
y(find(event_type=='exit_run'))=3.5;
y(find(event_type=='exit_social'))=3.5;
y(find(event_type=='exit_explore'))=3.5;


yy=fillmissing(y,'linear');

figure;a=1;
subplot(1,1,a),plot([x(1) x(end)],[3.5 3.5],'b');hold on%decisionpoint
subplot(1,1,a),plot([x(1) x(end)],[2 2],'b');%nest

subplot(1,1,a),plot(x,yy,'k');hold on
subplot(1,1,a),plot(x,y,'ko');hold on


%title(animal);
yticks([0.5 1.5 2 2.5 3 3.5 4 4.5 5]);
yticklabels({'frames','frames','nest','drink','run','decision point','food','social','explore'});
ylabel('area');
ylim([-40 5.5]);

for i=1:1:length(pel)
    subplot(1,1,a),plot([pel(i) pel(i)],[3.8 4.2],'y','LineWidth',3);
end
for i=1:1:length(drink)
    subplot(1,1,a),plot([drink(i) drink(i)],[2.4 2.8],'c','LineWidth',3);
end
for i=1:1:length(run)
    subplot(1,1,a),plot([run(i) run(i)],[2.8 3.2],'m','LineWidth',3);
end

good_cells = [1:23];
j=0;
for i=good_cells
   subplot(1,1,a),plot((z_scored_raw(:,i)/5-j),'Color',[i/size(raw_f,2) 1-i/size(raw_f,2) i/size(raw_f,2)],'LineWidth',1);hold on;
   j=j+1;
end


