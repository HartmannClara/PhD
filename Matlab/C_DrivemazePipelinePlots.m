clearvars -except events_averg ci_data; close all;
set(0,'defaultAxesFontSize',18);
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
        meanV =mean(events_averg.all.(fn1{j})(15:30,:),1);
        events_averg.all.(char(strcat(fn1{j},mel))) = meanV;
    end 
% averages with different event windows
events_averg.all.soc_averg = mean(events_averg.all.soc(10:25,:),1);
events_averg.all.drink_averg = mean(events_averg.all.drink(30:45,:),1);
events_averg.all.eat_averg = mean(events_averg.all.eat(30:45,:),1);
events_averg.all.home_averg = mean(events_averg.all.be(30:90,:),1);
events_averg.all.expl_averg = mean(events_averg.all.expl(10:25,:),1);
events_averg.all.drink_antic_averg = mean(events_averg.all.drink(14:29,:),1);
events_averg.all.eat_antic_averg = mean(events_averg.all.eat(14:29,:),1);

%% binary
% std of the population to compare to 
std_ed = std(events_averg.all.ed(1:10,:), 1);
std_ef = std(events_averg.all.ef(1:10,:), 1);
std_es = std(events_averg.all.es(1:10,:), 1);
std_ee = std(events_averg.all.ee(1:10,:), 1);
std_er = std(events_averg.all.er(1:10,:), 1);
std_be = std(events_averg.all.be(1:10,:), 1);
%figure; plot(std_ed);hold on;plot(std_ef);plot(std_es);plot(std_es);plot(std_ee);plot(std_er);plot(std_be);
%cll =30;% plots the average over all cells
%figure; plot(events_averg.all.ed(:,cll), LineWidth=2);hold on;plot(events_averg.all.ef(:,cll), LineWidth=2);plot(events_averg.all.es(:,cll), LineWidth=2);plot(events_averg.all.ee(:,cll), LineWidth=2);plot(events_averg.all.er(:,cll), LineWidth=2);plot(events_averg.all.be(1:31,cll), LineWidth=2);

fact= 3; % deciding factor, X times the standard deviationof the population in the baseline window
events=[];
column= 0;
fn3 = fieldnames(events_averg.all);
for hh=21:43
    column= column+1;
    for kk=1:59
        if ((hh == 21) || (hh ==27) || (hh ==29) || (hh ==36)|| (hh ==37)|| (hh ==42)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ed)% all events compared to the std of entry drink
            events(column,kk) = 1;
            %kk = kk + 1;
        elseif  ((hh == 21) || (hh ==27) || (hh ==29) || (hh ==36)|| (hh ==37)|| (hh ==42)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ed))
               events(column,kk) = 2;
  
        elseif  ((hh == 22)||(hh == 28)||(hh == 30)|| (hh ==38)|| (hh ==39)|| (hh ==43)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ef)
               events(column,kk) = 1;
        elseif  ((hh == 22)||(hh == 28)||(hh == 30)|| (hh ==38)|| (hh ==39)|| (hh ==43)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ef))
               events(column,kk) = 2;    

        elseif ( (hh == 23)||(hh == 24)||(hh == 31)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_es)
               events(column,kk) = 1;
        elseif ( (hh == 23)||(hh == 24)||(hh == 31)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_es))
               events(column,kk) = 2; 


        elseif  ((hh == 25)||(hh == 32)||(hh == 35))& (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ee)
               events(column,kk) = 1;
        elseif  ((hh == 25)||(hh == 32)||(hh == 35)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ee))
               events(column,kk) = 2; 

        elseif ( (hh == 26)||(hh == 33)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_er)
               events(column,kk) = 1;
        elseif ( (hh == 26)||(hh == 33) ) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_er))
               events(column,kk) = 2; 

        elseif ( (hh == 40)||(hh == 41)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_be)
               events(column,kk) = 1;
        elseif ( (hh == 40)||(hh == 41))  & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_be))
               events(column,kk) = 2;        

        else    
            events(column,kk) = 0;
        end    
    end
    
end

[r,c]=find(events==1);%finds ones
[r2,c2]=find(events ==2);%finds 2

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
noEventCells = find((sum(events,1)) == 0);
specializedCells = find((sum(events,1)) == 1);%cells which only have one function
%% piechart event cells, specialized and rest
figure;label={' no event cells', 'specialized cells',' event cells' };pieCells = [size(noEventCells,2), size(specializedCells,2), 59-(size(noEventCells,2) + size(specializedCells,2))];
pie(pieCells, '%.1f%%');lgd = legend(label);
%%compare nr of cells before entry vs consum
entryCellsEx = unique(cat(1, ed_cells, ef_cells, es_cells,  er_cells, be_cells));%,ed_cells_inh, ef_cells_inh, es_cells_inh,  er_cells_inh, be_cells_inh));
entryCellsInh = unique(cat(1, ed_cells_inh, ef_cells_inh, es_cells_inh, er_cells_inh, be_cells_inh));
consCellsEx = unique(cat(1, drink_cells, eat_cells, soc_cells,  run_cells, home_cells));%, drink_cells_inh, eat_cells_inh, soc_cells_inh, run_cells_inh, home_cells_inh));
consCellsInh = unique(cat(1, drink_cells_inh, eat_cells_inh, soc_cells_inh, run_cells_inh, home_cells_inh));
y = [size(entryCellsInh,1) size(entryCellsEx,1) ;size(consCellsInh,1) size(consCellsEx,1) ];
bar(y, 'stacked');xticklabels({'\fontsize{10}entry cells','\fontsize{10}consumption cells'}); lbl = {'inhibited', 'excited'}; lgd = legend(lbl);
%% 

win= 30;
lim = [-2 2];
limaverg = [-0.3 0.5];
f7=figure;% averages across 2all sessions
figure(f7);hold on;
        subplot(2,3,1),imagesc([-win:0],[1:size(events_averg.all.ed,2)],events_averg.all.ef');xline(-20,'--b');
        title(['\fontsize{10}enter drink',' trials: 97-160']);clim(lim)
        xticklabels({' ',' ',' ',''});ylabel('\fontsize{15}Cells');
        subplot(2,3,2),imagesc([-win:0],[1:size(events_averg.all.ef_nc,2)],events_averg.all.ef_nc');xline(-20,'--b');
        title(['\fontsize{10}enter drink nc',' trials: 9-73']);clim(lim)
        xticklabels({' ',' ',' ',''});ylabel('\fontsize{15}Cells');
        subplot(2,3,3),imagesc([-win:0],[1:size(events_averg.all.ef_c,2)],events_averg.all.ef_c');xline(-20,'--b');
        title(['\fontsize{10}enter drink c',' trials: 34-62']);clim(lim)
        xticklabels({' ',' ',' ',''});ylabel('\fontsize{15}Cells');

