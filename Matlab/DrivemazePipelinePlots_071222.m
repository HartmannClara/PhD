clearvars -except events_averg ci_data; close all;
set(0,'defaultAxesFontSize',18);

%%
fn=fieldnames(events_averg);me = {'_m'};
%loop through the fields to make mean all bsl trials
for ii=1: numel(fn)
    fn1=fieldnames(events_averg.(fn{ii}));
        for j=17:numel(fn1)
        events_averg.(fn{ii}).(char(strcat(fn1{j},me))) = nanmean(events_averg.(fn{ii}).(fn1{j}),3);
        end
end    

%% concatenating over averg from all cells
fn1=fieldnames(events_averg.g5);
    for j=33:numel(fn1)
        events_averg.all.(strrep((fn1{j}), '_averg_bsl_m', '')) = cat(2,events_averg.g5.(fn1{j}),events_averg.g2.(fn1{j}),events_averg.g4.(fn1{j}));
    end    
events_averg.all.be = events_averg.all.blockend; events_averg.all = rmfield(events_averg.all, "blockend");
%% averages around events (gives one value for a specific window around the event to compare to bsl)
events_averg.all.ed_averg = mean(events_averg.all.ed(15:30,:),1);
events_averg.all.ef_averg = mean(events_averg.all.ef(15:30,:),1);
events_averg.all.es_averg = mean(events_averg.all.es(15:30,:),1);
events_averg.all.soc_averg = mean(events_averg.all.soc(15:30,:),1);
events_averg.all.ee_averg = mean(events_averg.all.ee(15:30,:),1);
events_averg.all.er_averg = mean(events_averg.all.er(15:30,:),1);
events_averg.all.drink_averg = mean(events_averg.all.drink(30:45,:),1);
events_averg.all.drink_antic_averg = mean(events_averg.all.drink(15:30,:),1);
events_averg.all.eat_averg = mean(events_averg.all.eat(30:45,:),1);
events_averg.all.eat_antic_averg = mean(events_averg.all.eat(15:30,:),1);
events_averg.all.be_averg = mean(events_averg.all.be(15:30,:),1);
events_averg.all.home_averg = mean(events_averg.all.be(30:90,:),1);
events_averg.all.exd_averg = mean(events_averg.all.exd(15:30,:),1);
events_averg.all.exf_averg = mean(events_averg.all.exf(15:30,:),1);
events_averg.all.exs_averg = mean(events_averg.all.exs(15:30,:),1);
events_averg.all.exe_averg = mean(events_averg.all.exe(15:30,:),1);
events_averg.all.exr_averg = mean(events_averg.all.exr(15:30,:),1);
events_averg.all.run_averg = mean(events_averg.all.run(15:30,:),1);
events_averg.all.expl_averg = mean(events_averg.all.expl(15:30,:),1);

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

fact= 1.5; % deciding factor, 1.5 times the standard deviationof the population in the baseline window
events=[];
column= 0;
fn3 = fieldnames(events_averg.all);
for hh=17:35
    column= column+1;
    for kk=1:59
        if ((hh == 17) || (hh ==23) || (hh ==24) || (hh ==29)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ed)
            events(column,kk) = 1;
            %kk = kk + 1;
        elseif  ((hh == 17) || (hh ==23) || (hh ==24) || (hh ==29)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ed))
               events(column,kk) = 2;
  
        elseif  ((hh == 18)||(hh == 25)||(hh == 26)||(hh == 30)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ef)
               events(column,kk) = 1;
        elseif  ((hh == 18)||(hh == 25)||(hh == 26)||(hh == 30)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ef))
               events(column,kk) = 2;    

        elseif ( (hh == 19)||(hh == 20)||(hh == 31)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_es)
               events(column,kk) = 1;
        elseif  ((hh == 19)||(hh == 20)||(hh == 31)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_es))
               events(column,kk) = 2; 


        elseif  ((hh == 21)||(hh == 32))& (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ee)
               events(column,kk) = 1;
        elseif  ((hh == 21)||(hh == 32)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ee))
               events(column,kk) = 2; 

        elseif ( (hh == 22)||(hh == 33)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_er)
               events(column,kk) = 1;
        elseif ( (hh == 22)||(hh == 33) ) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_er))
               events(column,kk) = 2; 

        elseif ( (hh == 27)||(hh == 28)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_be)
               events(column,kk) = 1;
        elseif ( (hh == 27)||(hh == 28))  & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_be))
               events(column,kk) = 2;        

        else    
            events(column,kk) = 0;
        end    
    end
    
end

[r,c]=find(events==1);%finds ones
[r2,c2]=find(events ==2);%finds 2

fn3 = fieldnames(events_averg.all);
fn4 = fn3(17:35,1);
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
drink_antic_cells =c(find(r==8));drink_cells =c(find(r==7));eat_antic_cells =c(find(r==10));eat_cells =c(find(r==9));
be_cells =c(find(r==11));home_cells =c(find(r==12));soc_cells=c(find(r==4));run_cells=c(find(r==18));expl_cells=c(find(r==19));
exd_cells =c(find(r==13));exf_cells =c(find(r==14));exs_cells =c(find(r==15));exe_cells =c(find(r==16));exr_cells =c(find(r==17));
keep1=ismember(drink_antic_cells, ed_cells);keep2=ismember(eat_antic_cells, ef_cells);
drink_antic_cells = drink_antic_cells(find(keep1 ==0));eat_antic_cells = eat_antic_cells(find(keep2 ==0));%without entry cells

ed_cells_inh =c2(find(r2==1));ef_cells_inh =c2(find(r2==2));es_cells_inh =c2(find(r2==3));ee_cells_inh =c2(find(r2==5));er_cells_inh =c2(find(r2==6));
drink_antic_cells_inh =c2(find(r2==8));drink_cells_inh =c2(find(r2==7));eat_antic_cells_inh =c2(find(r2==10));eat_cells_inh =c2(find(r2==9));
be_cells_inh =c2(find(r2==11));home_cells_inh =c2(find(r2==12));soc_cells_inh=c2(find(r2==4));run_cells_inh=c2(find(r2==18));expl_cells_inh=c(find(r==19));
exd_cells_inh =c2(find(r2==13));exf_cells_inh =c2(find(r2==14));exs_cells_inh =c2(find(r2==15));exe_cells_inh =c2(find(r2==16));exr_cells_inh =c2(find(r2==17));
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
specializedCells = find((sum(events,1)) == 1);
%%compare nr of cells before entry vs consum
entryCellsEx = unique(cat(1, ed_cells, ef_cells, es_cells, ee_cells, er_cells, be_cells,ed_cells_inh, ef_cells_inh, es_cells_inh, ee_cells_inh, er_cells_inh, be_cells_inh));
%entryCellsInh = unique(cat(1, ed_cells_inh, ef_cells_inh, es_cells_inh, ee_cells_inh, er_cells_inh, be_cells_inh));
consCellsEx = unique(cat(1, drink_cells, eat_cells, soc_cells, expl_cells, run_cells, home_cells, drink_cells_inh, eat_cells_inh, soc_cells_inh, run_cells_inh, home_cells_inh));
%consCellsInh = unique(cat(1, drink_cells_inh, eat_cells_inh, soc_cells_inh, run_cells_inh, home_cells_inh));
y = [size(entryCellsEx,1), size(consCellsEx,1)];
bar(y);xticklabels({'entry cells','consumption cells'});
