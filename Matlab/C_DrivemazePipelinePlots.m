clearvars -except events_averg events_averg_shuf ci_data perc95 perc5 Perc95 Perc5; close all;
set(0,'defaultAxesFontSize',14);

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
events_averg.all.be = events_averg.all.blockend; 


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

events_averg.all = rmfield(events_averg.all, "blockend_averg");
events_averg.all = rmfield(events_averg.all, "run_averg");

events_averg.all = rmfield(events_averg.all, "ed_nc_averg");
events_averg.all = rmfield(events_averg.all, "ed_c_averg");
events_averg.all = rmfield(events_averg.all, "ef_nc_averg");
events_averg.all = rmfield(events_averg.all, "ef_c_averg");

%% binary
% std of the population to compare to 
% std_ed = std(events_averg.all.ed(1:10,:), 1);
% std_ef = std(events_averg.all.ef(1:10,:), 1);
% std_es = std(events_averg.all.es(1:10,:), 1);
% std_ee = std(events_averg.all.ee(1:10,:), 1);
% std_er = std(events_averg.all.er(1:10,:), 1);
% std_be = std(events_averg.all.be(1:10,:), 1);

%fact= 5; % deciding factor, X times the standard deviation of the population in the baseline window
% modify perc files so they match data 

%remove blockend, arrays dont match otherwise
perc95=perc95([1:8 10:19],:); Perc95=Perc95([1:8 10:19],:);
perc5=perc5([1:8 10:19],:); Perc5=Perc5([1:8 10:19],:);


% compare to shuffled thresholds
events=[];
column= 0;
fn3 = fieldnames(events_averg.all);start= 21;
for hh=22:39
    column= column+1;
    for kk=1:59
        if events_averg.all.(fn3{hh})(1,kk) > perc95(hh-start,kk)
            events(column,kk) = 1;
        elseif events_averg.all.(fn3{hh})(1,kk) < perc5(hh-start,kk)  
            events(column,kk) = 2;
        else    
            events(column,kk) = 0;
        end
    end
end



[r,c]=find(events==1);%finds ones excited cells
[r2,c2]=find(events ==2);%finds 2 inhibited cells

fn3 = fieldnames(events_averg.all);
fn4 = fn3(22:39,1);
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
drink_antic_cells =c(find(r==17));drink_cells =c(find(r==7));eat_antic_cells =c(find(r==18));eat_cells =c(find(r==8));
be_cells =c(find(r==15));home_cells =c(find(r==16));soc_cells=c(find(r==4));expl_cells=c(find(r==14));%run_cells=c(find(r==14))
exd_cells =c(find(r==9));exf_cells =c(find(r==10));exs_cells =c(find(r==11));exe_cells =c(find(r==12));exr_cells =c(find(r==13));
%ed_c_cells =c(find(r==16));ef_c_cells =c(find(r==18));ed_nc_cells =c(find(r==15));ef_nc_cells =c(find(r==17));
keep1=ismember(drink_antic_cells, ed_cells);keep2=ismember(eat_antic_cells, ef_cells);
drink_antic_cells = drink_antic_cells(find(keep1 ==0));eat_antic_cells = eat_antic_cells(find(keep2 ==0));%without entry cells

ed_cells_inh =c2(find(r2==1));ef_cells_inh =c2(find(r2==2));es_cells_inh =c2(find(r2==3));ee_cells_inh =c2(find(r2==5));er_cells_inh =c2(find(r2==6));
drink_antic_cells_inh =c2(find(r2==17));drink_cells_inh =c2(find(r2==7));eat_antic_cells_inh =c2(find(r2==18));eat_cells_inh =c2(find(r2==8));
be_cells_inh =c2(find(r2==15));home_cells_inh =c2(find(r2==16));soc_cells_inh=c2(find(r2==4));expl_cells_inh=c2(find(r2==14));%run_cells_inh=c2(find(r2==14))
exd_cells_inh =c2(find(r2==9));exf_cells_inh =c2(find(r2==10));exs_cells_inh =c2(find(r2==11));exe_cells_inh =c2(find(r2==12));exr_cells_inh =c2(find(r2==13));
%ed_c_cells_inh =c(find(r2==16));ef_c_cells_inh =c(find(r2==18));ed_nc_cells_inh =c(find(r2==15));ef_nc_cells_inh =c(find(r2==17));
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

dd = unique(specScore);NrSpecCells= dd; 
for cNr = 1:size(NrSpecCells,2)
    NrSpecCells(2,cNr) = size(find(specScore == dd(cNr)),2);
end    



% plot specialization score
figure; bar(NrSpecCells(2,:));hold on;
xticks([1:size(dd,2)]);xticklabels(NrSpecCells(1,:));xlabel('\fontsize{12}Nr of encoded events');
ylabel('\fontsize{12}Nr of cells');hold off;title('\fontsize{11}Specialization score')

%% spec score without exits?

%% piechart event cells, specialized and no event cells
figure;label={' \fontsize{12}no event cells',' \fontsize{12}event cells', '\fontsize{12}specialized cells' };pieCells = [size(noEventCells,2), 59-(size(noEventCells,2) + size(specializedCells,2)), size(specializedCells,2)];
pie(pieCells, '%.1f%%');lgd = legend(label);

%% venn diagrams
A = 1:59;
%drink module
setListData2={A ed_cells drink_antic_cells drink_cells};
setLabels2=[ " " ;"entry drink ";"drink antic";"drink"];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);
hh.ShowIntersectionCounts = true;

size(find(ismember(ed_cells, drink_antic_cells, "rows")),1)
size(find(ismember(ed_cells, drink_cells, "rows")),1)
size(find(ismember(drink_cells, drink_antic_cells, "rows")),1)


%eat module
size(find(ismember(ef_cells, eat_antic_cells, "rows")),1)
size(find(ismember(ef_cells, eat_cells, "rows")),1)
size(find(ismember(eat_cells, eat_antic_cells, "rows")),1)




%% bar plots 
% entries
Entries = [size(ed_cells,1) size(ed_cells_inh,1); size(ef_cells,1) size(ef_cells_inh,1); size(es_cells,1)  size(es_cells_inh,1);size(er_cells,1) size(er_cells_inh,1);size(ee_cells,1) size(ee_cells_inh,1);size(be_cells,1) size(be_cells_inh,1)];
X = categorical({'drink','food','social','run','explore','end'});X = reordercats(X,{'drink','food','social','run','explore', 'end'});
figure; bar(X,Entries,'stacked'); ylabel('\fontsize{12}Nr of cells');
% overlap
edef = size(find(ismember(ed_cells, ef_cells, "rows")),1);
edes = size(find(ismember(ed_cells, es_cells, "rows")),1);
eder = size(find(ismember(ed_cells, er_cells, "rows")),1);
edee = size(find(ismember(ed_cells, ee_cells, "rows")),1);
edbe = size(find(ismember(ed_cells, be_cells, "rows")),1);

EntriesD = [(size(ed_cells,1)-edef) edef (size(ef_cells,1)-edef); (size(ed_cells,1)-edes) edes (size(es_cells,1)-edes); (size(ed_cells,1)-eder) eder (size(er_cells,1)-eder);(size(ed_cells,1)-edee) edee (size(ee_cells,1)-edee);(size(ed_cells,1)-edbe) edbe (size(be_cells,1)-edbe)]
X = categorical({'drink vs food','drink vs social','drink vs run', 'drink vs expl', 'drink vs be'});X = reordercats(X,{'drink vs food','drink vs social','drink vs run', 'drink vs expl', 'drink vs be'});
figure; bar(X,EntriesD,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([0 30]);

% overlap greater than chance?
A= size(ed_cells,1)/59;
B= size(be_cells,1)/59;
ABfound = edbe/59
AnBexpected = (A*B)

% enter food
efed = size(find(ismember(ef_cells, ed_cells, "rows")),1);
efes = size(find(ismember(ef_cells, es_cells, "rows")),1);
efer = size(find(ismember(ef_cells, er_cells, "rows")),1);
efee = size(find(ismember(ef_cells, ee_cells, "rows")),1);
efbe = size(find(ismember(ef_cells, be_cells, "rows")),1);

EntriesF = [(size(ef_cells,1)-efed) efed (size(ed_cells,1)-efed);(size(ef_cells,1)-efes) efes (size(es_cells,1)-efes); (size(ef_cells,1)-efer) efer (size(er_cells,1)-efer); (size(ef_cells,1)-efee) efee (size(ee_cells,1)-efee); (size(ef_cells,1)-efbe) efbe (size(ee_cells,1)-efbe)]
X = categorical({'food vs drink','food vs social','food vs run','food vs explore','food vs be'});X = reordercats(X,{'food vs drink','food vs social','food vs run','food vs explore','food vs be'});
figure; bar(X,EntriesF,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([0 30]);

% overlap greater than chance?
A= size(ef_cells,1)/59;
B= size(ee_cells,1)/59;
ABfound = efee/59
AnBexpected = (A*B)

%enter social
esed = size(find(ismember(es_cells, ed_cells, "rows")),1);
esef = size(find(ismember(es_cells, ef_cells, "rows")),1);
esbe = size(find(ismember(es_cells, be_cells, "rows")),1);
eser = size(find(ismember(es_cells, er_cells, "rows")),1);
esee = size(find(ismember(es_cells, ee_cells, "rows")),1);

EntriesS = [ (size(es_cells,1)-esed) esed (size(ed_cells,1)-esed);(size(es_cells,1)-esef) esef (size(ef_cells,1)-esef);(size(es_cells,1)-eser) eser (size(er_cells,1)-eser); (size(es_cells,1)-esee) esee (size(ee_cells,1)-esee);(size(es_cells,1)-esbe) esbe (size(be_cells,1)-esbe)]
X = categorical({'soc vs drink','soc vs food','soc vs run','soc vs explore','soc vs be'});X = reordercats(X,{'soc vs drink','soc vs food','soc vs run','soc vs explore','soc vs be'});
figure; bar(X,EntriesS,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([0 30]);
% overlap greater than chance?
A= size(es_cells,1)/59;
B= size(be_cells,1)/59;
ABfound = esbe/59
AnBexpected = (A*B)

%% consumption %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cons = [size(drink_cells,1) size(drink_cells_inh,1); size(eat_cells,1) size(eat_cells_inh,1); size(soc_cells,1)  size(soc_cells_inh,1);size(expl_cells,1) size(expl_cells_inh,1);size(home_cells,1) size(home_cells_inh,1)];
X = categorical({'drink','eat','social','explore','home'});X = reordercats(X,{'drink','eat','social','explore','home'});
figure; bar(X,Cons,'stacked'); ylabel('\fontsize{12}Nr of cells');
% overlap drink 
drinkEat = size(find(ismember(drink_cells, eat_cells, "rows")),1);
drinkSoc = size(find(ismember(drink_cells, soc_cells, "rows")),1);
drinkExpl = size(find(ismember(drink_cells, expl_cells, "rows")),1);
drinkHome = size(find(ismember(drink_cells, home_cells, "rows")),1);

ConsAntic = [(size(drink_cells,1)-drinkEat) drinkEat (size(eat_cells,1)-drinkEat); (size(drink_cells,1)-drinkSoc) drinkSoc (size(soc_cells,1)-drinkSoc);(size(drink_cells,1)-drinkExpl) drinkExpl (size(expl_cells,1)-drinkExpl);(size(drink_cells,1)-drinkHome) drinkHome (size(home_cells,1)-drinkHome)]
X = categorical({'drink vs eat','drink vs soc','drink vs expl','drink vs home'});X = reordercats(X,{'drink vs eat','drink vs soc','drink vs expl','drink vs home'});
figure; bar(X,ConsAntic,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([0 40]);


% overlap greater than chance?
A= size(drink_cells,1)/59;
B= size(home_cells,1)/59;
ABfound = drinkHome/59
AnBexpected = (A*B)

% Eat 
eatDrink = size(find(ismember(eat_cells, drink_cells, "rows")),1);
eatSoc = size(find(ismember(eat_cells, soc_cells, "rows")),1);
eatExpl = size(find(ismember(eat_cells, expl_cells, "rows")),1);
eatHome = size(find(ismember(eat_cells, home_cells, "rows")),1);

ConsAntic = [(size(eat_cells,1)-eatDrink) eatDrink (size(drink_cells,1)-eatDrink); (size(eat_cells,1)-eatSoc) eatSoc (size(soc_cells,1)-eatSoc);(size(eat_cells,1)-eatExpl) eatExpl (size(expl_cells,1)-eatExpl);(size(eat_cells,1)-eatHome) eatHome (size(home_cells,1)-eatHome)]
X = categorical({'eat vs eat','eat vs soc','eat vs expl','eat vs home'});X = reordercats(X,{'eat vs eat','eat vs soc','eat vs expl','eat vs home'});
figure; bar(X,ConsAntic,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([0 40]);


% overlap greater than chance?
A= size(eat_cells,1)/59;
B= size(home_cells,1)/59;
ABfound = eatHome/59
AnBexpected = (A*B)

% Social 
socDrink = size(find(ismember(soc_cells, drink_cells, "rows")),1);
socEat = size(find(ismember(soc_cells, eat_cells, "rows")),1);
socExpl = size(find(ismember(soc_cells, expl_cells, "rows")),1);
socHome = size(find(ismember(soc_cells, home_cells, "rows")),1);

ConsAntic = [(size(soc_cells,1)-socDrink) socDrink (size(drink_cells,1)-socDrink); (size(soc_cells,1)-socEat) socEat (size(eat_cells,1)-socEat);(size(soc_cells,1)-socExpl) socExpl (size(expl_cells,1)-socExpl);(size(soc_cells,1)-socHome) socHome (size(home_cells,1)-socHome)]
X = categorical({'soc vs soc','soc vs eat','soc vs expl','soc vs home'});X = reordercats(X,{'soc vs soc','soc vs eat','soc vs expl','soc vs home'});
figure; bar(X,ConsAntic,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([0 40]);


% overlap greater than chance?
A= size(soc_cells,1)/59;
B= size(home_cells,1)/59;
ABfound = socHome/59
AnBexpected = (A*B)

%Inhibitory 

drinkEat = size(find(ismember(drink_cells, eat_cells, "rows")),1);
drinkInhEatInh = size(find(ismember(drink_cells_inh, eat_cells_inh, "rows")),1);
drinkInhEat= size(find(ismember(drink_cells_inh, eat_cells, "rows")),1);
EatInhDrink= size(find(ismember(eat_cells_inh, drink_cells, "rows")),1);

ConsAntic = [(size(drink_cells,1)-drinkEat) drinkEat (size(eat_cells,1)-drinkEat); (size(drink_cells_inh,1)-drinkInhEatInh) drinkInhEatInh (size(eat_cells_inh,1)-drinkInhEatInh); (size(drink_cells_inh,1)-drinkInhEat) drinkInhEat (size(eat_cells,1)-drinkInhEat); (size(drink_cells,1)-EatInhDrink) EatInhDrink (size(eat_cells_inh,1)-EatInhDrink)]
X = categorical({'drink vs eat','drink inh vs eat inh','drink inh vs eat','drink vs eat inh'});X = reordercats(X,{'drink vs eat','drink inh vs eat inh','drink inh vs eat','drink vs eat inh'});
figure; bar(X,ConsAntic,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([0 40]);


% overlap greater than chance?
A= size(drink_cells,1)/59;
B= size(eat_cells_inh,1)/59;
ABfound = EatInhDrink/59
AnBexpected = (A*B)


% antic cells
drinkAnticDrink = size(find(ismember(drink_antic_cells, drink_cells, "rows")),1);
eatAnticEat = size(find(ismember(eat_antic_cells, eat_cells, "rows")),1);
drinkAnticEatAntic = size(find(ismember(eat_antic_cells, drink_antic_cells, "rows")),1);

ConsAntic = [(size(drink_antic_cells,1)-drinkAnticDrink) drinkAnticDrink (size(drink_cells,1)-drinkAnticDrink); (size(eat_antic_cells,1)-eatAnticEat) eatAnticEat (size(eat_cells,1)-eatAnticEat); (size(eat_antic_cells,1)-drinkAnticEatAntic) drinkAnticEatAntic (size(drink_antic_cells,1)-drinkAnticEatAntic)]; 
X = categorical({'drink antic vs drink','eat antic vs eat ', 'eat antic vs drink antic'});X = reordercats(X,{'drink antic vs drink','eat antic vs eat ','eat antic vs drink antic'});
figure; bar(X,ConsAntic,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([0 25]);


% overlap greater than chance?
A= size(eat_antic_cells,1)/59;
B= size(drink_antic_cells,1)/59;
ABfound = drinkAnticEatAntic/59
AnBexpected = (A*B)





% social
socialEat = size(find(ismember(soc_cells, eat_cells, "rows")),1);
socialDrink = size(find(ismember(soc_cells, drink_cells, "rows")),1);

ConsAntic = [(size(soc_cells,1)-socialEat) socialEat (size(eat_cells,1)-socialEat); (size(soc_cells,1)-socialDrink) socialDrink (size(drink_cells,1)-socialDrink)]; 
X = categorical({'social vs eat','social vs drink '});X = reordercats(X,{'social vs eat','social vs drink '});
figure; bar(X,ConsAntic,'stacked'); ylabel('\fontsize{12}Nr of cells');%ylim([0 25]);


exploreEat = size(find(ismember(expl_cells, eat_cells, "rows")),1);
exploreDrink = size(find(ismember(expl_cells, drink_cells, "rows")),1);
socialExplore = size(find(ismember(soc_cells, expl_cells, "rows")),1);

ConsAntic = [(size(soc_cells,1)-socialEat) socialEat (size(eat_cells,1)-socialEat); (size(soc_cells,1)-socialDrink) socialDrink (size(drink_cells,1)-socialDrink)]; 
X = categorical({'social vs eat','social vs drink '});X = reordercats(X,{'social vs eat','social vs drink '});
figure; bar(X,ConsAntic,'stacked'); ylabel('\fontsize{12}Nr of cells');%ylim([0 25]);

%% scatterplots
 f1= figure; 
 figure(f1);hold on,
    subplot(2,2,1);scatter(events_averg.all.eat_averg,events_averg.all.drink_averg,	"g", LineWidth=1);
    xlabel('\fontsize{12}Ca+-response eat');ylabel('\fontsize{12}drink ');hold on;plot([0 2],[0 0],"k");plot([0 0],[0 2],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-2, 0.5],	"m");
    xlim([-1 2]);  %xlim([-1 1])
    scatter(([events_averg.all.eat_averg(1,(eat_cells(:,1)))]),([events_averg.all.drink_averg(1,(eat_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.eat_averg(1,(eat_cells_inh(:,1)))]),([events_averg.all.drink_averg(1,(eat_cells_inh(:,1)))]),'b', LineWidth=1);
    scatter(([events_averg.all.eat_averg(1,(drink_cells(:,1)))]),([events_averg.all.drink_averg(1,(drink_cells(:,1)))]),"*",'r');scatter(([events_averg.all.eat_averg(1,(drink_cells_inh(:,1)))]),([events_averg.all.drink_averg(1,(drink_cells_inh(:,1)))]),"*",'b');
    
    subplot(2,2,2);scatter(events_averg.all.eat_antic_averg,events_averg.all.drink_antic_averg,	"g", LineWidth=1);hold on;plot([0 2],[0 0],"k");plot([0 0],[0 2],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response eat anticipatory');ylabel('\fontsize{12}drink anticipatory ');
    xlim([-1 2]);
    scatter(([events_averg.all.eat_antic_averg(1,(eat_antic_cells(:,1)))]),([events_averg.all.drink_antic_averg(1,(eat_antic_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.eat_antic_averg(1,(eat_antic_cells_inh(:,1)))]),([events_averg.all.drink_antic_averg(1,(eat_antic_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.eat_antic_averg(1,(drink_antic_cells(:,1)))]),([events_averg.all.drink_antic_averg(1,(drink_antic_cells(:,1)))]),"*",'r');scatter(([events_averg.all.eat_antic_averg(1,(drink_antic_cells_inh(:,1)))]),([events_averg.all.drink_antic_averg(1,(drink_antic_cells_inh(:,1)))]),"*",'b');
    
    subplot(2,2,3);scatter(events_averg.all.eat_averg,events_averg.all.soc_averg,"g", LineWidth=1);hold on;plot([0 2],[0 0],"k");plot([0 0],[0 2],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-2, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response eat');ylabel('\fontsize{12}social contact');
    xlim([-1 2]);
    scatter(([events_averg.all.eat_averg(1,(eat_cells(:,1)))]),([events_averg.all.soc_averg(1,(eat_cells(:,1)))]),'r',  LineWidth=1);scatter(([events_averg.all.eat_averg(1,(eat_cells_inh(:,1)))]),([events_averg.all.soc_averg(1,(eat_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.eat_averg(1,(soc_cells(:,1)))]),([events_averg.all.soc_averg(1,(soc_cells(:,1)))]),"*",'r');scatter(([events_averg.all.eat_averg(1,(soc_cells_inh(:,1)))]),([events_averg.all.soc_averg(1,(soc_cells_inh(:,1)))]),"*",'b');
    
%     subplot(2,2,4);scatter(events_averg.all.be_averg,events_averg.all.home_averg,"g", LineWidth=1);hold on;plot([0 2],[0 0],"k");plot([0 0],[0 3],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
%     xlabel('\fontsize{12}Ca+-response block end');ylabel('\fontsize{12}homecage');
%     xlim([-1 2]);
%     scatter(([events_averg.all.be_averg(1,(be_cells(:,1)))]),([events_averg.all.home_averg(1,(be_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.be_averg(1,(be_cells_inh(:,1)))]),([events_averg.all.home_averg(1,(be_cells_inh(:,1)))]),'b',LineWidth=1);
%     scatter(([events_averg.all.be_averg(1,(home_cells(:,1)))]),([events_averg.all.home_averg(1,(home_cells(:,1)))]),"*",'r');scatter(([events_averg.all.be_averg(1,(home_cells_inh(:,1)))]),([events_averg.all.home_averg(1,(home_cells_inh(:,1)))]),"*",'b');

    subplot(2,2,4);scatter(events_averg.all.be_averg,events_averg.all.home_averg,"g", LineWidth=1);hold on;plot([0 2],[0 0],"k");plot([0 0],[0 3],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response block end');ylabel('\fontsize{12}homecage');
    xlim([-1 2]);
    scatter(([events_averg.all.be_averg(1,(be_cells(:,1)))]),([events_averg.all.home_averg(1,(be_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.be_averg(1,(be_cells_inh(:,1)))]),([events_averg.all.home_averg(1,(be_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.be_averg(1,(home_cells(:,1)))]),([events_averg.all.home_averg(1,(home_cells(:,1)))]),"*",'r');scatter(([events_averg.all.be_averg(1,(home_cells_inh(:,1)))]),([events_averg.all.home_averg(1,(home_cells_inh(:,1)))]),"*",'b');

f2= figure; lims = [-1 1];
 figure(f2);hold on;
    subplot(2,2,1);scatter(events_averg.all.ed_averg,events_averg.all.ef_averg,	"g", LineWidth=1);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");%plot(X,Y)
    xlabel('\fontsize{12}Ca+-response enter drink');ylabel('\fontsize{12}enter food ');
    xlim([-1 1.5]); 
    scatter(([events_averg.all.ed_averg(1,(ed_cells(:,1)))]),([events_averg.all.ef_averg(1,(ed_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.ed_averg(1,(ed_cells_inh(:,1)))]),([events_averg.all.ef_averg(1,(ed_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.ed_averg(1,(ef_cells(:,1)))]),([events_averg.all.ef_averg(1,(ef_cells(:,1)))]),"*",'r');scatter(([events_averg.all.ed_averg(1,(ef_cells_inh(:,1)))]),([events_averg.all.ef_averg(1,(ef_cells_inh(:,1)))]),"*",'b');
    legend({'ed cells ex','ed cells inh','ef cells ex' , 'ef cells inh'});
    subplot(2,2,2);scatter(events_averg.all.ed_averg,events_averg.all.es_averg,	"g", LineWidth=1);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter drink');ylabel('\fontsize{12}enter social ');
    xlim([-1 1.5]);
    scatter(([events_averg.all.ed_averg(1,(ed_cells(:,1)))]),([events_averg.all.es_averg(1,(ed_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.ed_averg(1,(ed_cells_inh(:,1)))]),([events_averg.all.er_averg(1,(ed_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.ed_averg(1,(es_cells(:,1)))]),([events_averg.all.es_averg(1,(es_cells(:,1)))]),"*",'r');scatter(([events_averg.all.ed_averg(1,(es_cells_inh(:,1)))]),([events_averg.all.es_averg(1,(es_cells_inh(:,1)))]),"*",'b');


    subplot(2,2,3);scatter(events_averg.all.ed_averg,events_averg.all.er_averg,	"g", LineWidth=1);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter drink');ylabel('\fontsize{12}enter run');
    xlim([-1 1.5]);
    scatter(([events_averg.all.ed_averg(1,(ed_cells(:,1)))]),([events_averg.all.er_averg(1,(ed_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.ed_averg(1,(ed_cells_inh(:,1)))]),([events_averg.all.er_averg(1,(ed_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.ed_averg(1,(er_cells(:,1)))]),([events_averg.all.er_averg(1,(er_cells(:,1)))]),"*",'r');scatter(([events_averg.all.ed_averg(1,(er_cells_inh(:,1)))]),([events_averg.all.er_averg(1,(er_cells_inh(:,1)))]),"*",'b');


    subplot(2,2,4);scatter(events_averg.all.ed_averg,events_averg.all.be_averg,	"g", LineWidth=1);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter drink');ylabel('\fontsize{12}block end');
    xlim([-1 1.5]);   
    scatter(([events_averg.all.ed_averg(1,(ed_cells(:,1)))]),([events_averg.all.be_averg(1,(ed_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.ed_averg(1,(ed_cells_inh(:,1)))]),([events_averg.all.ef_averg(1,(ed_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.ed_averg(1,(be_cells(:,1)))]),([events_averg.all.be_averg(1,(be_cells(:,1)))]),"*",'r');scatter(([events_averg.all.be_averg(1,(be_cells_inh(:,1)))]),([events_averg.all.be_averg(1,(be_cells_inh(:,1)))]),"*",'b');



%% compare nr of cells before entry vs consum

ExConsC = unique(cat(1,c(r == 4),c(r == 7),c(r == 8), c(r == 16)));
ExEntryC = unique(cat(1,c(r == 1),c(r == 2),c(r == 3), c(r == 5) ,c(r == 6) ,c (r == 15)));
isboth = ismember(ExConsC,ExEntryC);ExBoth = ExConsC(isboth == 1);
isExConsConly = ismember(ExConsC, ExBoth); ExConsConly = ExConsC(isExConsConly == 0); 
isExEntryConly = ismember(ExEntryC, ExBoth); ExEntryConly = ExEntryC(isExEntryConly == 0);

InhConsC = unique(cat(1,c2(r2 == 4),c2(r2 == 7),c2(r2 == 8), c2(r2 == 16)));
InhEntryC = unique(cat(1,c2(r2 == 1),c2(r2 == 2),c2(r2 == 3), c2(r2 == 5) ,c2(r2 == 6) ,c2 (r2 == 15)));
isboth = ismember(InhConsC,InhEntryC);InhBoth = InhConsC(isboth == 1);
isInhConsConly = ismember(InhConsC, InhBoth); InhConsConly = InhConsC(isInhConsConly == 0); 
isInhEntryConly = ismember(InhEntryC, InhBoth); InhEntryConly = InhEntryC(isInhEntryConly == 0);

ystacked = [size(InhEntryConly,1) size(InhBoth,1) (size(InhConsConly,1) + size(drink_antic_cells_inh,1) + size(eat_antic_cells_inh,1)) ;size(ExEntryConly,1) size(ExBoth,1) (size(ExConsConly,1)+ size(drink_antic_cells,1) + size(eat_antic_cells,1))];
figure; bar(ystacked,'stacked'); hold on; xticklabels({'Inh', 'Exc'});clear ylabel;ylabel('\fontsize{12}Nr of cells');
label={' \fontsize{12}Appetitive', '\fontsize{12} Both',' \fontsize{12}Consummatory'};lgd = legend(label);

% overlap greater than chance?
A= size(InhEntryC,1)/59;
B= (size(InhConsC,1) + size(drink_antic_cells_inh,1) + size(eat_antic_cells_inh,1))/59;
ABobserved = size(InhBoth,1)/59
AnBexpected = (A*B)

ABobserved*59
AnBexpected*59

A= size(ExEntryC,1)/59;
B= (size(ExConsC,1)+ size(drink_antic_cells,1) + size(eat_antic_cells,1))/59;
ABobserved = size(ExBoth,1)/59
AnBexpected = (A*B)



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
Cells = [14 16 5 57 6 37];

figure; plot(events_averg.all.(fnAll{Event}) (:, Cells), LineWidth=2);hold on; xlim([1 30]); 
xticks([1 5 10 15 20 25 30]); xticklabels({'-3','-2.5' ,'-2','-1.5' ,'-1' , '-0.5', '0'});xlabel('\fontsize{12}time from event [s]'); ylabel('\fontsize{12}averg Ca-activity (z-scored, s.d.)');
%plot([1 10], [-0.9 -0.9],"k", LineWidth=3);plot([15 30], [-0.9 -0.9],"k", LineWidth=3);
%text(2.5,-0.8,'baseline');text(16,-0.8,'event window');
legend({'\fontsize{10}c14','\fontsize{10}c16', '\fontsize{10}c5','\fontsize{10}c57','\fontsize{10}c6','\fontsize{10}c37'}, Location="northwest");

% for shuffles



%% raw traces of cells
enter = 2337;%ef 2361 pellet 2303 leave social
c1 = ci_data.bsl.g5.d_2022_09_08(17).session.raw_f_res(enter-200:enter+100,4);z_c1 = (c1-nanmean(c1(1:50,1)))/nanstd(c1(1:50,1));
c2 = ci_data.bsl.g5.d_2022_09_08(17).session.raw_f_res(enter-200:enter+100,8);z_c2 = (c2-nanmean(c2(1:50,1)))/nanstd(c2(1:50,1));
c3 = ci_data.bsl.g5.d_2022_09_08(17).session.raw_f_res(enter-200:enter+100,6);z_c3 = (c3-nanmean(c3(1:50,1)))/nanstd(c3(1:50,1));
c4 = ci_data.bsl.g5.d_2022_09_08(17).session.raw_f_res(enter-200:enter+100,2);z_c4 = (c4-nanmean(c4(1:50,1)))/nanstd(c4(1:50,1));
%x-nanmean(b)/nanstd(b)
exmpl_ed = table(smooth(zscore(c1)), smooth(zscore(c2)), smooth(zscore(c3)), smooth(zscore(c4)));
%figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
%plot(smooth(zscore(c40_ed)));hold on;plot(smooth(zscore(c52_ed)));plot(smooth(zscore(c56_ed)));plot(smooth(zscore(c54_ed)));
p1=plot(smooth(z_c1));hold on;p2=plot(smooth(z_c2));p3=plot(smooth(z_c3));%p4=plot(smooth(z_c54));
p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;%p4.LineWidth=1;
xl= xline(200,'-','entry food');xl.LineWidth = 1; 
xxl = xline(200+24,'-','cons food');xxl.LineWidth = 1;
xxxl = xline(200-34,'-','exit social');xxxl.LineWidth = 1;
xticks([0 50 100 150 200 250 300]);xticklabels({'-20','-15','-10','-5','0','5','10'});xlabel('\fontsize{12}time [s]');xlim([0 300]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');xlim([0 300]);
legend({'\fontsize{12}c4 (eat ex)','\fontsize{12}c8 (eat ex)', '\fontsize{12}c6 (ef/eat inh)'}, Location="northwest");


enter = 1606;
c1 = ci_data.bsl.g4.d_2022_09_02(4).session.raw_f_res(enter-100:enter+100,4);z_c1 = (c1-nanmean(c1(1:50,1)))/nanstd(c1(1:50,1));
c2 = ci_data.bsl.g4.d_2022_09_02(4).session.raw_f_res(enter-100:enter+100,14);z_c2 = (c2-nanmean(c2(1:50,1)))/nanstd(c2(1:50,1));
c3 = ci_data.bsl.g4.d_2022_09_02(4).session.raw_f_res(enter-100:enter+100,16);z_c3 = (c3-nanmean(c3(1:50,1)))/nanstd(c3(1:50,1));
c4 = ci_data.bsl.g4.d_2022_09_02(4).session.raw_f_res(enter-100:enter+100,20);z_c4 = (c4-nanmean(c4(1:50,1)))/nanstd(c4(1:50,1));
%x-nanmean(b)/nanstd(b)
exmpl_ed = table(smooth(zscore(c1)), smooth(zscore(c2)), smooth(zscore(c3)), smooth(zscore(c4)));
%figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
%plot(smooth(zscore(c40_ed)));hold on;plot(smooth(zscore(c52_ed)));plot(smooth(zscore(c56_ed)));plot(smooth(zscore(c54_ed)));
p1=plot(smooth(z_c1));hold on;p2=plot(smooth(z_c2));p3=plot(smooth(z_c3));p4=plot(smooth(z_c4));
p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;p4.LineWidth=1;
xl= xline(199,'-','entry drink');xl.LineWidth = 1; 
xxl= xline(199+20,'-','drink');xxl.LineWidth = 1;
xxxl= xline(200+263,'-','exit drink');xxxxl.LineWidth = 1;
xticks([0 50 100 150 200 250 300]);xticklabels({'-20','-15','-10','-5','0','5','10'});xlabel('\fontsize{12}time [s]')
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');xlim([0 300]);
legend({'\fontsize{12}c40','\fontsize{12}c50', '\fontsize{12}c52', '\fontsize{12}c56'}, Location="northwest");



enter = 6847;s=5;pluswin =500;
c1 = ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res(enter-100:enter+pluswin,4);z_c1 = (c1-nanmean(c1(1:50,1)))/nanstd(c1(1:50,1));
c2 = ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res(enter-100:enter+pluswin,14);z_c2 = (c2-nanmean(c2(1:50,1)))/nanstd(c2(1:50,1));
c3 = ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res(enter-100:enter+pluswin,16);z_c3 = (c3-nanmean(c3(1:50,1)))/nanstd(c3(1:50,1));
c4 = ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res(enter-100:enter+pluswin,20);z_c4 = (c4-nanmean(c4(1:50,1)))/nanstd(c4(1:50,1));
%x-nanmean(b)/nanstd(b)
exmpl_ed = table(smooth(zscore(c1)), smooth(zscore(c2)), smooth(zscore(c3)), smooth(zscore(c4)));
%figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
%plot(smooth(zscore(c40_ed)));hold on;plot(smooth(zscore(c52_ed)));plot(smooth(zscore(c56_ed)));plot(smooth(zscore(c54_ed)));
p1=plot(smooth(z_c1));hold on;p2=plot(smooth(z_c2));p3=plot(smooth(z_c3));p4=plot(smooth(z_c4));
p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;p4.LineWidth=1;
xl= xline(100,'-','entry drink');xl.LineWidth = 1; 
xxl= xline(100+35,'-','drink');xxl.LineWidth = 1;
xxxl= xline(100+165,'-','exit drink');xxxl.LineWidth = 1;
xxxxl= xline(100+186,'-','entry food');xxxxl.LineWidth = 1;
xxxxxl= xline(100+390,'-','eat');xxxxxl.LineWidth = 1;
xticks([0 100 200 300 400 500 600]);xticklabels({'0','10','20','30','40','50','60'});xlabel('\fontsize{12}time [s]')
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');xlim([0 600]);
legend({'\fontsize{12}c4','\fontsize{12}c14', '\fontsize{12}c16', '\fontsize{12}c20'}, Location="northwest");



enter = 6847;s=5;pluswin =50;minuswin = 50;
c1 = ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res(enter-minuswin:enter+pluswin,4);z_c1 = (c1-nanmean(c1(1:50,1)))/nanstd(c1(1:50,1));
c2 = ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res(enter-minuswin:enter+pluswin,14);z_c2 = (c2-nanmean(c2(1:50,1)))/nanstd(c2(1:50,1));
c3 = ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res(enter-minuswin:enter+pluswin,16);z_c3 = (c3-nanmean(c3(1:50,1)))/nanstd(c3(1:50,1));
c4 = ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res(enter-minuswin:enter+pluswin,20);z_c4 = (c4-nanmean(c4(1:50,1)))/nanstd(c4(1:50,1));
%x-nanmean(b)/nanstd(b)
exmpl_ed = table(smooth(zscore(c1)), smooth(zscore(c2)), smooth(zscore(c3)), smooth(zscore(c4)));
%figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
%plot(smooth(zscore(c40_ed)));hold on;plot(smooth(zscore(c52_ed)));plot(smooth(zscore(c56_ed)));plot(smooth(zscore(c54_ed)));
p1=plot(smooth(z_c1));hold on;p2=plot(smooth(z_c2));p3=plot(smooth(z_c3));p4=plot(smooth(z_c4));
p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;p4.LineWidth=1;
xl= xline(50,'-','entry drink');xl.LineWidth = 1; 
xxl= xline(50+35,'-','drink');xxl.LineWidth = 1;
%xxxl= xline(100+165,'-','exit drink');xxxl.LineWidth = 1;
%xxxxl= xline(100+186,'-','entry food');xxxxl.LineWidth = 1;
%xxxxxl= xline(100+390,'-','eat');xxxxxl.LineWidth = 1;
xticks([0 10 20 30 40 50 60 70 80 90 100 110 120 130 140 150]);xticklabels({'-5','-4','-3','-2','-1','0','1','2','3','4','5'});xlabel('\fontsize{12}time [s]')
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');xlim([0 100]);
%legend({'\fontsize{12}c4','\fontsize{12}c14', '\fontsize{12}c16', '\fontsize{12}c20'}, Location="northwest");



%% imgsc plots of individual cells
cellNr =7; tit = ['\fontsize{14}social cell'];%30 42
cell = events_averg.g4.soc_averg_bsl(:,cellNr,:); 
cNr = reshape(cell,[],size(cell,3));
cNr = cNr(:,1:end);

%lim = [-1 9];3
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
f7=figure;% averages across all sessions
figure(f7);hold on;
subplot(2,6,1),imagesc([-win:0],[1:size(events_averg.all.ed,2)],events_averg.all.ed');xline(-20,'--b');
        title(['\fontsize{10}enter drink',' trials: 62-88']);clim(lim);
        xticklabels({' ',' ',' ',''});ylabel('\fontsize{15}Cells');
        sp1= subplot(2,6,7); stdshade(events_averg.all.ed',0.3,[0 0 1]);set(sp1, 'Position', [0.13,0.43,0.102,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-scored)');
        ylim(limaverg);
subplot(2,6,2),imagesc([-win:0],[1:size(events_averg.all.ef,2)],events_averg.all.ef');xline(-20,'--b');
        title(['\fontsize{10}enter feed',' trials: 97-160']);clim(lim)
        xticklabels({' ',' ',' ',''});
        sp2= subplot(2,6,8), stdshade(events_averg.all.ef',0.3,[0 0 1]);set(sp2, 'Position', [0.2645,0.43,0.102,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
subplot(2,6,3),imagesc([-win:0],[1:size(events_averg.all.es,2)],events_averg.all.es');xline(-20,'--b');
        title(['\fontsize{10}enter social',' trials: 75-148']);clim(lim);
         xticklabels({' ',' ',' ',''});xlabel('\fontsize{15}time [s]');
        sp3= subplot(2,6,9); stdshade(events_averg.all.es',0.3,[0 0 1]);set(sp3, 'Position', [0.3991,0.43,0.102,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
subplot(2,6,4),imagesc([-win:0],[1:size(events_averg.all.ee,2)],events_averg.all.ee');xline(-20,'--b');
        title(['\fontsize{10}enter explore',' trials: 20-33']);clim(lim);
         xticklabels({' ',' ',' ',''});
        sp4= subplot(2,6,10); stdshade(events_averg.all.ee',0.3,[0 0 1]);set(sp4, 'Position', [0.5336,0.43,0.102,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
subplot(2,6,5),imagesc([-win:0],[1:size(events_averg.all.er,2)],events_averg.all.er');xline(-20,'--b');
        title(['\fontsize{10}enter run',' trials: 8-25']);clim(lim);
        xticklabels({'-3','-2','-1',' '});
        sp5= subplot(2,6,11); stdshade(events_averg.all.er',0.3,[0 0 1]);set(sp5, 'Position', [0.6682,0.43,0.102,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
subplot(2,6,6),imagesc([-win:0],[1:size(events_averg.all.be(1:31,:),2)],events_averg.all.be(1:31,:)');xline(-20,'--b');
        title(['\fontsize{10}block end',' trials: 10-71']);clim(lim);
        xticks([-30 -20 -10 0]); xticklabels({' ',' ',' ',' '});
        sp6= subplot(2,6,12); stdshade((events_averg.all.be(1:31,:))',0.3,[0 0 1]);set(sp6, 'Position', [0.8027,0.43,0.102,0.15]);
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

lim = [-3 3];limaverg2 = [-0.3 0.5];
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

%% ethogram plots
s=7;
event_type=table2array(ci_data.bsl.g5.d_2022_09_27(s).session.events(:,4));
raw_f_res =(ci_data.bsl.g5.d_2022_09_27(s).session.raw_f_res);
df_f_trace = dFoF(raw_f_res);%df/f
z_scored_raw= zscore(df_f_trace,0,1);
%z_scored_raw= zscore(raw_f_res,0,1);

x = (table2array(ci_data.bsl.g5.d_2022_09_27(s).session.events(:,5)));
y=NaN(size(x));
frameplot=raw_f_res;


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
   subplot(1,1,a),plot((z_scored_raw(:,i)/7-j),'Color',[i/size(raw_f_res,2) 1-i/size(raw_f_res,2) i/size(raw_f_res,2)],'LineWidth',1);hold on;
   j=j+1;
end












