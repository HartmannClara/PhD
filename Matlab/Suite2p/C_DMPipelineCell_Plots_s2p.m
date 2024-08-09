clearvars -except events_averg events_averg_data events_averg_ch ci_data perc95 perc5 Perc95 Perc5 Specs; close all;
set(0,'defaultAxesFontSize',14);
%% Plotting pipeline 
%load("events_averg_s2p.mat");
fn=fieldnames(events_averg);
% make new classification list for exits according to next entries and time
fnS = fieldnames(Specs);
for ii = 1:4
    exf_entry =[];exr_entry =[];exe_entry =[];exs_entry =[];exd_entry =[];
    exf_time = [];exr_time = [];exe_time = [];exs_time = [];exd_time = [];
    for j = 1: size(Specs.(fnS{ii}).DecTimes,1)
        exits = Specs.(fnS{ii}).DecTimes{j,2}(:,1);
        for jj = 1:numel(exits)
            nextEntry = []; TimeinD = [];
            exit = exits(jj);
            if  exit == "exit_feed"
                nextEntry = Specs.(fnS{ii}).DecTimes{j,2}(jj,2);
                TimeinD = Specs.(fnS{ii}).DecTimes{j,1}(jj,1);
                exf_entry =[exf_entry nextEntry];
                exf_time =[exf_time TimeinD];
                
            elseif  exit == "exit_run" 
                nextEntry = Specs.(fnS{ii}).DecTimes{j,2}(jj,2);
                TimeinD = Specs.(fnS{ii}).DecTimes{j,1}(jj,1);
                exr_entry =[exr_entry nextEntry];
                exr_time =[exr_time TimeinD];

            elseif  exit == "exit_explore"
                nextEntry = Specs.(fnS{ii}).DecTimes{j,2}(jj,2);
                TimeinD = Specs.(fnS{ii}).DecTimes{j,1}(jj,1);
                exe_entry =[exe_entry nextEntry];
                exe_time =[exe_time TimeinD];

            elseif  exit == "exit_social"
                nextEntry = Specs.(fnS{ii}).DecTimes{j,2}(jj,2);
                TimeinD = Specs.(fnS{ii}).DecTimes{j,1}(jj,1);
                exs_entry =[exs_entry nextEntry];
                exs_time =[exs_time TimeinD];

            elseif  exit == "exit_drink"
                nextEntry = Specs.(fnS{ii}).DecTimes{j,2}(jj,2);
                TimeinD = Specs.(fnS{ii}).DecTimes{j,1}(jj,1);
                exd_entry =[exd_entry nextEntry];
                exd_time =[exd_time TimeinD];

            end
        end
    end
    Specs.(fnS{ii}).ExitsSort.exf_entry =exf_entry';Specs.(fnS{ii}).ExitsSort.exf_time =exf_time';
    Specs.(fnS{ii}).ExitsSort.exr_entry =exr_entry';Specs.(fnS{ii}).ExitsSort.exr_time =exr_time';
    Specs.(fnS{ii}).ExitsSort.exe_entry =exe_entry';Specs.(fnS{ii}).ExitsSort.exe_time =exe_time';
    Specs.(fnS{ii}).ExitsSort.exs_entry =exs_entry';Specs.(fnS{ii}).ExitsSort.exs_time =exs_time';
    Specs.(fnS{ii}).ExitsSort.exd_entry =exd_entry';Specs.(fnS{ii}).ExitsSort.exd_time =exd_time;
end    

% now I am only interested in the exits which are long enough to not be
% included in the next entry window -> Deciscion time > 6s (3s enty window
% and 3s exit window, no overlap)
% fnE=fieldnames(events_averg.g5);
% remF = [1:27 34:38];
% EventExits = events_averg.g5; EventExits = rmfield(EventExits,fnE(remF));
% exfSnips=[];exrSnips=[];exeSnips=[];exsSnips=[];exdSnips=[];beSnips=[];
% for ii = 1:size(EventExits.exf_averg_bsl,3)
%     snip = EventExits.exf_averg_bsl(:,:,ii);
%     if Specs.g5.ExitsSort.exf_entry(ii,1) == "enter_food"
%         exfSnips = [exfSnips snip];
%     elseif Specs.g5.ExitsSort.exf_entry(ii,1) == "enter_run"
%         exrSnips = [exrSnips snip];
%     elseif Specs.g5.ExitsSort.exf_entry(ii,1) == "enter_explore"
%         exeSnips = [exeSnips snip];
%     elseif Specs.g5.ExitsSort.exf_entry(ii,1) == "enter_social"
%         exsSnips = [exsSnips snip];
%     elseif Specs.g5.ExitsSort.exf_entry(ii,1) == "enter_drink"
%         exdSnips = [exdSnips snip];
%     elseif Specs.g5.ExitsSort.exf_entry(ii,1) == "imaging_stop"
%         beSnips = [beSnips snip];    
%     end
% end    
% 


%% separating out the neuropil trace
% % separate average neuropil trace out 
% fn=fieldnames(events_averg);me = {'_NP'};
% %loop through the fields to make mean all bsl trials
% for ii=1: numel(fn)
%     fn1=fieldnames(events_averg.(fn{ii}));
%         for j=20:numel(fn1)
%             EV = events_averg.(fn{ii}).(fn1{j});
%             if isempty(EV) ~= true
%                events_averg.(fn{ii}).(char(strcat(fn1{j},me)))= EV(:,end,:); %collect traces of last cell == mean neuropil signal of all cells
%             else    
%                events_averg.(fn{ii}).(char(strcat(fn1{j},me))) = [];
%             end
%         end
% end   
% %remove neuropil from cells
% for ii=1: numel(fn)
%     fn1=fieldnames(events_averg.(fn{ii}));
%         for j=1:38
%             EV = events_averg.(fn{ii}).(fn1{j});
%             if isempty(EV) ~= true
%                events_averg.(fn{ii}).(fn1{j})= EV(:,1:end-1,:); %all cells but last one
%             else    
%                j=j+1;
%             end
%         end
% end   
% save("events_NPsort_s2p_.mat", "events_averg");
%% plotting neuropil traces during eating

% fnEA = fieldnames(events_averg);
% anim = 4;
% fnEAan = fieldnames(events_averg.(fnEA{anim}));
% fnEAanimal = fnEAan(20:38,1);% baselined individual traces
% fnEAanimalNP = fnEAan(39:end,1);
% Event = 8; CellID =17;
% %
% Cell1 = events_averg.(fnEA{anim}).(fnEAanimal{Event})(:,CellID,:);
% C1 = reshape(Cell1 ,[],size(Cell1,3));
% %NP
% Cell1NP = events_averg.(fnEA{anim}).(fnEAanimalNP{Event});
% C1NP = reshape(Cell1NP ,[],size(Cell1NP,3));
% %entry
% %figure; plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
% %xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);
% %ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');%title(['\fontsize{10}entry drink inhibited cell']);
% 
% %cons
% figure; plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
% xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1', '0', '1', '2', '3'});xlabel('\fontsize{12}time from event[s]');xlim([1 61]);
% ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');title(['\fontsize{10} Cell']);
% 
% %figure; plot(C1NP, LineWidth=1);hold on; plot((mean(C1NP,2)),'k' ,LineWidth=4);
% %xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1', '0', '1', '2', '3'});xlabel('\fontsize{12}time from event[s]');xlim([1 61]);
% %ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');title(['\fontsize{10} Neuropil']);

%% get cell numbers for each animal
fn=fieldnames(events_averg);
for ii=1:4
    fn1=fieldnames(events_averg.(fn{ii}));
        for j=1:19
        Cellnumbers(j, ii) = size(events_averg.(fn{ii}).(fn1{j}),2);
        end
end
Rownames = {(fn1(1:19))};
Cellnumbers = array2table(Cellnumbers,'VariableNames',{'g5','g4','g2','g12'},'RowNames',Rownames{1,1});% How many trials there are of every event typer per animal
%% get trial numbers for each animal
for ii=1: 4
    fn1=fieldnames(events_averg.(fn{ii}));
        for j=1:19
        Trialnumbers(j, ii) = size(events_averg.(fn{ii}).(fn1{j}),3);
        end
end
Trialnumbers(Trialnumbers==1) = 0;
Rownames = {(fn1(1:19))};
Trialnumbers = array2table(Trialnumbers,'VariableNames',{'g5','g4','g2','g12'},'RowNames',Rownames{1,1});% How many trials there are of every event typer per animal


%% makes an average over all trials per cell
me = {'_m'};
%loop through the fields to make mean all bsl trials
for ii=1: numel(fn)
    fn1=fieldnames(events_averg.(fn{ii}));
        for j=20:numel(fn1)
        events_averg.(fn{ii}).(char(strcat(fn1{j},me))) = nanmean(events_averg.(fn{ii}).(fn1{j}),3);
        end
end  
%% concatenating over averg from all cells (all cells together)
fn1=fieldnames(events_averg.g5);
    for j=39:54
        events_averg.all.(strrep((fn1{j}), '_averg_bsl_m', '')) = cat(2,events_averg.g5.(fn1{j}),events_averg.g4.(fn1{j}),events_averg.g2.(fn1{j}),events_averg.g12.(fn1{j}));
    end
    for j=55:numel(fn1)
        events_averg.all.(strrep((fn1{j}), '_averg_nc_bsl_m', '_nc')) = cat(2,events_averg.g5.(fn1{j}),events_averg.g4.(fn1{j}),events_averg.g2.(fn1{j}),events_averg.g12.(fn1{j}));
    end
% rename blockend to be!!!!!!!!!!!!!! 
events_averg.all.be = events_averg.all.blockend;
fna= fieldnames(events_averg.all);
newAllorder = [1:8 20 10:19];
for i = 1:numel(newAllorder)
    NewAll.(fna{newAllorder(i)}) = events_averg.all.(fna{newAllorder(i)});
end    
events_averg.all = NewAll;

%% averages around events (gives one value for a specific window around the event to compare to bsl)
fn1=fieldnames(events_averg.all);mel = {'_averg'};
    for j=1:numel(fn1)
        meanV =mean(events_averg.all.(fn1{j})(16:30,:),1); %16:30 
        events_averg.all.(char(strcat(fn1{j},mel))) = meanV;
    end 
% averages with different event windows
events_averg.all.soc_averg = mean(events_averg.all.soc(20:35,:),1);%20:35
events_averg.all.drink_averg = mean(events_averg.all.drink(31:45,:),1);%31:45
events_averg.all.eat_averg = mean(events_averg.all.eat(31:45,:),1);%30:45
events_averg.all.run_averg = mean(events_averg.all.run(31:45,:),1);%31:45
events_averg.all.home_averg = mean(events_averg.all.be(30:89,:),1);%30:89
events_averg.all.expl_averg = mean(events_averg.all.expl(20:35,:),1);%20:35
events_averg.all.drink_antic_averg = mean(events_averg.all.drink(16:30,:),1);%16:30
events_averg.all.eat_antic_averg = mean(events_averg.all.eat(16:30,:),1);%16:30

%events_averg.all = rmfield(events_averg.all, "run_averg");

% add missing cells to events with missing trials
% ze35=zeros(1,35);
% ze72=zeros(1,72);
% ze24=zeros(1,24);
% ze17=zeros(1,17);
% events_averg.all.er_averg = cat(2,events_averg.all.er_averg(:,1:74), ze35, events_averg.all.er_averg(:,75:122),ze72, ze24,ze17);
% events_averg.all.run_averg = cat(2,events_averg.all.run_averg(:,1:74), ze35, events_averg.all.run_averg(:,75:122),ze72, ze24,ze17);
% 
% events_averg.all.ed_nc_averg = cat(2,events_averg.all.ed_nc_averg(:,1:74),ze35, events_averg.all.ed_nc_averg(:,75:122),ze72, ze24,events_averg.all.ed_nc_averg(:,123:end));
% 
% events_averg.all.ee_averg = cat(2,events_averg.all.ee_averg(:,1:74),ze35,events_averg.all.ee_averg(:,75:end));
% events_averg.all.expl_averg = cat(2,events_averg.all.expl_averg(:,1:74),ze35,events_averg.all.expl_averg(:,75:end));



%% binary
% std of the population to compare to 
std_ed = std(events_averg.all.ed(1:10,:), 1);
std_ef = std(events_averg.all.ef(1:10,:), 1);
std_es = std(events_averg.all.es(1:10,:), 1);
std_ee = std(events_averg.all.ee(1:10,:), 1);
std_er = std(events_averg.all.er(1:10,:), 1);
std_be = std(events_averg.all.be(1:10,:), 1);

fn3 = fieldnames(events_averg.all);fact=1;
events=[];
column= 0;
for hh=1:(numel(fn3))
    column= column+1;
    for kk=1:size()
        if events_averg.all.(fn3{hh})(1,kk) == 0 %if trial does not exist, dont mark as cell
            events(column,kk) = 3;
        elseif events_averg.all.(fn3{hh})(1,kk) > (perc95(hh,kk)*fact)
            events(column,kk) = 1;
        elseif events_averg.all.(fn3{hh})(1,kk) < (perc5(hh,kk)*fact) 
            events(column,kk) = 2;
        elseif events_averg.all.(fn3{hh})(1,kk) ~= 0  % if its anything else but exact 0 (no trials)   
            events(column,kk) = 0;
        end
    end
end
% modify order of data files so they match Perc files
fn1=fieldnames(events_averg.all);
fn= fn1(20:41);
fn3 = cat(1,fn(1:5), fn(7:8), fn(10:14),fn(16),fn(18:19),fn(9),fn(20:22),fn(6),fn(15));
for i = 1:numel(fn3)
    events_averg.allSort.(fn3{i})= events_averg.all.(fn3{i});
end    


fnallSort = fieldnames(events_averg.allSort);

% compare to shuffled thresholds
perc95 = table2array(Perc95);
perc5 = table2array(Perc5);

fn3 = fieldnames(events_averg.allSort);fact=2;
events=[];
column= 0;
for hh=1:(numel(fn3))
    column= column+1;
    for kk=1:115
        if events_averg.allSort.(fn3{hh})(1,kk) == 0 %if trial does not exist, dont mark as cell
            events(column,kk) = 3;
        elseif events_averg.allSort.(fn3{hh})(1,kk) > (perc95(hh,kk)*fact)
            events(column,kk) = 1;
        elseif events_averg.allSort.(fn3{hh})(1,kk) < (perc5(hh,kk)*fact) 
            events(column,kk) = 2;
        elseif events_averg.allSort.(fn3{hh})(1,kk) ~= 0  % if its anything else but exact 0 (no trials)   
            events(column,kk) = 0;
        end
    end
end

[r,c]=find(events==1);%finds ones excited cells
[r2,c2]=find(events ==2);%finds 2 inhibited cells

fn4 = fieldnames(events_averg.allSort);
for ev = 1:numel(fn4)
    events_averg.cells.(char(strrep((fn4{ev}), '_averg', '_cells'))) = c(find(r==ev));
    events_averg.cells.(char(strrep((fn4{ev}), '_averg', '_cells_inh'))) = c2(find(r2==ev));
end

%extremeCells = events_averg.cells;
%find neuropil cells and remove from overall cell count
events_averg.cellsNP = events_averg.cells;
NpCells = [28 49 61 115];
NpC = [12 66]; %only g2 and g12 trials run and enter run
fn =fieldnames(events_averg.cells);
for i= 1:numel(fn)
    evt = events_averg.cells.(fn{i});
    if i <39 
        [X,Y] = ismember(evt,NpCells);
        evt(Y(X)) = [];
        events_averg.cells.(fn{i}) = evt;
    elseif i > 38
        [X,Y] = ismember(evt,NpC);
        evt(Y(X)) = [];
        events_averg.cells.(fn{i}) = evt;
    end    
end 



% as vars in workspace
ed_cells =c(find(r==1));ef_cells =c(find(r==2));es_cells =c(find(r==3));ee_cells =c(find(r==5));er_cells =c(find(r==20));
drink_antic_cells =c(find(r==18));drink_cells =c(find(r==6));eat_antic_cells =c(find(r==19));eat_cells =c(find(r==7));
be_cells =c(find(r==16));home_cells =c(find(r==17));soc_cells=c(find(r==4));expl_cells=c(find(r==13));run_cells=c(find(r==21));
%exd_cells =c(find(r==8));exf_cells =c(find(r==9));exs_cells =c(find(r==10));exe_cells =c(find(r==11));exr_cells =c(find(r==12));
ef_nc_cells=c(find(r==14));er_nc_cells=c(find(r==15));


ed_cells_inh =c2(find(r2==1));ef_cells_inh =c2(find(r2==2));es_cells_inh =c2(find(r2==3));ee_cells_inh =c2(find(r2==5));er_cells_inh =c2(find(r2==20));
drink_antic_cells_inh =c2(find(r2==18));drink_cells_inh =c2(find(r2==6));eat_antic_cells_inh =c2(find(r2==19));eat_cells_inh =c2(find(r2==7));
be_cells_inh =c2(find(r2==16));home_cells_inh =c2(find(r2==17));soc_cells_inh=c2(find(r2==4));expl_cells_inh=c2(find(r2==13));run_cells_inh=c2(find(r2==21));
%exd_cells_inh =c2(find(r2==8));exf_cells_inh =c2(find(r2==9));exs_cells_inh =c2(find(r2==10));exe_cells_inh =c2(find(r2==11));exr_cells_inh =c2(find(r2==12));
ef_nc_cells_inh=c2(find(r2==14));er_nc_cells_inh=c2(find(r2==15));


find(events_averg.allSort.eat_averg>0.5);

%% percentages

events_averg.cellswithExit = events_averg.cells;
%remove exits from cell list
fnC = fieldnames(events_averg.cells);
remF = [15:24];
events_averg.cells = rmfield(events_averg.cells, fnC(remF));
all=115;
fn1=fieldnames(events_averg.cells);
for j=1:numel(fn1)
    events_averg.perc.(fn1{j}) = [(size(events_averg.cells.(fn1{j}),1)/all)*100 ; size(events_averg.cells.(fn1{j}),1)];
end
%% modify events array for counting, remove exits
events2 = events;events2(events2 == 3) =0;% replace all 3 (non existant trials)
events2(8:12,:) = [];
noEventCells = find((sum(events2,1)) == 0); size(noEventCells,2);

%% Specialization score
specScore = [];
for evnt = 1:size(events2,2)
    specNr = 0;
    for cl = 1:size(events2,1)
        if events2(cl,evnt) ~= 0 
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



% % plot specialization score
% figure; bar(NrSpecCells(2,:));hold on;
% xticks([1:size(dd,2)]);xticklabels(NrSpecCells(1,:));xlabel('\fontsize{12}Nr of encoded events');
% ylabel('\fontsize{12}Nr of cells');hold off;title('\fontsize{11}Specialization score')
% ax = gca;ax.FontSize = 10;
% %specializedCells = (find(specScore == 1));


%% piechart event cells, specialized and no event cells
figure;label={' \fontsize{12}no event cells',' \fontsize{12}event cells', '\fontsize{12}specialized cells' };pieCells = [NrSpecCells(2,1), 115-(NrSpecCells(2,1) + NrSpecCells(2,2)), NrSpecCells(2,2)];
pie(pieCells, '%.1f%%');lgd = legend(label);

%% idendity of specialized cells

specializedCells = (find(specScore == 1));

for i = 1:numel(specializedCells);
    eventtype = find(events2(:,specializedCells(1,i)) ~= 0);
    specializedCells(2,i) = eventtype;
end
%event numbers
%eventlist = fieldnames(events_averg.allSort);eventlist(8:12,:) = [];

specC = unique(specializedCells(2,:));
for ii = 1:numel(specC)
    spec{:,ii} = specializedCells(1,find(specializedCells(2,:)==specC(1,ii)));% cells which encode events in specC
end

%%plot

% plot specialization score
figure; 
subplot(1,2,1);bar(NrSpecCells(2,:));hold on;
xticks([1:size(dd,2)]);xticklabels(NrSpecCells(1,:));xlabel('\fontsize{12}Nr of encoded events');
ylabel('\fontsize{12}Nr of cells');title('\fontsize{11}Specialization score')
ax = gca;ax.FontSize = 10;
%y= [size(spec{1,1},2) size(spec{1,2},2) size(spec{1,3},2) size(spec{1,4},2) size(spec{1,5},2) size(spec{1,6},2)  size(spec{1,8},2) size(spec{1,9},2) (size(spec{1,10},2)+size(spec{1,7},2))];

y=[7 4 4 3 2 1 1 1 1];
x=[1:numel(specC)-1] 
subplot(1,2,2); bar(x,y);
%;xticklabels({'enter drink', 'enter food', 'social', 'drink', 'eat','explore','nest' , 'eat antic', 'enter run' });
xticklabels({'nest','social','explore', 'eat', 'enter run','enter food','eat antic','enter drink','drink' });
%ylabel('\fontsize{12}Nr of cells');
title('\fontsize{11}Specialized Cells')
ylim([0 8]);
ax = gca;ax.FontSize = 10;
hold off;

%pi chart with only specialized cells
figure;label={'nest','social','explore', 'eat', 'enter run','enter food','eat antic','enter drink','drink' };pieCells = [7, 4, 4, 3, 2, 1, 1, 1, 1];
pie(pieCells, '%.1f%%');lgd = legend(label);


%% venn diagrams
% A = 1:115;
% %excitatory
exCons = unique(cat(1,eat_cells, drink_cells, expl_cells, soc_cells, home_cells, run_cells));
InhConsum = unique(cat(1,eat_cells_inh, drink_cells_inh, expl_cells_inh, soc_cells_inh, home_cells_inh, run_cells_inh));
BOTH = unique(cat(1,exCons,InhConsum));
totalcons = unique(cat(1,BOTH, exCons, InhConsum));size(totalcons);%98-56 =42!!




BothApp = unique(cat(1,exApp,InhApp));
totalapp = unique(cat(1,BothApp, exApp, InhApp));size(totalapp);

overl = ismember(BOTH,BothApp);size(find(overl ==1));;

%
overl=ismember(exCons,exApp);size(find(overl ==1))%13 cells are both ex Cons and exApp
overl=ismember(exCons,InhConsum);size(find(overl ==1))%23 cella are both exCons but also InhCons 
overl=ismember(exApp,InhConsum);size(find(overl ==1))%20 cells are both exApp and Inh cons
% 13+23-20 = 16cells to subtract from 

% exApp = unique(cat(1,ef_cells, ed_cells, ee_cells, es_cells, be_cells, er_cells));
% InhApp = unique(cat(1,ef_cells_inh, ed_cells_inh, ee_cells_inh, es_cells_inh, be_cells_inh, er_cells_inh));
% 
% setListData2={BOTH BothApp};
% setLabels2=[ " ";" "];
% hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);
% hh.ShowIntersectionCounts = true;

% 
% %drink module
% setListData2={A ed_cells drink_antic_cells drink_cells};
% setLabels2=[ " " ;"entry drink ";"drink antic";"drink"];
% hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);
% hh.ShowIntersectionCounts = true;
% 
% size(find(ismember(ed_cells, drink_antic_cells, "rows")),1)
% size(find(ismember(ed_cells, drink_cells, "rows")),1)
% size(find(ismember(drink_cells, drink_antic_cells, "rows")),1)
% %eat module
% size(find(ismember(ef_cells, eat_antic_cells, "rows")),1)
% size(find(ismember(ef_cells, eat_cells, "rows")),1)
% size(find(ismember(eat_cells, eat_antic_cells, "rows")),1)

%% bar plots 
% entries
Entries = [size(ed_cells,1) size(ed_cells_inh,1); size(ef_cells,1) size(ef_cells_inh,1); size(es_cells,1)  size(es_cells_inh,1);size(er_cells,1) size(er_cells_inh,1);size(ee_cells,1) size(ee_cells_inh,1);size(be_cells,1) size(be_cells_inh,1)];
X = categorical({'drink','food','social','run','explore','end'});X = reordercats(X,{'drink','food','social','run','explore', 'end'});
figure; bar(X,Entries,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([1 60]);

Entries = [size(ed_cells,1) size(ed_cells_inh,1); size(ef_cells,1) size(ef_cells_inh,1); size(es_cells,1)  size(es_cells_inh,1);size(er_cells,1) size(er_cells_inh,1);size(ee_cells,1) size(ee_cells_inh,1);size(be_cells,1) size(be_cells_inh,1)];
ylabel({'Entries'});
yticks([1:6]); yticklabels({'drink','food','social','run','explore','end'});
figure; b1= barh(Entries,'stacked');ylim([0.5 6.5]);xlim([0 26]);%


%% consumption %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Cons = [size(drink_cells,1) size(drink_cells_inh,1); size(eat_cells,1) size(eat_cells_inh,1); size(soc_cells,1)  size(soc_cells_inh,1); size(run_cells,1)  size(run_cells_inh,1)     ;size(expl_cells,1) size(expl_cells_inh,1);size(home_cells,1) size(home_cells_inh,1)];
X = categorical({'drink','eat','social','run','explore','home'});X = reordercats(X,{'drink','eat','social','run','explore','home'});
figure; bar(X,Cons,'stacked'); ylabel('\fontsize{12}Nr of cells');ylim([1 60]);

yticks([1:6]); yticklabels({'drink','food','social','run','explore','end'});
figure; b1= barh(Cons,'stacked');ylim([0.5 6.5]);xlim([0 61]);%

%% compare nr of cells before entry vs consum


ExConsC = unique(cat(1,eat_cells, drink_cells, expl_cells, soc_cells, home_cells, run_cells));
InhConsC = unique(cat(1,eat_cells_inh, drink_cells_inh, expl_cells_inh, soc_cells_inh, home_cells_inh, run_cells_inh));
ExEntryC = unique(cat(1,ef_cells, ed_cells, ee_cells, es_cells, be_cells, er_cells, eat_antic_cells, drink_antic_cells));
InhEntryC = unique(cat(1,ef_cells_inh, ed_cells_inh, ee_cells_inh, es_cells_inh, be_cells_inh, er_cells_inh,eat_antic_cells_inh, drink_antic_cells_inh));


Conall = unique(cat(1,ExConsC, InhConsC));size(Conall,1)%all consuption incl overlap
Appall =unique(cat(1,ExEntryC, InhEntryC));size(Appall,1)%all appetitive including overlap

overlap= ismember(Conall,Appall); overl = Conall(find(overlap ==0));size(overl,1)%all consumption without appetitive
overlap2= ismember(Appall,Conall); overl2 = Appall(find(overlap2 ==0));size(overl2,1)%all appetitve without consumption

exConsonly = ismember(ExConsC, InhConsC);exConsonly = ExConsC(find(exConsonly ==0));size(exConsonly,1)% exit cons without inh cons
exConsonly2 = ismember(exConsonly, Appall);exConsonly2 = exConsonly(find(exConsonly2 ==0));size(exConsonly2,1)% exit cons without inh cons and without app




setListData2={1:111 Conall Appall };
setLabels2=[ " " ;"consumption ";"appetitive"];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);
hh.ShowIntersectionCounts = true;

exCons = unique(cat(1,eat_cells, drink_cells, expl_cells, soc_cells, home_cells, run_cells));
InhConsum = unique(cat(1,eat_cells_inh, drink_cells_inh, expl_cells_inh, soc_cells_inh, home_cells_inh, run_cells_inh));
BOTH = unique(cat(1,exCons,InhConsum));
size(find(ismember(exCons, InhConsum) == 1)); 
totalcons = unique(cat(1,BOTH, exCons, InhConsum));size(totalcons);%98-56 =42!!

exApp = unique(cat(1,ef_cells, ed_cells, ee_cells, es_cells, be_cells, er_cells, eat_antic_cells, drink_antic_cells));
InhApp = unique(cat(1,ef_cells_inh, ed_cells_inh, ee_cells_inh, es_cells_inh, be_cells_inh, er_cells_inh, eat_antic_cells_inh,drink_antic_cells_inh));
size(find(ismember(exApp, InhApp) == 1))

 x = 2020; 
 y = [98-69 69 72-69]; 
 b = bar(x,y,"stacked");hold on
 x2= 2021
 y2= [size(InhConsum,1)-23 23 56];
 b2 = bar(x2,y2,"stacked");
 x3= 2022
 y3= [size(InhApp,1)-29 29 size(exApp,1)-29 ];
 b3 = bar(x3,y3,"stacked");

pie(y2);
pie(y3);

%
pie only


isboth = ismember(ExConsC,ExEntryC);ExBoth = ExConsC(isboth == 1);
isExConsConly = ismember(ExConsC, ExBoth); ExConsConly = ExConsC(isExConsConly == 0); 

ExConsnoApp = ismember(ExConsConly, InhBoth); ExConsnoApp = ExConsnoApp(ExConsnoApp == 0);% cells which are ex cons but not inh cons
%and also not inh app!!

isExEntryConly = ismember(ExEntryC, ExBoth); ExEntryConly = ExEntryC(isExEntryConly == 0);

isboth = ismember(InhConsC,InhEntryC);InhBoth = InhConsC(isboth == 1);
isInhConsConly = ismember(InhConsC, InhBoth); InhConsConly = InhConsC(isInhConsConly == 0); 
isInhEntryConly = ismember(InhEntryC, InhBoth); InhEntryConly = InhEntryC(isInhEntryConly == 0);

%ystacked = [size(InhEntryConly,1) size(InhBoth,1) (size(InhConsConly,1)) ;size(ExEntryConly,1) size(ExBoth,1) (size(ExConsConly,1))];
ystacked = [(size(InhEntryConly,1) + size(ExEntryConly,1)) (size(InhBoth,1)+ size(ExBoth,1)) (size(InhConsConly,1) + size(ExConsConly,1))]; 
figure; bar(ystacked,'stacked'); hold on; %yticklabels({'Inh', 'Exc'});clear ylabel;xlabel('\fontsize{12}Nr of cells');
label={' \fontsize{12}Appetitive', '\fontsize{12} Both',' \fontsize{12}Consummatory'};lgd = legend(label);


%venn R
%Inh: A =115 B =24 C=11 both =30
%Exc: A =115 B =11 C=42 both =37




% overlap greater than chance?
A= size(InhEntryC,1)/115;
B= size(InhConsC,1)/115;
ABobserved = size(InhBoth,1)/115;
AnBexpected = (A*B)

ABobserved*115
AnBexpected*115

A= size(ExEntryC,1)/115;
B= size(ExConsC,1)/115;
ABobserved = size(ExBoth,1)/115
AnBexpected = (A*B)



%% all trials from one cell
fnEA = fieldnames(events_averg);
anim = 4;
fnEAan = fieldnames(events_averg.(fnEA{anim}));
fnEAanimal = fnEAan(20:38,1);% baselined individual traces
fnEAanimalNP = fnEAan(39:end,1);
Event = 8; CellID =15;
%
Cell1 = events_averg.(fnEA{anim}).(fnEAanimal{Event})(:,CellID,:);
C1 = reshape(Cell1 ,[],size(Cell1,3));
%NP
Cell1NP = events_averg.(fnEA{anim}).(fnEAanimalNP{Event});
C1NP = reshape(Cell1NP ,[],size(Cell1NP,3));
%entry
% figure; plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
% xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time from event[s]');xlim([1 31]);
% ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');%title(['\fontsize{10}entry drink inhibited cell']);

%cons
figure; plot(C1, LineWidth=1);hold on; plot((mean(C1,2)),'k' ,LineWidth=4);
xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1', '0', '1', '2', '3'});xlabel('\fontsize{12}time from event[s]');xlim([1 61]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');%title(['\fontsize{10}entry drink inhibited cell']);

%% averg traces from multiple cells
fnAll = fieldnames(events_averg.all);
Event = 8;
%Cells = [1 35 38 66 76];
Cells = [35 38 66 76 114];

figure; plot(events_averg.all.(fnAll{Event}) (:, Cells), LineWidth=2);hold on; xlim([1 30]); 
xticks([1 5 10 15 20 25 30]); xticklabels({'-3','-2' ,'-1','0' ,'1' , '2', '3'});xlabel('\fontsize{12}time from event [s]'); ylabel('\fontsize{12}averg Ca-activity (z-scored, s.d.)');
%plot([1 10], [-0.9 -0.9],"k", LineWidth=3);plot([15 30], [-0.9 -0.9],"k", LineWidth=3);
%text(2.5,-0.8,'baseline');text(16,-0.8,'event window');
legend({'\fontsize{10}c5','\fontsize{10}c4', '\fontsize{10}c3','\fontsize{10}c2','\fontsize{10}c1'}, Location="northwest");

%% averg traces of event cells
%cells = [excited no event inhibited soc cell]
event = "ed_averg_bsl"; tit = ['\fontsize{14}enter drink'];cells =[50 64 36 24];
%event = "ef_averg_bsl"; tit = ['\fontsize{14}enter food'];cells =[14 64 53 24];
%event = "es_averg_bsl"; tit = ['\fontsize{14}enter social'];cells =[50 64 18 24];
%event = "er_averg_bsl"; tit = ['\fontsize{14}enter run'];cells =[66 64 78 24];
%event = "ee_averg_bsl"; tit = ['\fontsize{14}enter explore'];cells =[74 6 24];%no inhib explore cells
%event = "blockend_averg_bsl"; tit = ['\fontsize{14}block end'];cells =[102 64 66 24];
% 
%event = "drink_averg_bsl"; tit = ['\fontsize{14}drink'];cells =[51 64 102 24];
%event = "eat_averg_bsl"; tit = ['\fontsize{14}eat'];cells =[76 64 14 24];
%event = "soc_averg_bsl"; tit = ['\fontsize{14}social'];cells =[82 64 50 24];
%event = "expl_averg_bsl"; tit = ['\fontsize{14}explore'];cells =[47 64 19 24];
%event = "run_averg_bsl"; tit = ['\fontsize{14}run'];cells =[33+49 64+49 38+49];%24 does not exist
%event = "blockend_averg_bsl"; tit = ['\fontsize{10}nest'];cells =[52 64 66 24];

color = { 'r', 'k', 'b', 'w'};

for i =1:numel(cells)
    cellNr =cells(i) ; 
    if cellNr < 29
        animal = 'g5';plus = 0;
    elseif cellNr < 50
        animal = "g4";plus = 28;
    elseif  cellNr <  62
        animal = "g2";plus = 49;
    elseif  cellNr > 61
        animal = "g12";plus = 61;
    end
cell = events_averg.(char(animal)).(char(event))(:,cellNr-plus,:); 
%cell =events_averg.(char(animal)).(char(event))(21:101,cellNr-plus,:);%BLOCKEND
cNr = reshape(cell,[],size(cell,3));
 
sp1 = stdshade(cNr',0.3,color{i});hold on;

        %ylimits = [-3 6];xlimits = [1 61];
        %xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');

        %ylimits = [-3 6];xlimits = [1 41];
        %xticks([1 11 21 31 41]);xticklabels({'-1','0','1','2','3'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');

        %ylimits = [-3 6];xlimits = [1 81];
        %xticks([1 11 21 31 41 51 61 71 81]);xticklabels({'-1','0','1','2','3', '4', '5','6', '7'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');
        
        ylimits = [-3 6];xlimits = [1 31];
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)'); 
        
        ax = gca;ax.FontSize = 10;
        ylim(ylimits);xlim(xlimits);
        title(tit);
        
end
hold off;

%Position settings
%nest 0.13,0.425,0.1360,0.15
%drink 0.13,0.425,0.102,0.15
%social 0.13,0.425,0.068,0.15
%entries 0.13,0.425,0.051,0.15

%% raw traces of cells
enter = 5217;%ef 2361 pellet 2303 leave social
s= 9; cells = [107-61 101-61 78-61];
raw_f_res =ci_data.bsl.g12.d_2023_03_01(s).session.raw_f_res;
Fneu = ci_data.bsl.g12.d_2023_03_01(s).session.Fneu;
subneuro_trace = subneuro(raw_f_res,Fneu,0.8); % subtracts neuropil trace from raw trace by factor x subneuro(rawF, Fneuropil, factor)
subneuro_trace(:,(size(subneuro_trace,2)+1)) = mean(Fneu,2);%last "cell" is mean neuropil trace
df_f_trace = dFoF(subneuro_trace);%df/f
raw_f= zscore(df_f_trace,0,1);% zscored using sample sd

c1 = df_f_trace(enter-200:enter+100,cells(1));z_c1 = (c1-nanmean(c1(1:50,1)))/nanstd(c1(1:50,1));
c2 = df_f_trace(enter-200:enter+100,cells(2));z_c2 = (c2-nanmean(c2(1:50,1)))/nanstd(c2(1:50,1));
c3 = df_f_trace(enter-200:enter+100,cells(3));z_c3 = (c3-nanmean(c3(1:50,1)))/nanstd(c3(1:50,1));

z_c1 = raw_f(enter-200:enter+100,cells(1));
z_c2 = raw_f(enter-200:enter+100,cells(2));
z_c3 = raw_f(enter-200:enter+100,cells(3));

%c4 = df_f_trace(enter-200:enter+100,cells(4));z_c4 = (c4-nanmean(c4(1:50,1)))/nanstd(c4(1:50,1));
%x-nanmean(b)/nanstd(b)
%exmpl_ed = table(smooth(zscore(c1)), smooth(zscore(c2)), smooth(zscore(c3)), smooth(zscore(c4)));
%figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
%plot(smooth(zscore(c40_ed)));hold on;plot(smooth(zscore(c52_ed)));plot(smooth(zscore(c56_ed)));plot(smooth(zscore(c54_ed)));
p1=plot(smooth(z_c1),'LineWidth',5);hold on;p2=plot(smooth(z_c2),'LineWidth',5);p3=plot(smooth(z_c3),'LineWidth',5); %p4=plot(smooth(z_c4));
p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;%p4.LineWidth=1;
xl= xline(200-25,'-','entry food');xl.LineWidth = 1; 
xxl = xline(200+19,'-','eat');xxl.LineWidth = 1;
xxxl = xline(200-93,'-','drink');xxxl.LineWidth = 1;
%xxxxl = xline(200-45,'-','exit drink');xxxxl.LineWidth = 1;
xticks([0+25 50+25 100+25 150+25 200+25 250+25 300+25]);xticklabels({'-15','-10','-5','0','5','10','15'});xlabel('\fontsize{12}time [s]');xlim([0 300]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');ylim([-2 9]);
legend({'\fontsize{12}c1','\fontsize{12}c2', '\fontsize{12}c3'}, Location="northwest");

%% imgsc plots of individual cells
%cellNr =100 ; tit = ['\fontsize{14}drink cell'];%30 42
cellNr =76 ; tit = ['\fontsize{14}eat cell'];%30 42
if cellNr < 29
    animal = 'g5';plus = 0;
elseif cellNr < 50
    animal = "g4";plus = 28;
elseif  cellNr <  62
    animal = "g2";plus = 49;
elseif  cellNr > 61
    animal = "g12";plus = 61;
end    

cell = events_averg.(char(animal)).eat_averg_bsl(:,cellNr-plus,:); 
%cell = events_averg.(char(animal)).ed_averg_bsl([21:101],cellNr-plus,:);
cNr = reshape(cell,[],size(cell,3));
%cNr = cNr(:,1:end);

lim = [-2 8];
fA = figure;
figure(fA);
        s1= subplot(2,1,1);imagesc([1 61],1:size(cNr,2),cNr');clim(lim);xlim([1 61]);hold on;
        title(tit);xticklabels({'','','','','','',''});ylabel('\fontsize{10}Trials');
        set(s1, 'Position', [0.13,0.35,0.28,0.53]); c = colorbar(s1,'Position',[0.43,0.35,0.01,0.53]);
        
        s2= subplot(2,1,2);stdshade(cNr',0.3,[0 0 1]);xlim([1 61]);ylim([-2 4]);
        xticks([1 11 21 31 41 51 61])
        xticklabels({'-3','-2','-1','0', '1', '2','3'});xlabel('\fontsize{12}time [s]')
        ylabel('\fontsize{7}Ca-activity (z-score)');
        set(s2, 'Position', [0.13,0.14,0.28,0.2]);
        

%% imgsc plots of averages
% sort cells acording to amplitude in event window for plotting
SortedEvents = events_averg.allSort;
fn= fieldnames(SortedEvents);
% give all cells a label %sort cells according to value in event window

for i= 1:(numel(fn))
    SortedEvents.(fn{i})=SortedEvents.(fn{i})(SortedEvents.(fn{i}) ~= 0);% removes non existant tirals again
    label = [1:size(SortedEvents.(fn{i}),2)];
    clear temp order;
    SortedEvents.(fn{i})(2,:) = label;
    [temp, order] = sort(SortedEvents.(fn{i})(1,:));% sorted according to evt window value 
    SortedEvents.(fn{i}) = SortedEvents.(fn{i})(:,order);
    SortedEventsList(i,1:(size(SortedEvents.(fn{i}),2)))= SortedEvents.(fn{i})(2,:);
end

% sort averg traces according to SortedEventsList
events_averg.allAsc = events_averg.all;
fn = fieldnames(events_averg.allAsc);
events_averg.allAsc = rmfield(events_averg.allAsc,fn(20:41));events_averg.allAsc = rmfield(events_averg.allAsc,fn(17));
events_averg.allAsc.home = events_averg.allAsc.be;
events_averg.allAsc.drink_antic = events_averg.allAsc.drink;
events_averg.allAsc.eat_antic = events_averg.allAsc.eat;

fn= fieldnames(events_averg.allAsc);
fn3 = cat(1,fn(1:5), fn(7:8), fn(10:14),fn(16:18),fn(9),fn(19),fn(20:21),fn(6),fn(15));
 for i = 1:numel(fn3)
    events_averg.allAscend.(fn3{i})= events_averg.allAsc.(fn3{i}); % now allAscend matches SortedEventsList
 end    
% move neuropil traces to beginning of SortedEventsList
NpCells = [28 49 61 115];
NpC = [12 66]; %only g2 and g12 trials
for i= 1:21
    evt = SortedEventsList(i,:);
    if i <20 
        [X,Y] = ismember(evt,NpCells);
        evt(Y(X)) = [];
        NPbegin = cat(2,NpCells,evt);% Np cells at beginning
        SortedEventsList(i,:) = NPbegin;
    elseif i > 19
        [X,Y] = ismember(evt,NpC);
        evt(Y(X)) = [];
        NPbegin = cat(2,NpC,evt);% Np cells at beginning
        SortedEventsList(i,:) = NPbegin;
    end    
end

%sort allAllscend in ascending order according to SortedEventsList
fn =fieldnames(events_averg.allAscend);
for i = 1:(numel(fn))
    tosort = events_averg.allAscend.(fn{i});
    if i <20
        for j = 1:(size(SortedEventsList(i,:),2))
        sorted(:,j) = tosort(:,(SortedEventsList(i,j)));
        end 
    elseif i > 19
        for k = 1:66
        sorted(:,k) = tosort(:,(SortedEventsList(i,k)));
        end 
    end    
    events_averg.allAsc.(fn{i}) = sorted;% averaged traces of all cells sorted in ascending order
    clear sorted
end    


%% imgsc entries
win= 30;
lim = [-3 3];
limaverg = [-1 1];
f7=figure;% averages across all sessions
figure(f7);hold on;
sp11 = subplot(2,6,1),imagesc([-win:0],[1:size(events_averg.allAsc.ed,2)],events_averg.allAsc.ed');xline(-20,'--b');yline(4.5,'-w');left=0.13; width=0.051; btw= 0.0325;
        %trials = sort(table2array(Trialnumbers(1,:)));trials(trials == 0) = []; 
        set(sp11, 'Position',[left,0.58,width,0.34]);
        %title(['\fontsize{10}enter drink',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);
        clim(lim);
        xticklabels({' ',' ',' ',''});ylabel('\fontsize{15}Cells'); 
        
        sp1= subplot(2,6,7); stdshade(events_averg.allAsc.ed',0.3,[0 0 1]);set(sp1, 'Position', [left,0.425,width,0.15]);hold on
            %stdshade(events_averg_ch.allAsc.ed',0.3,[1 1 1]);
            xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-scored)');
            ylim(limaverg);hold off

sp22 = subplot(2,6,2),imagesc([-win:0],[1:size(events_averg.allAsc.ef,2)],events_averg.allAsc.ef');xline(-20,'--b');yline(4.5,'-w');left =left+width+btw;
        %trials = sort(table2array(Trialnumbers(2,:)));trials(trials == 0) = []; 
        set(sp22, 'Position',[left,0.58,width,0.34]);
        %title(['\fontsize{10}enter food',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);
        clim(lim);
        xticklabels({' ',' ',' ',''});
        sp2= subplot(2,6,8); stdshade(events_averg.allAsc.ef',0.3,[0 0 1]);set(sp2, 'Position', [left,0.425,width,0.15]);hold on
        %stdshade(events_averg_ch.allAsc.ef',0.3,[1 1 1]);hold off
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xtickangle(0);
        ylim(limaverg);
sp33= subplot(2,6,3),imagesc([-win:0],[1:size(events_averg.allAsc.es,2)],events_averg.allAsc.es');xline(-20,'--b');yline(4.5,'-w');left =left+width+btw;
        %trials = sort(table2array(Trialnumbers(3,:)));trials(trials == 0) = []; 
        set(sp33, 'Position',[left,0.58,width,0.34]);
        %title(['\fontsize{10}enter social',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);
        clim(lim);
        xticklabels({' ',' ',' ',''});xlabel('\fontsize{15}time [s]');
        sp3= subplot(2,6,9); stdshade(events_averg.allAsc.es',0.3,[0 0 1]);set(sp3, 'Position', [left,0.425,width,0.15]);hold on
        %stdshade(events_averg_ch.allAsc.es',0.3,[1 1 1]);hold off
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xtickangle(0);
        ylim(limaverg);
sp44=subplot(2,6,4),imagesc([-win:0],[1:size(events_averg.allAsc.ee,2)],events_averg.allAsc.ee');xline(-20,'--b');yline(4.5,'-w');left =left+width+btw;
        %trials = sort(table2array(Trialnumbers(5,:)));trials(trials == 0) = []; 
        set(sp44, 'Position',[left,0.58,width,0.34]);
        %title(['\fontsize{10}enter explore',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);
        clim(lim);
         xticklabels({' ',' ',' ',''});
        sp4= subplot(2,6,10); stdshade(events_averg.allAsc.ee',0.3,[0 0 1]);set(sp4, 'Position', [left,0.425,width,0.15]);hold on
        %stdshade(events_averg_ch.allAsc.ee',0.3,[1 1 1]);hold off
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xtickangle(0);
        ylim(limaverg);
sp55=subplot(2,6,5),imagesc([-win:0],[1:size(events_averg.allAsc.er,2)],events_averg.allAsc.er');xline(-20,'--b');yline(2.5,'-w');left =left+width+btw;
        %trials = sort(table2array(Trialnumbers(6,:)));trials(trials == 0) = []; 
        set(sp55, 'Position',[left,0.58,width,0.34]);
        title(['\fontsize{10}enter run',' trials: 26 ']);clim(lim);
        xticklabels({'-3','-2','-1',' '});
        sp5= subplot(2,6,11); stdshade(events_averg.allAsc.er',0.3,[0 0 1]);set(sp5, 'Position', [left,0.425,width,0.15]);hold on
        %stdshade(events_averg_ch.allAsc.er',0.3,[1 1 1]);hold off
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xtickangle(0);
        ylim(limaverg);
sp66=subplot(2,6,6),imagesc([-win:0],[1:size(events_averg.allAsc.be(1:31,:),2)],events_averg.allAsc.be(1:31,:)');xline(-20,'--b');yline(4.5,'-w');left =left+width+btw;
        %trials = sort(table2array(Trialnumbers(9,:)));trials(trials == 0) = []; 
        set(sp66, 'Position',[left,0.58,width,0.34]);
        %title(['\fontsize{10}block end',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);
        clim(lim);
        xticks([-30 -20 -10 0]); xticklabels({' ',' ',' ',' '});
        sp6= subplot(2,6,12); stdshade((events_averg.allAsc.be(1:31,:))',0.3,[0 0 1]);set(sp6, 'Position', [left,0.425,width,0.15]);hold on
        %stdshade((events_averg_ch.allAsc.be(1:31,:))',0.3,[1 1 1]);hold off
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xtickangle(0);
        ylim(limaverg);
        h = axes(f7,'visible','off');left =left+width+btw;
        c = colorbar(h,'Position',[left 0.425 0.01 0.49]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 4;
        caxis(h,lim);


%% imgsc consumption
win= 30;
lim = [-3 3];limaverg2= [-1 1];
f9=figure;% averages across all sessions [xLeft, yBottom, width, height]
figure(f9);hold on;

sp1 = subplot(2,6,1);imagesc([-win:win],[1:size(events_averg.allAsc.drink,2)],events_averg.allAsc.drink');xline(1,'--w');xline(-20,'--b');yline(4.5,'-w');left=0.13; btw =0.0325; width1 = 0.0170;
        trials = sort(table2array(Trialnumbers(7,:)));trials(trials == 0) = [];set(sp1, 'Position',[left,0.58,width1*6,0.34]);
        title(['\fontsize{10}drink',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);clim(lim)
        xticks([1 11 21 31 41 51 61]); xticklabels({'-3','-2','-1','0','1','2','3'});;ylabel('\fontsize{15}Cells');xlabel('\fontsize{15}time [s]');
        sp7= subplot(2,6,7); stdshade(events_averg.allAsc.drink',0.3,[0 0 1]);set(sp7, 'Position',[left,0.425,width1*6,0.15]);hold on
        stdshade(events_averg_ch.allAsc.drink',0.3,[1 1 1]);hold off
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');
        ylim(limaverg2);yticks([-1 -0.5 0 0.5 1]);yticklabels({'','-0.5','0','0.5',''});xtickangle(0);

sp2 =subplot(2,6,2);imagesc([-win:win],[1:size(events_averg.allAsc.eat,2)],events_averg.allAsc.eat');xline(1,'--w');yline(4.5,'-w');xline(-20,'--b');left =left+width1*6+btw;
        trials = sort(table2array(Trialnumbers(8,:)));trials(trials == 0) = [];set(sp2, 'Position',[left,0.58,width1*6,0.34]);
        title(['\fontsize{10}eat',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);clim(lim)
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});
        sp8= subplot(2,6,8); stdshade(events_averg.allAsc.eat',0.3,[0 0 1]);set(sp8, 'Position', [left,0.425,width1*6,0.15]);hold on
        stdshade(events_averg_ch.allAsc.eat',0.3,[1 1 1]);hold off
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});
        ylim(limaverg2);yticks([-1 -0.5 0 0.5 1]);yticklabels({'','-0.5','0','0.5',''});xtickangle(0);

sp3 =subplot(2,6,3);imagesc([1:41],[1:size(events_averg.allAsc.soc,2)],events_averg.allAsc.soc');yline(4.5,'-w');xline(11,'--b'); xline(11.2,'--w');left =left+width1*6+btw;
        trials = sort(table2array(Trialnumbers(4,:)));trials(trials == 0) = [];set(sp3, 'Position',  [left,0.58,width1*4,0.34]);
        title(['\fontsize{10}social',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);clim(lim)
        xticks([1 11 21 31 41]); xticklabels({' ',' ','1','2',' '});
        sp9= subplot(2,6,9); stdshade(events_averg.allAsc.soc',0.3,[0 0 1]);set(sp9, 'Position',  [left,0.425,width1*4,0.15]);hold on
        stdshade(events_averg_ch.allAsc.soc',0.3,[1 1 1]);hold off
        xlim([1 41]) ;xticks([1 11 21 31 41]);xticklabels({'-1','0','1','2','3'}); 
        ylim(limaverg2);yticks([-1 -0.5 0 0.5 1]);yticklabels({'','-0.5','0','0.5',''});xtickangle(0);


sp4 =subplot(2,6,4);imagesc([1:41],[1:size(events_averg.allAsc.expl,2)],events_averg.allAsc.expl');yline(4.5,'-w');xline(11,'--b'); xline(11.2,'--w');left =left+width1*4+btw;
        trials = sort(table2array(Trialnumbers(16,:)));trials(trials == 0) = [];set(sp4, 'Position',[left,0.58,width1*4,0.34]);
        title(['\fontsize{10}explore',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);clim(lim); 
        xticks([1 11 21 31 41]); xticklabels({' ',' ','1','2',' '});
        sp10= subplot(2,6,10); stdshade(events_averg.allAsc.expl',0.3,[0 0 1]);set(sp10, 'Position',[left,0.425,width1*4,0.15]);hold on
        stdshade(events_averg_ch.allAsc.expl',0.3,[1 1 1]);hold off
        xlim([1 41]) ;xticks([1 11 21 31 41]);xticklabels({'-1','0','1','2','3'}); 
    ylim(limaverg2);yticks([-1 -0.5 0 0.5 1]);yticklabels({'','-0.5','0','0.5',''});  xtickangle(0);  


sp5 =subplot(2,6,5);imagesc([-win:win],[1:size(events_averg.allAsc.run(:,[2 13:end]),2)],events_averg.allAsc.run(:,[2 13:end])');xline(1,'--w');yline(1.5,'-w');xline(-20,'--b');left =left+width1*4+btw;
        trials = sort(table2array(Trialnumbers(15,:)));trials(trials == 0) = [];set(sp5, 'Position',[left,0.58,width1*6,0.34]);
        title(['\fontsize{10}run',' trials: ', num2str(trials(1,1))]);clim(lim)
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});
        sp11= subplot(2,6,11); stdshade(events_averg.allAsc.run(:,[2 13:end])',0.3,[0 0 1]);set(sp11, 'Position', [left,0.425,width1*6,0.15]);hold on
        stdshade(events_averg_ch.allAsc.run',0.3,[1 1 1]);hold off
        xlim([1 61]);xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});
        ylim(limaverg2);yticks([-1 -0.5 0 0.5 1]);yticklabels({'','-0.5','0','0.5',''});xtickangle(0);

sp6 =subplot(2,6,6);imagesc([1:81],1:size(events_averg.allAsc.home,2),events_averg.allAsc.home(21:101,:)');yline(4.5,'-w');xline(11,'--w');left =left+width1*6+btw;
        trials = sort(table2array(Trialnumbers(9,:)));trials(trials == 0) = [];set(sp6, 'Position',[left,0.58,width1*8,0.34]);
        title(['\fontsize{10}nest',' trials: ', num2str(trials(1,1)), '-' , num2str(trials(1,end)) ]);clim(lim)
        xlim([1 81]);xticks([1 11 21 31 41 51 61 71 81]);xticklabels({' ','2','4','6','8',' ',' ',' ',' '});
        sp12= subplot(2,6,12); stdshade((events_averg.allAsc.home(21:101,:))',0.3,[0 0 1]);set(sp12, 'Position', [left,0.425,width1*8,0.15]);hold on
        stdshade((events_averg_ch.allAsc.home(21:101,:))',0.3,[1 1 1]);hold off
        xlim([1 81]);xticks([1 11 21 31 41 51 61 71 81]);xticklabels({'-1','0','1','2','3','4','5','6','7'});linkaxes([sp6, sp12],'x');
        ylim(limaverg2); yticks([-1 -0.5 0 0.5 1]);yticklabels({'','-0.5','0','0.5',''});xtickangle(0);

        h = axes(f9,'visible','off');left =left+width1*8+btw;
        c = colorbar(h,'Position',[left 0.425 0.01 0.49]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 4;
        caxis(h,lim);        


%% entries consumption vs non consumption trials
lim = [-3 3];limaverg2 = [-0.3 0.5];
f10=figure;
figure(f10);hold on;
subplot(2,4,1),imagesc([-win:0],[1:size(events_averg.allAsc.ed_nc,2)],events_averg.allAsc.ed_nc');
        title(['\fontsize{10}enter no drink',' trials: 42']);clim(lim);
        xticks([1 11 21 31 ]); xticklabels({'-3','-2','-1','0'});;ylabel('\fontsize{15}Cells');xlabel('\fontsize{15}time [s]');
        sp7= subplot(2,4,5); stdshade(events_averg.allAsc.ed_nc',0.3,[0 0 1]);set(sp7, 'Position', [0.13,0.43,0.1566,0.15]);
        xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');
        ylim(limaverg2);
subplot(2,4,2),imagesc([-win:0],[1:size(events_averg.allAsc.ed_c,2)],events_averg.allAsc.ed_c');
        title(['\fontsize{10}enter drink',' trials: 46-82']);clim(lim);
        xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});
        sp8= subplot(2,4,6), stdshade(events_averg.allAsc.ed_c',0.3,[0 0 1]);set(sp8, 'Position', [0.3361,0.43,0.1566,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg2);
subplot(2,4,3),imagesc([-win:0],[1:size(events_averg.allAsc.ef_nc,2)],events_averg.allAsc.ef_nc');
        title(['\fontsize{10}enter no feed',' trials: 9-73']);clim(lim);
        xticks([1 11 21 31]); xticklabels({' ','1','2',' '});
        sp9= subplot(2,4,7); stdshade(events_averg.allAsc.ef_nc',0.3,[0 0 1]);set(sp9, 'Position', [0.5422,0.43,0.1566,0.15]);
        xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg2);
subplot(2,4,4),imagesc([-win:0],[1:size(events_averg.allAsc.ef_c,2)],events_averg.allAsc.ef_c');
        title(['\fontsize{10}enter feed',' trials: 34-97']);clim(lim);
        xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});
        sp10= subplot(2,4,8); stdshade(events_averg.allAscAsc.ef_c',0.3,[0 0 1]);set(sp10, 'Position', [0.7484,0.43,0.1566,0.15]);
       xticks([1 11 21 31 ]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg2);
        h = axes(f10,'visible','off');
        c = colorbar(h,'Position',[0.93 0.45 0.02 0.45]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 2.5;
        caxis(h,lim);        


%% ethogram plots
s=9;animal = "g12";
event_type=table2array(ci_data.bsl.g12.d_2023_03_01(s).session.events(:,4));
raw_f_res =(ci_data.bsl.g12.d_2023_03_01(s).session.raw_f_res);
Fneu = ci_data.bsl.g12.d_2023_03_01(s).session.Fneu;

subneuro_trace = subneuro(raw_f_res,Fneu,0.8); % subtracts neuropil trace from raw trace by factor x subneuro(rawF, Fneuropil, factor)
subneuro_trace(:,(size(subneuro_trace,2)+1)) = mean(Fneu,2);%last "cell" is mean neuropil trace
df_f_trace = dFoF(subneuro_trace);%df/f
raw_f= zscore(df_f_trace,0,1);% zscored using sample sd

x = (table2array(ci_data.bsl.g12.d_2023_03_01(s).session.events(:,5)));
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
ylim([-5 5.5]);
xticks([0 2400 4800 7200 9600 12000]);
xticklabels({'0','5','10','15','20','25'});
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

good_cells = [1 78-61 101-61 107-61 54];
j=0;
for i=good_cells
   subplot(1,1,a),plot((raw_f([4 4 4 4:end],i)/6-j),'Color',[i/size(raw_f,2) 1-i/size(raw_f,2) i/size(raw_f,2)],'LineWidth',1);hold on;
   j=j+1;
end

subplot(1,1,a),plot([300 600],[1 1],'k','LineWidth',2);%0.5s
subplot(1,1,a),plot([300 300],[0.98 0.98+(2/6)],'k','LineWidth',2);%2Std
ylabel("Ca-activity (z-scored)");

%% Averg Time spent in pods (total seconds and percentage?)
fn = fieldnames(Specs);
for ii= 1:4
    for j =1:size(Specs.(fn{ii}).PodTimes,2)
        Times=[];relTimes=[];
        for jj = 1:size(Specs.(fn{ii}).PodTimes,1)
            %actual seconds
            Times = cat(1,Times, Specs.(fn{ii}).PodTimes{jj,j});
            %proportions of whole block time
            TiS = Specs.(fn{ii}).PodTimes{jj,j};TinProp = TiS/Specs.(fn{ii}).PodTimes{jj,6};
            relTimes = cat(1,relTimes, TinProp);
        end
        Sec{:,j}= Times;
        Prop{:,j} = relTimes;
    end
    SEC{ii,:}=Sec;% collect data from all animals 
    PROP{ii,:}=Prop;
end

fn = fieldnames(Specs);
for ii= 1:4
    DecTimes=[];relDecTimes = [];
    for j =1:size(Specs.(fn{ii}).DecTimes,1)
        %Desicion Times
        DecTimes = cat(1,DecTimes, Specs.(fn{ii}).DecTimes{j,1});
        DT = Specs.(fn{ii}).DecTimes{j,1};DTprop = DT/Specs.(fn{ii}).PodTimes{j,6};
        relDecTimes=cat(1,relDecTimes, DTprop);
    end
    SEC{ii,1}{1,7}=DecTimes;% collect data from all animals 
    PROP{ii,1}{1,7}=relDecTimes;
end
secprop = {'food' 'run' 'explore' 'social' 'drink' 'whole block' 'desicion'};

%% trial stats plots
% total time spent in maze and pods
% total time
%g5 g4 g2 g12

%mean drink, food, social, explore, run, desicion
drink = [ mean(SEC{1,1}{1,5}) mean(SEC{2,1}{1,5}) mean(SEC{3,1}{1,5}) mean(SEC{4,1}{1,5})];
meand= mean(drink);
food = [ mean(SEC{1,1}{1,1}) mean(SEC{2,1}{1,1}) mean(SEC{3,1}{1,1}) mean(SEC{4,1}{1,1})];
meanf= mean(food);
social = [ mean(SEC{1,1}{1,4}) mean(SEC{2,1}{1,4}) mean(SEC{3,1}{1,4}) mean(SEC{4,1}{1,4})];
means= mean(social);
explore = [ mean(SEC{1,1}{1,3}) mean(SEC{2,1}{1,3}) mean(SEC{3,1}{1,3}) mean(SEC{4,1}{1,3})];
meane= mean(explore);
run = [ mean(SEC{1,1}{1,2}) mean(SEC{2,1}{1,2}) mean(SEC{3,1}{1,2}) mean(SEC{4,1}{1,2})];
meanr= mean(run);
decs = [ mean(SEC{1,1}{1,7}) mean(SEC{2,1}{1,7}) mean(SEC{3,1}{1,7}) mean(SEC{4,1}{1,7})];
meandec= mean(decs);
totalNr = [size(SEC{1,1}{1,6},1) size(SEC{2,1}{1,6},1)/2 size(SEC{3,1}{1,6},1) size(SEC{4,1}{1,6},1)/2];
meanNr = mean(totalNr);

figure % averg time in each pod  (total Nr of blocks in a day, , points are animals
hold on
y = [ meand meanf means meane meanr meandec];
yy = [drink food social explore run decs];
bar(y)
colors = ["r" "g" "b" "y"];
for i = 1:4
    plot([1 2 3 4 5 6 ],[drink(i) food(i) social(i) explore(i) run(i) decs(i)],[colors(i) + '.'],'MarkerSize', 25);
end
xticks([1 2 3 4 5 6 ]);xticklabels({, 'water', 'food', 'social', 'explore', 'run', 'decs'});
title('Averg time spent in maze');
ylabel('time [s]');

hold off
%% relative time in maye
%mean drink, food, social, explore, run, desicion
drink = [ mean(PROP{1,1}{1,5}) mean(PROP{2,1}{1,5}) mean(PROP{3,1}{1,5}) mean(PROP{4,1}{1,5})];
meand= mean(drink);
food = [ mean(PROP{1,1}{1,1}) mean(PROP{2,1}{1,1}) mean(PROP{3,1}{1,1}) mean(PROP{4,1}{1,1})];
meanf= mean(food);
social = [ mean(PROP{1,1}{1,4}) mean(PROP{2,1}{1,4}) mean(PROP{3,1}{1,4}) mean(PROP{4,1}{1,4})];
means= mean(social);
explore = [ mean(PROP{1,1}{1,3}) mean(PROP{2,1}{1,3}) mean(PROP{3,1}{1,3}) mean(PROP{4,1}{1,3})];
meane= mean(explore);
run = [ mean(PROP{1,1}{1,2}) mean(PROP{2,1}{1,2}) mean(PROP{3,1}{1,2}) mean(PROP{4,1}{1,2})];
meanr= mean(run);
decs = [ mean(PROP{1,1}{1,7}) mean(PROP{2,1}{1,7}) mean(PROP{3,1}{1,7}) mean(PROP{4,1}{1,7})];
meandec= mean(decs);
totalNr = [size(PROP{1,1}{1,6},1) size(PROP{2,1}{1,6},1)/2 size(PROP{3,1}{1,6},1) size(PROP{4,1}{1,6},1)/2];
meanNr = mean(totalNr);

figure % averg time in each pod  (total Nr of blocks in a day, , points are animals
hold on
y = [ meand meanf means meane meanr meandec];
yy = [drink food social explore run decs];
bar(y)
colors = ["r" "g" "b" "y"];
for i = 1:4
    plot([1 2 3 4 5 6 ],[drink(i) food(i) social(i) explore(i) run(i) decs(i)],[colors(i) + '.'],'MarkerSize', 25);
end
xticks([1 2 3 4 5 6 ]);xticklabels({, 'water', 'food', 'social', 'explore', 'run', 'decs'});
title('Relative time spent in maze');
ylabel('Proportion of time');

hold off
%%
%Nr of blocks, duration of blocks
totalNr = [size(SEC{1,1}{1,6},1) size(SEC{2,1}{1,6},1)/2 size(SEC{3,1}{1,6},1) size(SEC{4,1}{1,6},1)/2];
meanNr = mean(totalNr);

totalDur = [mean(SEC{1,1}{1,6},1) mean(SEC{2,1}{1,6},1)/2 mean(SEC{3,1}{1,6},1) mean(SEC{4,1}{1,6},1)/2];
meanDur = mean(totalDur);

figure
hold on
bar(1,meanDur);
for i = 1:4
    plot([1],[totalDur(i)],[colors(i) + '.'],'MarkerSize', 25);
end
ylim([0 710]);
yticks([120 240 360 480 600 720]);yticklabels({'2', '4', '6', '8', '10', '12'});
title('Avergage duration of blocks');
xticks([ ]);xticklabels({''});
ylabel('time [min]');
hold off

figure
hold on
bar(1,meanNr);
for i = 1:4
    plot([1],[totalNr(i)],[colors(i) + '.'],'MarkerSize', 25);
end
%ylim([0 0]);
%yticks([120 240 360 480 600 720]);yticklabels({'2', '4', '6', '8', '10', '12'});
title('Nr of blocks per day');
xticks([ ]);xticklabels({''});
ylabel('Nr of blocks');
hold off


%% Nr of trials per session

%mean drink, food, social, explore, run, desicion entries per session
drink = [ size(SEC{1,1}{1,5},1)/size(SEC{1,1}{1,6},1) size(SEC{2,1}{1,5},1)/size(SEC{2,1}{1,6},1) size(SEC{3,1}{1,5},1)/size(SEC{3,1}{1,6},1) size(SEC{4,1}{1,5},1)/size(SEC{4,1}{1,6},1)];
meand= mean(drink);
food = [ size(SEC{1,1}{1,1},1)/size(SEC{1,1}{1,6},1) size(SEC{2,1}{1,1},1)/size(SEC{2,1}{1,6},1) size(SEC{3,1}{1,1},1)/size(SEC{3,1}{1,6},1) size(SEC{4,1}{1,1},1)/size(SEC{4,1}{1,6},1)];
meanf= mean(food);
social = [ size(SEC{1,1}{1,4},1)/size(SEC{1,1}{1,6},1) size(SEC{2,1}{1,4},1)/size(SEC{2,1}{1,6},1) size(SEC{3,1}{1,4},1)/size(SEC{3,1}{1,6},1) size(SEC{4,1}{1,4},1)/size(SEC{4,1}{1,6},1)];
means= mean(social);
explore = [ size(SEC{1,1}{1,3},1)/size(SEC{1,1}{1,6},1) size(SEC{2,1}{1,3},1)/size(SEC{2,1}{1,6},1) size(SEC{3,1}{1,3},1)/size(SEC{3,1}{1,6},1) size(SEC{4,1}{1,3},1)/size(SEC{4,1}{1,6},1)];
meane= mean(explore);
run = [ size(SEC{1,1}{1,2},1)/size(SEC{1,1}{1,6},1) size(SEC{2,1}{1,2},1)/size(SEC{1,1}{1,6},1) size(SEC{3,1}{1,2},1)/size(SEC{1,1}{1,6},1) size(SEC{4,1}{1,2},1)/size(SEC{1,1}{1,6},1)];
meanr= mean(run);
decs = [ size(SEC{1,1}{1,7},1)/size(SEC{1,1}{1,6},1) size(SEC{2,1}{1,7},1)/size(SEC{1,1}{1,6},1) size(SEC{3,1}{1,7},1)/size(SEC{1,1}{1,6},1) size(SEC{4,1}{1,7},1)/size(SEC{1,1}{1,6},1)];
meandec= mean(decs);

figure 
hold on
y = [ meand meanf means meane meanr ];
yy = [drink food social explore run ];
bar(y)
colors = ["r" "g" "b" "y"];
for i = 1:4
    plot([1 2 3 4 5 ],[drink(i) food(i) social(i) explore(i) run(i)],[colors(i) + '.'],'MarkerSize', 25);
end
xticks([1 2 3 4 5 ]);xticklabels({, 'water', 'food', 'social', 'explore', 'run'});
title('Avergage pod visits per block');
ylabel('Nr of visits')
legend;

hold off

           









