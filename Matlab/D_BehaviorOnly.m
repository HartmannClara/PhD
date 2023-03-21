% analysis including behavior
%trials separated in runs (2+) or singles
clearvars -except ci_data entriesBehavior; close all; % Analysis of ci_data file, creates events_averg data file
% 06.03.2023
%% correcting out extra data (imaging data) 
% fn=fieldnames(ci_data_beh.bsl);
% for ii =3:5
%     fn1=fieldnames(ci_data_beh.bsl.(fn{ii}));
%     for j=1: numel(fn1)%day
%        fn2=fieldnames(ci_data_beh.bsl.(fn{ii}).(fn1{j}));
%        for k=1: numel(fn2)%session
%           for s=1:size(ci_data_beh.bsl.(fn{ii}).(fn1{j}),2)%access the data
%               fn3= fieldnames(ci_data_beh.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}));  
%               eve = ci_data_beh.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2});
%               ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).events=eve; 
%           end
%        end
%     end
% end
% ci_data_beh.bsl.g5 = ci_data.bsl.g5;
% ci_data_beh.bsl.g2 = ci_data.bsl.g2;
% ci_data_beh.bsl.g4 = ci_data.bsl.g4;



ci_data2 = ci_data;

%%

skipSessionFrames = 20; %nr of missed frames when sessions are skipped
skippedCounter = 0;FrameSkipped = {}; %counts how many sessions are skipped bc nr of missed frames is tooSwap hi C_matrix_averggh C_matrix_averg
skipDrink = 0; skippedDrinkTrials = {};% logs nr and id of drink trials skipped bc of faulty sensor
DrinkNonCons = {};FoodNonCons = {};RunNonCons = {};
%DrinkEntryAll = {};FoodEntryAll = {};RunEntryAll = {};SocEntryAll = {};ExpEntryAll={};
Entry_clsf = {};
nr_ef = {}; confus = zeros(6,6,6);nr_entries = zeros(1,6,6);
fn=fieldnames(ci_data.bsl);
%loop through the fields
for ii=1: numel(fn)%animal
    fn1=fieldnames(ci_data.bsl.(fn{ii}));
    Entry_clsf = {};
    DrinkEntryAll = {};FoodEntryAll = {};RunEntryAll = {};SocEntryAll = {};ExpEntryAll={};
    for j=1: numel(fn1)%day
        fn2=fieldnames(ci_data.bsl.(fn{ii}).(fn1{j}));
        for k=1: numel(fn2)%session
            for s=1:size(ci_data.bsl.(fn{ii}).(fn1{j}),2)%access the data
                fn3= fieldnames(ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}));  
                %missed_frames = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3});
                event_type = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1}).(4);
                hw_time = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1}).(6);  
                %win =30; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % check if lick sensor is ok, otherwise skip drink snips
                %check for consumption and non consumption enter drink trials
                ed_events = find(event_type=='enter_drink');
                ef_events = find(event_type=='enter_feed');
                er_events = find(event_type=='enter_run');
                ee_events = find(event_type=='enter_explore');
                es_events = find(event_type=='enter_social');
                be_events = find(event_type=='block_end');
                entries_total=sort([ef_events;er_events; ee_events; es_events; ed_events]);
                entries_total_be=sort([ef_events;er_events; ee_events; es_events; ed_events;be_events]);

                if numel(entries_total) > 0  && event_type((entries_total(1))+1) == "run" |  event_type((entries_total(1))+1) == "retrieve_pellet"| event_type((entries_total(1))+1) == "drink"
                    consum = true;
                else
                    consum = false;
                end
    
                count =1;
                for ev=1:size(entries_total,1)
                    if ev == size(entries_total,1) && count > 1 % before is same
                       clsf = "er" ;
                       if event_type((entries_total(ev))+1) == "run" |  event_type((entries_total(ev))+1) == "retrieve_pellet"| event_type((entries_total(ev))+1) == "drink"
                           consum = true;
                       else
                           consum = false;
                       end
                       e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(ev))) ,char(event_type((entries_total(ev)))), clsf, count, consum};
                       Entry_clsf = [Entry_clsf; e_clsf];
                       count = 1; 
                    elseif ev == size(entries_total,1) && count == 1 % before is different, last can only be a single
                       clsf = "s" ; ;%single
                       if event_type((entries_total(ev))+1) == "run" |  event_type((entries_total(ev))+1) == "retrieve_pellet"| event_type((entries_total(ev))+1) == "drink"
                           consum = true;
                       else
                           consum = false;
                       end
                       e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(ev))) ,char(event_type((entries_total(ev)))), clsf, (count-1), consum};
                       Entry_clsf = [Entry_clsf; e_clsf];   

                    elseif count == 1 && event_type((entries_total(ev))) ~= event_type((entries_total(ev+1))) %if next one is different its a single
                       clsf = "s" ; ;%single
                       if event_type((entries_total(ev))+1) == "run" |  event_type((entries_total(ev))+1) == "retrieve_pellet"| event_type((entries_total(ev))+1) == "drink"
                           consum = true;
                       else
                           consum = false;
                       end
                       e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(ev))) ,char(event_type((entries_total(ev)))), clsf, (count-1), consum};
                       Entry_clsf = [Entry_clsf; e_clsf];

                    elseif count == 1 && event_type((entries_total(ev))) == event_type((entries_total(ev+1))) %if next is the same but before is diff(count==1) its start of run
                       clsf = "sr" ; %1st of run 
                       e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(ev))) ,char(event_type((entries_total(ev)))), clsf, count, consum};
                       Entry_clsf = [Entry_clsf; e_clsf];
                       count = 2;
                       
                    elseif count > 1 && event_type((entries_total(ev))) == event_type((entries_total(ev+1)))% middle of a run %if the same as before
                       clsf ="r" ;
                       if event_type((entries_total(ev))+1) == "run" |  event_type((entries_total(ev))+1) == "retrieve_pellet"| event_type((entries_total(ev))+1) == "drink"
                           consum = true;
                       else
                           consum = false;
                       end
                       e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(ev))) ,char(event_type((entries_total(ev)))), clsf, count, consum};
                       Entry_clsf = [Entry_clsf; e_clsf];
                       count = count + 1;

                    elseif count > 1 && event_type((entries_total(ev))) ~= event_type((entries_total(ev+1)))%end of a run
                       clsf = "er";
                       if event_type((entries_total(ev))+1) == "run" |  event_type((entries_total(ev))+1) == "retrieve_pellet"| event_type((entries_total(ev))+1) == "drink"
                           consum = true;
                       else
                           consum = false;
                       end
                       e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(ev))) ,char(event_type((entries_total(ev)))), clsf, count, consum};
                       Entry_clsf = [Entry_clsf; e_clsf];
                       count = 1;                     
                    end
                end

                entry_types = ["enter_feed", "enter_run","enter_explore","enter_social","enter_drink", "block_end"];           
                
                for trial =1:(size(entry_types,2)-1)%without block end bc nothing comes after
                    trial_type = entry_types(trial);
                    for ent =1:size(entries_total_be,1)
                        if event_type((entries_total_be(ent))) == trial_type
                            nr_entries(1,trial,ii) = nr_entries(1,trial,ii)+1; % collects the total number of entries of that type
                            for et =1:size(entry_types,2)
                                if event_type((entries_total_be(ent+1))) == entry_types(et) %checks which is the next entry and collects it
                                    confus(trial,et,ii)=confus(trial,et,ii) + 1;
                                end
                            end
                        end
                    end 
                end
                %transitions collected per animal
                nr_entries(1,6,ii) = nr_entries(1,6,ii)+1; %collects nr of block_ends = nr of blocks 
                %collect first entry type
                for trial =1:(size(entry_types,2)-1)
                    trial_type = entry_types(trial);
                    if event_type((entries_total_be(1))) == trial_type
                       confus(6,trial,ii) = confus(6,trial,ii)+1;
                    end
                end 
            end   
        end 
    end
    % ensure that number of trials align
    %Entry_all = cell2table(cat(1,FoodEntryAll, RunEntryAll, ExpEntryAll, SocEntryAll, DrinkEntryAll),"VariableNames",["animal" "day" "session" "nr" "hw"]);
    %Entry_all = Entry_all(:,[1 2 3 5]);
    EntryClsf_2 = cell2table(Entry_clsf,"VariableNames",["animal" "day" "session" "hw" "type" "clasf" "count" "consum"]);%EntryCl = EntryClsf(:,1:4);
    %Lia = ismember(EntryCl,Entry_all);
    %EntryClsf_2 = EntryClsf(Lia,:);
    % saves in struct according to animal
    types = ["ef", "er","ee","es","ed"];  
    for ty = 1:size(types,2)
         animal =char(fn{ii});type = char(types{ty});
         Clsf_entry.(animal).(type) = EntryClsf_2(EntryClsf_2.type == entry_types(ty), 6:8);
    end
end

%%
 
% code to count singles and different type of runs
fn = fieldnames(Clsf_entry);
for ii=1: numel(fn)%animal
    fn1=fieldnames(Clsf_entry.(fn{ii}));
    for jj = 1:numel(fn1) %entry type
        SRevents(1,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "s" ),1);
        SRevents(2,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "er" & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,2)) == 2),1);
        SRevents(3,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "er" & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,2)) == 3),1);
        SRevents(4,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "er" & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,2)) == 4),1);
        SRevents(5,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "er" & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,2)) >= 5),1);% five or more
        SR_classif.(fn{ii}).(fn1{jj}) = SRevents;
    end
end   
fn = fieldnames(SR_classif.g5);
for ii=1: numel(fn)%entry type
    SRall = [];
    fn1=fieldnames(Clsf_entry);
    for jj = 1:numel(fn1)%animal 
       SReventsAll = SR_classif.(fn1{jj}).(fn{ii});
       SRall = [SRall, SReventsAll];
       SR_classif.all.(fn{ii}) = SRall; 
    end
end 
%% code to count consumption and non consumption
fn = fieldnames(Clsf_entry);
for ii=1: numel(fn)%animal
    fn1=fieldnames(Clsf_entry.(fn{ii}));
    for jj = 1:numel(fn1) %entry type
        CNCevents(1,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "s" & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,3)) == 1),1); %consumption singles
        CNCevents(2,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "s" & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,3)) == 0),1); %non c singles
        CNCevents(3,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "sr" & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,3)) == 1),1);%c 1st run
        CNCevents(4,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "sr" & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,3)) == 0),1);%nc 1st run
        CNCevents(5,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "r"  & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,3)) == 1),1);%c run
        CNCevents(6,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "r"  & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,3)) == 0),1);%nc  run 
        CNCevents(7,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "er"  & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,3)) == 1),1);%c end run
        CNCevents(8,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "er"  & table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,3)) == 0),1);%nc end run
        CNC_classif.(fn{ii}).(fn1{jj}) = CNCevents;
    end
end   
fn = fieldnames(CNC_classif.g5);
for ii=1: numel(fn)%entry type
    CNCall = [];
    fn1=fieldnames(Clsf_entry);
    for jj = 1:numel(fn1)%animal 
       CNCeventsAll = CNC_classif.(fn1{jj}).(fn{ii});
       CNCall = [CNCall, CNCeventsAll];
       CNC_classif.all.(fn{ii}) = CNCall; 
    end
end 

%% calculate correct confusion matrix values
C_matrix = [];
for a = 1:size(confus,3)
    for r = 1:6
        for c = 1:6
            C_matrix(r,c,a) = confus(r,c,a)/nr_entries(1,r,a);
        end  
    end
end
C_matrix(isnan(C_matrix))=0;

C_matrix_averg = [];C_matrix_av = 0;
for r = 1:6
       for c = 1:6
           for a=1:size(C_matrix,3)
                C_matrix_av =  C_matrix_av + C_matrix(r,c,a);
           end
           C_matrix_averg (r,c) = C_matrix_av/size(C_matrix,3);
           C_matrix_av = 0;
       end  
end

T_Cmatrix_av = array2table(C_matrix_averg,"RowNames",["ef" "er" "ee" "es" "ed" "h"], "VariableNames",["ef" "er" "ee" "es" "ed" "h"]);
entriesBehavior.bsl.nr_entries = nr_entries;
entriesBehavior.bsl.counts = confus;
entriesBehavior.bsl.matrix = C_matrix;
entriesBehavior.bsl.matrix_averg = C_matrix_averg;
entriesBehavior.bsl.matrix_tbl = T_Cmatrix_av;
entriesBehavior.bsl.SR_classif = SR_classif;
entriesBehavior.bsl.CNC_classif = CNC_classif;

%save("SR_classif.mat", "SR_classif");
save("entriesBehavior2.mat", "entriesBehavior");
%save("events_classif.mat", "events_classif");
%save("events_averg_120323.mat", "events_averg");

