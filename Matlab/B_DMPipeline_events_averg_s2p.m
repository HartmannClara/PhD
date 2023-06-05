clearvars -except ci_data; close all; 
% Analysis of ci_data file, creates events_averg data file
% 03/06/2023
% to do:
%% 
factor = 0.8; % factor by which the neuropil is subtracted 0.7-0.9
skipSessionFrames = 20; %nr of missed frames when sessions are skipped
skippedCounter = 0;FrameSkipped = {}; %counts how many sessions are skipped bc nr of missed frames is too high
skipDrink = 0; skippedDrinkTrials = {};% logs nr and id of drink trials skipped bc of faulty sensor
DrinkNonCons = {};FoodNonCons = {};RunNonCons = {};
DrinkEntryAll = {};FoodEntryAll = {};RunEntryAll = {};
fn=fieldnames(ci_data.bsl);
%loop through the fields
for ii=1: numel(fn)
    fn1=fieldnames(ci_data.bsl.(fn{ii}));
    ed_averg=[];ef_averg=[];es_averg=[];ee_averg=[];er_averg=[];drink_averg=[];eat_averg=[];blockend_averg=[];soc_averg=[];expl_averg =[];
    exd_averg=[];exf_averg=[];exs_averg=[];exe_averg=[];exr_averg=[];imstop_averg=[];run_averg=[];soc_averg_m=[];
    ed_averg_bsl=[];ef_averg_bsl=[];es_averg_bsl=[];ee_averg_bsl=[];er_averg_bsl=[];drink_averg_bsl=[];eat_averg_bsl=[];blockend_averg_bsl=[];soc_averg_bsl=[];expl_averg_bsl =[];
    exd_averg_bsl=[];exf_averg_bsl=[];exs_averg_bsl=[];exe_averg_bsl=[];exr_averg_bsl=[];imstop_averg_bsl=[];run_averg_bsl=[];soc_averg_bsl_m=[];
    PodFrames = [];
    for j=1: numel(fn1)
        fn2=fieldnames(ci_data.bsl.(fn{ii}).(fn1{j}));
        for k=1: numel(fn2)
            for s=1:size(ci_data.bsl.(fn{ii}).(fn1{j}),2)%access the data sessions
                fn3= fieldnames(ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}));
                missed_frames = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{4});
                if numel(missed_frames) >= skipSessionFrames
                    skippedCounter = skippedCounter + 1;
                    fr_skip = {fn(ii), fn1(j), s};
                    FrameSkipped = [FrameSkipped; fr_skip];
                    s = s+1;
                else
                latency =  ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(3); % latency to consume in ms   
                new_event_frames= ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(5);% in events column 5 (frame)
                event_type = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(4);
                raw_f_res=ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
                Fneu = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2});

                subneuro_trace = subneuro(raw_f_res,Fneu,factor); % subtracts neuropil trace from raw trace by factor x subneuro(rawF, Fneuropil, factor)
                subneuro_trace(:,(size(subneuro_trace,2)+1)) = mean(Fneu,2);%last "cell" is mean neuropil trace
                
                df_f_trace = dFoF(subneuro_trace);%df/f
                raw_f= zscore(df_f_trace,0,1);% zscored using sample sd
             
                %%
              
                % check if lick sensor is ok, otherwise skip drink snips
                ed_event_frames = new_event_frames(find(event_type=='enter_drink'));
                d_event_frames = new_event_frames(find(event_type=='drink'));
                ed_d_frames = cat(1,ed_event_frames, d_event_frames);   
                if size(ed_d_frames,1) == size(unique(ed_d_frames),1)
                    pp = 1:15;
                elseif size(ed_d_frames,1) > size(unique(ed_d_frames),1)
                    pp = [1:5 7:15]; skipDrink = skipDrink +1;% skips drink snips
                    skipped = {fn(ii), fn1(j), s};
                    skippedDrinkTrials = [skippedDrinkTrials; skipped];
                end   
                %check for consumption and non consumption enter drink trials
                ed_events = find(event_type=='enter_drink');
                for ev=1:size(ed_events,1)
                    if event_type((ed_events(ev)) + 1) ~= "drink"
                        ed_nc = {fn(ii), fn1(j), s, ev};
                        DrinkNonCons = [DrinkNonCons; ed_nc];
                    end
                end
                 %check for consumption and non consumption enter food trials
                ef_events = find(event_type=='enter_feed');
                for ev=1:size(ef_events,1)
                    if event_type((ef_events(ev)) + 1) ~= "retrieve_pellet"
                        ef_nc = {fn(ii), fn1(j), s, ev};
                        FoodNonCons = [FoodNonCons; ef_nc];
                    end
                end 
                 %check for consumption and non consumption enter run trials
                er_events = find(event_type=='enter_run');
                for ev=1:size(er_events,1)
                    if event_type((er_events(ev)) + 1) ~= "run"
                        er_nc = {fn(ii), fn1(j), s, ev};
                        RunNonCons = [RunNonCons; er_nc];
                    end
                end 

                 %check how long animals stay in pods and how long block is
                 ee_events = find(event_type=='enter_explore'); es_events = find(event_type=='enter_social');
                 exf_events = find(event_type=='exit_feed');exd_events = find(event_type=='exit_drink');exr_events = find(event_type=='exit_run');
                 exs_events = find(event_type=='exit_social');exe_events = find(event_type=='exit_explore');
                 istarts = find(event_type=='imaging_start');B_ends = find(event_type=='block_end');
                 AllEntries = {ef_events,er_events, ee_events,es_events,ed_events,istarts};AllEntriesList = sort(vertcat(AllEntries{:}));
                 AllExits = {exf_events,exr_events, exe_events,exs_events,exd_events,B_ends}; AllExitsList = sort(vertcat(AllExits{:}));

                for et=1:size(AllEntries,2)%event type
                    NrEt = AllEntries{1,et};
                    NrEx = AllExits{1,et};
                    clear PodFrames;
                    for ev = 1:size(NrEt,1) %nr of events
                        PodFrames(ev,1) =  (new_event_frames(NrEx(ev)) - new_event_frames(NrEt(ev)))/10; %frames in the pod divided by framerate
                    end
                    PodTimes{1,et} = PodFrames;
                end
                % legth of time in desicion zone
                AllEntriesList(1:2) = [];% remove imaging start and first entry
                AllEntriesList(end+1,1) =  AllExitsList(end,1);% add block end as last entry into something
                AllExitsList(end)= [];% remove block end from exits
                for et=1:size(AllExitsList,1)%event type
                    DecFrames(et,1) =  (new_event_frames(AllEntriesList(et)) - new_event_frames(AllExitsList(et)))/10;% check frames from exit to next entry
                    DecTypes(et,1) = event_type(AllExitsList(et)); % save what came before and what after
                    DecTypes(et,2) = event_type(AllEntriesList(et));
                end
                Desicion = {DecFrames, DecTypes};
                TimeSpecs = [TimeSpecs; PodTimes];%pod times: food, run, drink, explore, social, whole block
                DecSpecs = [DecSpecs; Desicion];% Decision zone: time, exit to entry types
                
                for p=pp
                    clear snips snips_bslined snips2 snips_bslined2
                    if p==1
                        t='enter_drink';
                        tt='enter drink';
                    elseif p==2
                        t='enter_feed';
                        tt='enter feed';
                    elseif p==3
                        t='enter_social';
                        tt='enter social';
                    elseif p==4
                        t='enter_explore';
                        tt='enter explore';
                    elseif p==5
                        t='enter_run';
                        tt='enter run';
                    elseif p==6
                        t='drink';
                        tt='drink';
                    elseif p==7
                        t='retrieve_pellet';   
                        tt='retrieve pellet';
                    elseif p==8
                        t='block_end';   
                        tt='block_end';      
                    elseif p==9
                        t='exit_drink';
                        tt='exit drink';
                    elseif p==10
                        t='exit_feed';
                        tt='exit feed';
                    elseif p==11
                        t='exit_social';
                        tt='exit social';
                    elseif p==12
                        t='exit_explore';
                        tt='exit explore';
                    elseif p==13
                        t='exit_run';
                        tt='exit run';
                    elseif p==14
                        t='imaging_stop';
                        tt='imaging stop';
                    elseif p==15
                        t='run';
                        tt='run';    
                    end
                    
                    trigs=new_event_frames(find(event_type==t));%trigs are the indexes of events of event type t
                    if size(trigs,1)==0 % if that event does not exist in the session, move on
                        p=p+1;
                    else
                        win = 30;% event window
                        fte = 25; %frames to entry/exit; how long it takes an animal to go from door to entry beam

                        for i=1:size(trigs,1)%loop through the events of type t
                            trig=trigs(i);
                            
                                if trig>win+1 & p==1 
                                    snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); %entry drink, entries take 3s/30 frames bf door entry (fte before entry beam break)
                                    entry_d_bsl = mean(snips(1:10,:,i),1);% baselined to 30-20s before door entry
                                    snips_bslined(:,:,i)=snips(:,:,i)-entry_d_bsl;
                                    ed_averg=cat(3,ed_averg,snips(:,:,i));
                                    ed_averg_bsl=cat(3,ed_averg_bsl,snips_bslined(:,:,i));
                                    ed_all = {fn(ii), fn1(j), s, i};
                                    DrinkEntryAll = [DrinkEntryAll; ed_all];
                                elseif trig>win+1 & p==2
                                    snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); % entry food
                                    entry_f_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
                                    ef_averg=cat(3,ef_averg,snips(:,:,i));
                                    ef_averg_bsl=cat(3,ef_averg_bsl,snips_bslined(:,:,i));
                                    ef_all = {fn(ii), fn1(j), s, i};
                                    FoodEntryAll = [FoodEntryAll; ef_all];
                                    %% consumptions without beam break
                                elseif trig>win+1 & p==3
                                    snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); % entry social
                                    entry_s_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-entry_s_bsl;
                                    es_averg=cat(3,es_averg,snips(:,:,i));
                                    es_averg_bsl=cat(3,es_averg_bsl,snips_bslined(:,:,i));
                    
                                    snips2(:,:,i)=raw_f(trig:trig+win,:); %social interact, takes 3s after beam break 
                                    entry_soc_bsl = mean(snips(20:30,:,i),1);%baselined to 1 sec bf entry beam break
                                    snips_bslined2(:,:,i)=snips2(:,:,i)-entry_soc_bsl;
                                    soc_averg=cat(3,soc_averg,snips2(:,:,i));
                                    soc_averg_bsl=cat(3,soc_averg_bsl,snips_bslined2(:,:,i));
                    
                                elseif trig>win+1 & p==4
                                    snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); %entry explore
                                    entry_e_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-entry_e_bsl;
                                    ee_averg=cat(3,ee_averg,snips(:,:,i));
                                    ee_averg_bsl=cat(3,ee_averg_bsl,snips_bslined(:,:,i));

                                    snips2(:,:,i)=raw_f(trig:trig+win,:); %explore activity, takes 3s after beam break
                                    entry_expl_bsl = mean(snips(20:30,:,i),1);%baselined to 1 sec bf entry beam break
                                    snips_bslined2(:,:,i)=snips2(:,:,i)-entry_expl_bsl;
                                    expl_averg=cat(3,expl_averg,snips2(:,:,i));
                                    expl_averg_bsl=cat(3,expl_averg_bsl,snips_bslined2(:,:,i));

                                elseif trig>win+1 & p==5
                                    snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); %entry run
                                    entry_r_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-entry_r_bsl;
                                    er_averg=cat(3,er_averg,snips(:,:,i));
                                    er_averg_bsl=cat(3,er_averg_bsl,snips_bslined(:,:,i));
                                    er_all = {fn(ii), fn1(j), s, i};
                                    RunEntryAll = [RunEntryAll; er_all];
                                % consumptions with beam breaks    
                                elseif trig>win+1 & p==6 %drink, takes 3s before and 3s after drink 
                                    snips(:,:,i)=raw_f(trig-win:trig+win,:);
                                    d_bsl = mean(snips(1:10,:,i),1); %baselined to first first 10 frames of 60 frame event window (2s before drinking)
                                    snips_bslined(:,:,i)=snips(:,:,i)-d_bsl;
                                    drink_averg=cat(3,drink_averg,snips(:,:,i));
                                    drink_averg_bsl=cat(3,drink_averg_bsl,snips_bslined(:,:,i));
                                elseif trig>win+1 & p==7 %pellet, takes 3s before and 3s after beam break 
                                    snips(:,:,i)=raw_f(trig-win:trig+win,:);
                                    f_bsl = mean(snips(1:10,:,i),1);%takes 1s right after entry before pellet
                                    snips_bslined(:,:,i)=snips(:,:,i)-f_bsl;
                                    eat_averg=cat(3,eat_averg,snips(:,:,i));
                                    eat_averg_bsl=cat(3,eat_averg_bsl,snips_bslined(:,:,i));
                                elseif trig>win+1 & p==8 %blockend
                                    snips(:,:,i)=raw_f(trig-win:trig+200,:);%blockend  from -30 to 200 frames bf/after blockend 23s, be at 31frames
                                    be_bsl = mean(snips(1:10,:,i),1);%baselined to 2 seconds before be
                                    snips_bslined(:,:,i)=snips(:,:,i)-be_bsl;
                                    blockend_averg=cat(3,blockend_averg,snips(:,:,i));
                                    blockend_averg_bsl=cat(3,blockend_averg_bsl,snips_bslined(:,:,i));
%% exits
                                   
                                elseif trig>win+1 & p==9%exit drink
                                    snips(:,:,i)=raw_f(trig-win:trig+win,:); %maybe its more interesting what happens after the door is crossed (30s before beam ) and when animal is in center(+win)
                                    exd_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-exd_bsl;
                                    exd_averg=cat(3,exd_averg,snips(:,:,i));
                                    exd_averg_bsl=cat(3,exd_averg_bsl,snips_bslined(:,:,i));
                                elseif trig>win+1 & p==10
                                    snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit food
                                    exf_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-exf_bsl;
                                    exf_averg=cat(3,exf_averg,snips(:,:,i));
                                    exf_averg_bsl=cat(3,exf_averg_bsl,snips_bslined(:,:,i));
                                elseif trig>win+1 & p==11
                                    snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit social
                                    exs_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-exs_bsl;
                                    exs_averg=cat(3,exs_averg,snips(:,:,i));
                                    exs_averg_bsl=cat(3,exs_averg_bsl,snips_bslined(:,:,i));
                                elseif trig>win+1 & p==12
                                    snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit explore
                                    exe_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-exe_bsl;
                                    exe_averg=cat(3,exe_averg,snips(:,:,i));
                                    exe_averg_bsl=cat(3,exe_averg_bsl,snips_bslined(:,:,i));
                                elseif trig>win+1 & p==13
                                    snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit run
                                    exr_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-exr_bsl;
                                    exr_averg=cat(3,exr_averg,snips(:,:,i));
                                    exr_averg_bsl=cat(3,exr_averg_bsl,snips_bslined(:,:,i));
                                elseif trig>win+1 & p==14 & trig<= size(raw_f,1)%imaging stop 
                                    snips(:,:,i)=raw_f(trig-win:trig,:);
                                    snips_bslined(:,:,i)=snips(:,:,i)-be_bsl;
                                    imstop_averg=cat(3,imstop_averg,snips(:,:,i));
                                    imstop_averg_bsl=cat(3,imstop_averg_bsl,snips_bslined(:,:,i));
                                elseif trig>win+1 & p==15 %run
                                    snips(:,:,i)=raw_f(trig-win:trig+win,:); % run
                                    r_bsl = mean(snips(1:10,:,i),1);%bsl 2s before run
                                    snips_bslined(:,:,i)=snips(:,:,i)-r_bsl;
                                    run_averg=cat(3,run_averg,snips(:,:,i));
                                    run_averg_bsl=cat(3,run_averg_bsl,snips_bslined(:,:,i));
                                           
                                end                        
                            end
                        end  
                end
                end
            end    
        end 
    end
    %%saving into events_averg and Specs structs
    Specs.(fn{ii}).PodTimes = TimeSpecs;
    Specs.(fn{ii}).DecTimes = DecSpecs;
    events_averg.(fn{ii}).ed_averg = ed_averg;
    events_averg.(fn{ii}).ef_averg = ef_averg;
    events_averg.(fn{ii}).es_averg = es_averg;
    events_averg.(fn{ii}).soc_averg = soc_averg;
    events_averg.(fn{ii}).ee_averg = ee_averg;
    events_averg.(fn{ii}).er_averg = er_averg;
    events_averg.(fn{ii}).drink_averg = drink_averg;
    events_averg.(fn{ii}).eat_averg = eat_averg;
    events_averg.(fn{ii}).blockend_averg = blockend_averg;
    events_averg.(fn{ii}).exd_averg = exd_averg;
    events_averg.(fn{ii}).exf_averg = exf_averg;
    events_averg.(fn{ii}).exs_averg = exs_averg;
    events_averg.(fn{ii}).exe_averg = exe_averg;
    events_averg.(fn{ii}).exr_averg = exr_averg;
    events_averg.(fn{ii}).run_averg = run_averg;
    events_averg.(fn{ii}).expl_averg = expl_averg;

    events_averg.(fn{ii}).ed_averg_bsl = ed_averg_bsl;
    events_averg.(fn{ii}).ef_averg_bsl = ef_averg_bsl;
    events_averg.(fn{ii}).es_averg_bsl = es_averg_bsl;
    events_averg.(fn{ii}).soc_averg_bsl = soc_averg_bsl;
    events_averg.(fn{ii}).ee_averg_bsl = ee_averg_bsl;
    events_averg.(fn{ii}).er_averg_bsl = er_averg_bsl;
    events_averg.(fn{ii}).drink_averg_bsl = drink_averg_bsl;
    events_averg.(fn{ii}).eat_averg_bsl = eat_averg_bsl;
    events_averg.(fn{ii}).blockend_averg_bsl = blockend_averg_bsl;
    events_averg.(fn{ii}).exd_averg_bsl = exd_averg_bsl;
    events_averg.(fn{ii}).exf_averg_bsl = exf_averg_bsl;
    events_averg.(fn{ii}).exs_averg_bsl = exs_averg_bsl;
    events_averg.(fn{ii}).exe_averg_bsl = exe_averg_bsl;
    events_averg.(fn{ii}).exr_averg_bsl = exr_averg_bsl;
    events_averg.(fn{ii}).run_averg_bsl = run_averg_bsl;
    events_averg.(fn{ii}).expl_averg_bsl = expl_averg_bsl; 
end


%% Cons and non consumption trials
DrinkEntryAll = cell2table(DrinkEntryAll, "VariableNames",["animal" "day" "session" "trial"]); FoodEntryAll = cell2table(FoodEntryAll, "VariableNames",["animal" "day" "session" "trial"]); RunEntryAll = cell2table(RunEntryAll, "VariableNames",["animal" "day" "session" "trial"]);
DrinkNonCons = cell2table(DrinkNonCons, "VariableNames",["animal" "day" "session" "trial"]); FoodNonCons = cell2table(FoodNonCons, "VariableNames",["animal" "day" "session" "trial"]); RunNonCons = cell2table(RunNonCons, "VariableNames",["animal" "day" "session" "trial"]);
% saves index (reffering to list of all entry trials in order) of non
% consummatory trials for water, food and run pod
%events_averg.non_cons.D_all = DrinkEntryAll;events_averg.non_cons.F_all = FoodEntryAll;events_averg.non_cons.R_all = RunEntryAll;

% split into animals 
G = findgroups(FoodEntryAll{:,1});
T_split = splitapply( @(varargin) varargin, FoodEntryAll , G);
subTables = cell(size(T_split, 1));% Allocate empty cell array fo sizxe equal to number of rows in T_Split
for i = 1:size(T_split, 1)% Create sub tables
subTables{i} = table(T_split{i, :}, 'VariableNames', ...
FoodEntryAll.Properties.VariableNames);
end
[Lia,efNCg2] = ismember(FoodNonCons,subTables{1,1});[Lia,efNCg4] = ismember(FoodNonCons,subTables{2,1});[Lia,efNCg5] = ismember(FoodNonCons,subTables{3,1});
efNCg2 = nonzeros(efNCg2);efNCg4 = nonzeros(efNCg4);efNCg5 = nonzeros(efNCg5); % these are the indexs of non consumption trials in events_averg. "animal".ef_averg(ef_averg_bsl)

G = findgroups(DrinkEntryAll{:,1});
T_split = splitapply( @(varargin) varargin, DrinkEntryAll , G);
subTables = cell(size(T_split, 1));% Allocate empty cell array fo sizxe equal to number of rows in T_Split
for i = 1:size(T_split, 1)% Create sub tables
subTables{i} = table(T_split{i, :}, 'VariableNames', ...
DrinkEntryAll.Properties.VariableNames);
end
[Lia,edNCg2] = ismember(DrinkNonCons,subTables{1,1});[Lia,edNCg4] = ismember(DrinkNonCons,subTables{2,1});[Lia,edNCg5] = ismember(DrinkNonCons,subTables{3,1});
edNCg2 = nonzeros(edNCg2);edNCg4 = nonzeros(edNCg4);edNCg5 = nonzeros(edNCg5); % these are the indexs of non consumption trials in events_averg. "animal".ed_averg(ed_averg_bsl)


G = findgroups(RunEntryAll{:,1});
T_split = splitapply( @(varargin) varargin, RunEntryAll , G);
subTables = cell(size(T_split, 1));% Allocate empty cell array fo sizxe equal to number of rows in T_Split
for i = 1:size(T_split, 1)% Create sub tables
subTables{i} = table(T_split{i, :}, 'VariableNames', ...
RunEntryAll.Properties.VariableNames);
end
[Lia,erNCg2] = ismember(RunNonCons,subTables{1,1});[Lia,erNCg4] = ismember(RunNonCons,subTables{2,1});[Lia,erNCg5] = ismember(RunNonCons,subTables{3,1});
erNCg2 = nonzeros(erNCg2);erNCg4 = nonzeros(erNCg4);erNCg5 = nonzeros(erNCg5); % these are the indexes of non consumption trials in events_averg. "animal".er_averg(er_averg_bsl)


%split in trials in events_averg
%enter drink
events_averg.g2.ed_nc_averg_bsl = events_averg.g2.ed_averg_bsl(:,:,edNCg2);% non consumption trials
consum = events_averg.g2.ed_averg_bsl; consum(:,:, edNCg2) =[];events_averg.g2.ed_c_averg_bsl= consum;%consumption trials
events_averg.g5.ed_nc_averg_bsl = events_averg.g5.ed_averg_bsl(:,:,edNCg5);% non consumption trials
consum = events_averg.g5.ed_averg_bsl; consum(:,:, edNCg5) =[];events_averg.g5.ed_c_averg_bsl = consum;%consumption trials
events_averg.g4.ed_nc_averg_bsl = events_averg.g4.ed_averg_bsl(:,:,edNCg4);% non consumption trials
consum = events_averg.g4.ed_averg_bsl; consum(:,:, edNCg4) =[];events_averg.g4.ed_c_averg_bsl = consum;%consumption trials

%enter food
events_averg.g2.ef_nc_averg_bsl = events_averg.g2.ef_averg_bsl(:,:,efNCg2);% non consumption trials
consum = events_averg.g2.ef_averg_bsl; consum(:,:, efNCg2) =[];events_averg.g2.ef_c_averg_bsl = consum;%consumption trials
events_averg.g5.ef_nc_averg_bsl = events_averg.g5.ef_averg_bsl(:,:,efNCg5);% non consumption trials
consum = events_averg.g5.ef_averg_bsl; consum(:,:, efNCg5) =[];events_averg.g5.ef_c_averg_bsl = consum;%consumption trials
events_averg.g4.ef_nc_averg_bsl = events_averg.g4.ef_averg_bsl(:,:,efNCg4);% non consumption trials
consum = events_averg.g4.ef_averg_bsl; consum(:,:, efNCg4) =[];events_averg.g4.ef_c_averg_bsl = consum;%consumption trials

%enter run

save("events_averg.mat", "events_averg");

