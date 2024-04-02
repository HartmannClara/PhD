clearvars -except ci_data events_averg_unshuf; close all; 
%%Shuffled data 
% Analysis of ci_data file, creates events_averg data file
% 10.01.2023
%% 
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
    for j=1: numel(fn1)
        fn2=fieldnames(ci_data.bsl.(fn{ii}).(fn1{j}));
        for k=1: numel(fn2)
            for s=1:size(ci_data.bsl.(fn{ii}).(fn1{j}),2)%access the data
                fn3= fieldnames(ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}));
                missed_frames = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3});
                if numel(missed_frames) >= skipSessionFrames
                    skippedCounter = skippedCounter + 1;
                    fr_skip = {fn(ii), fn1(j), s};
                    FrameSkipped = [FrameSkipped; fr_skip];
                    s = s+1;
                else
                new_event_frames= ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2}).(5);
                event_type = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2}).(4);
                raw_f_res=ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
                df_f_trace = dFoF(raw_f_res);%df/f
                raw_f= zscore(df_f_trace,0,1);% zscored using sample sd
                %smoothing?
                %raw_f = movmean(raw_f,3,1);% 3 point centered moving average
                raw_f_z =zscore(df_f_trace,0,1); 
                %%shuffle raw_f
                shuff = [];
                shufFac = randi([100 1000],1,200);
                for f = 1:(numel(shufFac))
                    shuff(:,:,f)= circshift(raw_f,(round((shufFac(f)*size(raw_f,1))/1000)),1);%shifts raw_f by random shufFac
                end
                %raw_f = nanmean(shuff,3);% averages over the iterations to make shuffled trace

                ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).raw_f_shuf = shuff;
                ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).raw_f = raw_f_z;
                
%                 ses = 6;cell = 6;
%                 figure; plot(ci_data.bsl.g5.d_2022_09_08(ses).session.raw_f(:,cell));hold on; 
%                 plot(ci_data.bsl.g5.d_2022_09_08(ses).session.raw_f_shuf(:,cell));
%                 
                %%
                for sh = 1:size(shuff,3)
                    raw_f = shuff(:,:,sh);
                    win =30;
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
                            for i=1:size(trigs,1)%loop through the events of type t
                                trig=trigs(i);
                                
                                    if trig>win+1 & p==1 
                                        snips(:,:,i)=raw_f(trig-win:trig,:); %entry drink, entries take 3s/30frames bf beam break
                                        entry_d_bsl = mean(snips(1:10,:,i),1);
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_d_bsl;
                                        ed_averg=cat(3,ed_averg,snips(:,:,i));
                                        ed_averg_bsl=cat(3,ed_averg_bsl,snips_bslined(:,:,i));
                                        ed_all = {fn(ii), fn1(j), s, i};
                                        DrinkEntryAll = [DrinkEntryAll; ed_all];
                                    elseif trig>win+1 & p==2
                                        snips(:,:,i)=raw_f(trig-win:trig,:); % entry food
                                        entry_f_bsl = mean(snips(1:10,:,i),1);
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
                                        ef_averg=cat(3,ef_averg,snips(:,:,i));
                                        ef_averg_bsl=cat(3,ef_averg_bsl,snips_bslined(:,:,i));
                                        ef_all = {fn(ii), fn1(j), s, i};
                                        FoodEntryAll = [FoodEntryAll; ef_all];
                                    elseif trig>win+1 & p==3
                                        snips(:,:,i)=raw_f(trig-win:trig,:); % entry social
                                        entry_s_bsl = mean(snips(1:10,:,i),1);
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_s_bsl;
                                        es_averg=cat(3,es_averg,snips(:,:,i));
                                        es_averg_bsl=cat(3,es_averg_bsl,snips_bslined(:,:,i));
                        
                                        snips2(:,:,i)=raw_f(trig:trig+win,:); %social interact, takes 3s after beam break
                                        snips_bslined2(:,:,i)=snips2(:,:,i)-entry_s_bsl;
                                        soc_averg=cat(3,soc_averg,snips2(:,:,i));
                                        soc_averg_bsl=cat(3,soc_averg_bsl,snips_bslined2(:,:,i));
                        
                                    elseif trig>win+1 & p==4
                                        snips(:,:,i)=raw_f(trig-win:trig,:); %entry explore
                                        entry_e_bsl = mean(snips(1:10,:,i),1);
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_e_bsl;
                                        ee_averg=cat(3,ee_averg,snips(:,:,i));
                                        ee_averg_bsl=cat(3,ee_averg_bsl,snips_bslined(:,:,i));
    
                                        snips2(:,:,i)=raw_f(trig:trig+win,:); %explore activity, takes 3s after beam break
                                        snips_bslined2(:,:,i)=snips2(:,:,i)-entry_e_bsl;
                                        expl_averg=cat(3,expl_averg,snips2(:,:,i));
                                        expl_averg_bsl=cat(3,expl_averg_bsl,snips_bslined2(:,:,i));
    
                                    elseif trig>win+1 & p==5
                                        snips(:,:,i)=raw_f(trig-win:trig,:); %entry run
                                        entry_r_bsl = mean(snips(1:10,:,i),1);
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_r_bsl;
                                        er_averg=cat(3,er_averg,snips(:,:,i));
                                        er_averg_bsl=cat(3,er_averg_bsl,snips_bslined(:,:,i));
                                        er_all = {fn(ii), fn1(j), s, i};
                                        RunEntryAll = [RunEntryAll; er_all];
                                    elseif trig>win+1 & p==6 %drink, takes 3s before and 3s after drink 
                                        snips(:,:,i)=raw_f(trig-win:trig+win,:);
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_d_bsl;
                                        drink_averg=cat(3,drink_averg,snips(:,:,i));
                                        drink_averg_bsl=cat(3,drink_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==7 %pellet, takes 3s before and 3s after beam break 
                                        snips(:,:,i)=raw_f(trig-win:trig+win,:);
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
                                        eat_averg=cat(3,eat_averg,snips(:,:,i));
                                        eat_averg_bsl=cat(3,eat_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==8 %blockend
                                        snips(:,:,i)=raw_f(trig-win:trig+100,:);%blockend and 10 seconds
                                        be_bsl = mean(snips(1:10,:,i),1);
                                        snips_bslined(:,:,i)=snips(:,:,i)-be_bsl;
                                        blockend_averg=cat(3,blockend_averg,snips(:,:,i));
                                        blockend_averg_bsl=cat(3,blockend_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==9
                                        snips(:,:,i)=raw_f(trig-win:trig+win,:); % exit drink
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_d_bsl;
                                        exd_averg=cat(3,exd_averg,snips(:,:,i));
                                        exd_averg_bsl=cat(3,exd_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==10
                                        snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit food
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
                                        exf_averg=cat(3,exf_averg,snips(:,:,i));
                                        exf_averg_bsl=cat(3,exf_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==11
                                        snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit social
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_s_bsl;
                                        exs_averg=cat(3,exs_averg,snips(:,:,i));
                                        exs_averg_bsl=cat(3,exs_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==12
                                        snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit explore
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_e_bsl;
                                        exe_averg=cat(3,exe_averg,snips(:,:,i));
                                        exe_averg_bsl=cat(3,exe_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==13
                                        snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit run
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_r_bsl;
                                        exr_averg=cat(3,exr_averg,snips(:,:,i));
                                        exr_averg_bsl=cat(3,exr_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==14 & trig<= size(raw_f,1)%imaging stop 
                                        snips(:,:,i)=raw_f(trig-win:trig,:);
                                        snips_bslined(:,:,i)=snips(:,:,i)-be_bsl;
                                        imstop_averg=cat(3,imstop_averg,snips(:,:,i));
                                        imstop_averg_bsl=cat(3,imstop_averg_bsl,snips_bslined(:,:,i));
                                    elseif trig>win+1 & p==15 %run
                                        snips(:,:,i)=raw_f(trig-win:trig+win,:); % run
                                        snips_bslined(:,:,i)=snips(:,:,i)-entry_r_bsl;
                                        run_averg=cat(3,run_averg,snips(:,:,i));
                                        run_averg_bsl=cat(3,run_averg_bsl,snips_bslined(:,:,i));
                                               
                                    end
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
                        end
                    end  
                end

                end
            end    
        end 
    end
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

save("events_averg_shuf_200.mat", "events_averg",'-v7.3');


