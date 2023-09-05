% analysis including behavior
%trials separated in runs (2+) or singles
clearvars -except ci_data; close all; 
% Analysis of ci_data file, creates events_averg data file
% 06.03.2023
%% 
skipSessionFrames = 20; %nr of missed frames when sessions are skipped
skippedCounter = 0;FrameSkipped = {}; %counts how many sessions are skipped bc nr of missed frames is too high
skipDrink = 0; skippedDrinkTrials = {};% logs nr and id of drink trials skipped bc of faulty sensor
DrinkNonCons = {};FoodNonCons = {};RunNonCons = {};
%DrinkEntryAll = {};FoodEntryAll = {};RunEntryAll = {};SocEntryAll = {};ExpEntryAll={};
Entry_clsf = {};
nr_ef = {}; confus = zeros(6,6,3);nr_entries = zeros(1,6,3);
fn=fieldnames(ci_data.bsl);
%loop through the fields
for ii=1: numel(fn)%animal
    fn1=fieldnames(ci_data.bsl.(fn{ii}));
    ed_averg=[];ef_averg=[];es_averg=[];ee_averg=[];er_averg=[];drink_averg=[];eat_averg=[];blockend_averg=[];soc_averg=[];expl_averg =[];
    exd_averg=[];exf_averg=[];exs_averg=[];exe_averg=[];exr_averg=[];imstop_averg=[];run_averg=[];soc_averg_m=[];
    ed_averg_bsl=[];ef_averg_bsl=[];es_averg_bsl=[];ee_averg_bsl=[];er_averg_bsl=[];drink_averg_bsl=[];eat_averg_bsl=[];blockend_averg_bsl=[];soc_averg_bsl=[];expl_averg_bsl =[];
    exd_averg_bsl=[];exf_averg_bsl=[];exs_averg_bsl=[];exe_averg_bsl=[];exr_averg_bsl=[];imstop_averg_bsl=[];run_averg_bsl=[];soc_averg_bsl_m=[];
    Entry_clsf = {};
    DrinkEntryAll = {};FoodEntryAll = {};RunEntryAll = {};SocEntryAll = {};ExpEntryAll={};
    for j=1: numel(fn1)%day
        fn2=fieldnames(ci_data.bsl.(fn{ii}).(fn1{j}));
        for k=1: numel(fn2)%session
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
                hw_time = ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2}).(6);
                raw_f_res=ci_data.bsl.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
                df_f_trace = dFoF(raw_f_res);%df/f
                raw_f= zscore(df_f_trace,0,1);% zscored using sample sd
                if size(new_event_frames,1) < 4
                    s = s+1;
                else    
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
                ef_events = find(event_type=='enter_feed');
                er_events = find(event_type=='enter_run');
                ee_events = find(event_type=='enter_explore');
                es_events = find(event_type=='enter_social');
                be_events = find(event_type=='block_end');
                entries_total=sort([ef_events;er_events; ee_events; es_events; ed_events]);
                entries_total_be=sort([ef_events;er_events; ee_events; es_events; ed_events;be_events]);

                if event_type((entries_total(1))+1) == "run" |  event_type((entries_total(1))+1) == "retrieve_pellet"| event_type((entries_total(1))+1) == "drink"
                    consum = true;
                else
                    consum = false;
                end
                e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(1))),char(event_type((entries_total(1)))),"s", consum};%add first entry to list, always a single
                Entry_clsf = [Entry_clsf; e_clsf];
        
                
                for ev=2:size(entries_total,1) % check for singles and runs( last entries can only be a certain type of run

                    if ev <= size(entries_total,1)-1 && event_type((entries_total(ev))) ~= event_type((entries_total(ev-1))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev+1)))
                        clsf = "s" ;%single
                    elseif ev <= size(entries_total,1)-1 && event_type((entries_total(ev))) ~= event_type((entries_total(ev-1))) && event_type((entries_total(ev))) == event_type((entries_total(ev+1)))
                        clsf = "sbr" ; %single before run  
                    
                    elseif ev <= size(entries_total,1)-1 && ev > 2 && event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev-2))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev+1)))%dublett
                       clsf = "r2" ;%duplett
                    elseif ev == 2  && event_type((entries_total(ev))) == event_type((entries_total(ev-1)))  && event_type((entries_total(ev))) ~= event_type((entries_total(ev+1)))%dublett
                       clsf = "r2" ;%duplett on second entry   
                    elseif ev > size(entries_total,1)-1 && event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev-2)))%dublett on last one
                       clsf = "r2" ;%duplett on last entry 

                    elseif  ev <= size(entries_total,1)-1 && ev > 3 && event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) == event_type((entries_total(ev-2))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev-3))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev+1)))%triplett
                       clsf = "r3" ;%triplett
                    elseif  ev == 3 && event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) == event_type((entries_total(ev-2))) &&  event_type((entries_total(ev))) ~= event_type((entries_total(ev+1)))%triplett
                       clsf = "r3" ;%triplett on 3rd entry   
                    elseif  ev > size(entries_total,1)-1 && event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) == event_type((entries_total(ev-2))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev-3)))%triplett
                       clsf = "r3" ;%triplett on last entry 

                    elseif  ev <= size(entries_total,1)-1 && ev > 4 && event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) == event_type((entries_total(ev-2))) && event_type((entries_total(ev))) == event_type((entries_total(ev-3))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev+1)))%quadruple
                       clsf = "r4" ;%quadruple
                    %elseif  ev == 4 && event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) == event_type((entries_total(ev-2))) && event_type((entries_total(ev))) == event_type((entries_total(ev-3))) && event_type((entries_total(ev))) ~= event_type((entries_total(ev+1)))%quadruple
                    %   clsf = "r4" ;%quadruple on 4th entry   
                    elseif  event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) == event_type((entries_total(ev-2))) && event_type((entries_total(ev))) == event_type((entries_total(ev-3))) %quadruple
                       clsf = "r4" ;%quadruple on last entry   


                    elseif  event_type((entries_total(ev))) == event_type((entries_total(ev-1))) && event_type((entries_total(ev))) == event_type((entries_total(ev+1)))%part of a run
                       clsf = "r" ; %part of run
                    else
                       clsf = "s" ;
                    end
                    if event_type((entries_total(ev))+1) == "run" |  event_type((entries_total(ev))+1) == "retrieve_pellet"| event_type((entries_total(ev))+1) == "drink"
                        consum = true;
                    else
                        consum = false;
                    end
                    e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(ev))) ,char(event_type((entries_total(ev)))), clsf, consum};
                    Entry_clsf = [Entry_clsf; e_clsf];
     
                end 

                for ev=2:size(entries_total,1)
                    if event_type((entries_total(ev))) ~= event_type((entries_total(ev-1))) && count == 1 %if not like one before and not end of a run its a single
                       clsf = "s" ;%single
                    elseif event_type((entries_total(ev))) ~= event_type((entries_total(ev-1))) && count > 1 %end of a run
                       clsf = ["r" num2str(count)];

    


                    elseif event_type((entries_total(ev))) == event_type((entries_total(ev-1)))%if the same as before
                       count = count + 1;



                       
                    elseif event_type((entries_total(ev))) ~= event_type((entries_total(ev-1))) && count > 1
                       clsf = "s" ;%single   
                    



                    end
                    e_clsf = {fn(ii), fn1(j), s, hw_time((entries_total(ev))) ,char(event_type((entries_total(ev)))), clsf, consum};
                    Entry_clsf = [Entry_clsf; e_clsf];
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
                    ev_index = hw_time(find(event_type==t));
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
                                    ed_all = {fn(ii), fn1(j), s, i, ev_index(i)};
                                    DrinkEntryAll = [DrinkEntryAll; ed_all];

                                elseif trig>win+1 & p==2
                                    snips(:,:,i)=raw_f(trig-win:trig,:); % entry food
                                    entry_f_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
                                    ef_averg=cat(3,ef_averg,snips(:,:,i));
                                    ef_averg_bsl=cat(3,ef_averg_bsl,snips_bslined(:,:,i));
                                    ef_all = {fn(ii), fn1(j), s, i, ev_index(i)};
                                    FoodEntryAll = [FoodEntryAll; ef_all];
                                elseif trig>win+1 & p==3
                                    snips(:,:,i)=raw_f(trig-win:trig,:); % entry social
                                    entry_s_bsl = mean(snips(1:10,:,i),1);
                                    snips_bslined(:,:,i)=snips(:,:,i)-entry_s_bsl;
                                    es_averg=cat(3,es_averg,snips(:,:,i));
                                    es_averg_bsl=cat(3,es_averg_bsl,snips_bslined(:,:,i));
                                    es_all = {fn(ii), fn1(j), s, i, ev_index(i)};
                                    SocEntryAll = [SocEntryAll; es_all];
                    
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
                                    ee_all = {fn(ii), fn1(j), s, i, ev_index(i)};
                                    ExpEntryAll = [ExpEntryAll; ee_all];

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
                                    er_all = {fn(ii), fn1(j), s, i, ev_index(i)};
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
    % ensure that number of trials align
    Entry_all = cell2table(cat(1,FoodEntryAll, RunEntryAll, ExpEntryAll, SocEntryAll, DrinkEntryAll),"VariableNames",["animal" "day" "session" "nr" "hw"]);
    Entry_all = Entry_all(:,[1 2 3 5]);
    EntryClsf = cell2table(Entry_clsf,"VariableNames",["animal" "day" "session" "hw" "type" "clasf" "consum"]);EntryCl = EntryClsf(:,1:4);
    Lia = ismember(EntryCl,Entry_all);
    EntryClsf_2 = EntryClsf(Lia,:);
    % saves in struct according to animal
    types = ["ef", "er","ee","es","ed"];  
    for ty = 1:size(types,2)
         animal =char(fn{ii});type = char(types{ty});
         Clsf_entry.(animal).(type) = EntryClsf_2(EntryClsf_2.type == entry_types(ty), 6:7);
    end
end



%% Cons and non consumption trials

edNCg2 = find(table2array(Clsf_entry.g2.ed(:,2)) == false);
edNCg5 = find(table2array(Clsf_entry.g5.ed(:,2)) == false);
edNCg4 = find(table2array(Clsf_entry.g4.ed(:,2)) == false);
efNCg2 = find(table2array(Clsf_entry.g2.ef(:,2)) == false);
efNCg5 = find(table2array(Clsf_entry.g5.ef(:,2)) == false);
efNCg4 = find(table2array(Clsf_entry.g4.ef(:,2)) == false);


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



%% split trials in singles and runs
% code to sort the imaging trials into runs and singles saved in
% events_classif
fn = fieldnames(Clsf_entry);
for ii=1: numel(fn)%animal
    fn1=fieldnames(Clsf_entry.(fn{ii}));
    for jj = 1:numel(fn1) %entry type
        Singles = find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "s");
        Sbr = find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "sbr");
        Singles = cat(1,Singles,Sbr);
        evnt = char(strcat(fn1(jj), "_averg_bsl"));
        events_classif.(fn{ii}).(fn1{jj}).singles = events_averg.(fn{ii}).(evnt)(:,:,Singles);
        Runs = events_averg.(fn{ii}).(evnt); Runs(:,:, Singles) =[];
        events_classif.(fn{ii}).(fn1{jj}).Runs =Runs;
    end
end   
% code to count singles and different type of runs
fn = fieldnames(Clsf_entry);
for ii=1: numel(fn)%animal
    fn1=fieldnames(Clsf_entry.(fn{ii}));
    for jj = 1:numel(fn1) %entry type
        SRevents(1,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "s"),1); % single is just s without sbr
        SRevents(2,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "r2"),1);
        SRevents(3,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "r3"),1);
        SRevents(4,1)= size(find(table2array(Clsf_entry.(fn{ii}).(fn1{jj})(:,1)) == "r4"),1);
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

%calculate correct confusion matrix values
C_matrix = [];
for a = 1:size(confus,3)
    for r = 1:6
        for c = 1:6
            C_matrix(r,c,a) = confus(r,c,a)/nr_entries(1,r,a);
        end  
    end
end
C_matrix_averg = [];C_matrix_av = 0;
for r = 1:6
       for c = 1:6
           for a=1:size(C_matrix,3)
                C_matrix_av =  C_matrix_av + C_matrix(r,c,a);
           end
           C_matrix_averg (r,c) = C_matrix_av/size(C_matrix,3)
           C_matrix_av = 0;
       end  
end

T_Cmatrix_av = array2table(C_matrix_averg,"RowNames",["ef" "er" "ee" "es" "ed" "h"], "VariableNames",["ef" "er" "ee" "es" "ed" "h"]);
entriesBehavior.nr_entries = nr_entries;
entriesBehavior.counts = confus;
entriesBehavior.matrix = C_matrix;
entriesBehavior.matrix_averg = C_matrix_averg;
entriesBehavior.matrix_tbl = T_Cmatrix_av;




%save("entriesBehavior.mat", "entriesBehavior");
%save("events_classif.mat", "events_classif");
%save("events_averg.mat", "events_averg");

