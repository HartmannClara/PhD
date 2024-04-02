%%analzsis of adlib data 
%1.11.23
%takes in ci_data.adlib_corr with corrected missed frames
clear all
close all
load('ci_data_adlib_corr.mat');
%%
factor = 0.9; % factor by which the neuropil is subtracted 0.7-0.9
skipSessionFrames = 30; %nr of missed frames when sessions are skipped
skippedCounter = 0;FrameSkipped = {};DecSkipped = {}; %counts how many sessions are skipped bc nr of missed frames is too high
skipDrink = 0; skippedDrinkTrials = {};% logs nr and id of drink trials skipped bc of faulty sensor
DrinkEntryNC = {};FoodEntryNC = {};RunEntryNC = {};
DrinkEntryAll = {};FoodEntryAll = {};RunEntryAll = {};DrinkAll = {};EatAll = {};
FoodExitAll ={};RunExitAll ={};ExpExitAll ={};SocExitAll ={};DrinkExitAll ={}; ImgStopAll = {};

Adlib_pellets = [];%collects snips and timestamps of pellet episodes
Adlib_drinks = [];%collects snips and timestamps of drink episodes

fn=fieldnames(ci_data.adlib);
%loop through the fields
for ii=1: numel(fn)
    fn1=fieldnames(ci_data.adlib.(fn{ii}));
    ed_averg=[];ef_averg=[];drink_averg=[];eat_averg=[];
    exd_averg=[];exf_averg=[];
    ed_averg_bsl=[];ef_averg_bsl=[];drink_averg_bsl=[];eat_averg_bsl=[];
    exd_averg_bsl=[];exf_averg_bsl=[];
    TimeSpecs = [];DecSpecs = [];
    %ed_averg_nc=[];ef_averg_nc=[];er_averg_nc=[];ed_averg_nc_bsl=[];ef_averg_nc_bsl=[];er_averg_nc_bsl=[];
    for j=1: numel(fn1)%day
        fn2=fieldnames(ci_data.adlib.(fn{ii}).(fn1{j}));
        for k=1: numel(fn2)
            for s=1:size(ci_data.adlib.(fn{ii}).(fn1{j}),2)%access the data sessions
                fn3= fieldnames(ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}));
                missed_frames = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).missed_frames;
                if numel(missed_frames) >= skipSessionFrames
                    skippedCounter = skippedCounter + 1;
                    fr_skip = {fn(ii), fn1(j), s};
                    FrameSkipped = [FrameSkipped; fr_skip];
                    s = s+1;
                else
                latency =  ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(3); % latency to consume in ms   
                new_event_frames= ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(5);% in events column 5 (frame)
                event_type = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{3}).(4);
                raw_f_res=ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
                Fneu = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2});

                subneuro_trace = subneuro(raw_f_res,Fneu,factor); % subtracts neuropil trace from raw trace by factor x subneuro(rawF, Fneuropil, factor)
                subneuro_trace(:,(size(subneuro_trace,2)+1)) = mean(Fneu,2);%last "cell" is mean neuropil trace
                
                
                df_f_trace = dFoF(subneuro_trace);%df/f
                raw_f= zscore(df_f_trace,0,1);% zscored using sample sd
                mean_raw_f = mean(raw_f(:,1:end-1),2);
                raw_f= [raw_f,mean_raw_f];

             %% check if session is usable
              % collect all event types in single variables 
                ee_events = find(event_type=='enter_explore'); es_events = find(event_type=='enter_social');
                exf_events = find(event_type=='exit_feed');exd_events = find(event_type=='exit_drink');exr_events = find(event_type=='exit_run');
                exs_events = find(event_type=='exit_social');exe_events = find(event_type=='exit_explore');
                istarts = find(event_type=='imaging_start');B_ends = find(event_type=='block_end');
                ed_events = find(event_type=='enter_drink');ef_events = find(event_type=='enter_feed');er_events = find(event_type=='enter_run');

                AllEntries = {ef_events,er_events, ee_events,es_events,ed_events,istarts};AllEntriesList = sort(vertcat(AllEntries{:}));
                AllExits = {exf_events,exr_events, exe_events,exs_events,exd_events,B_ends}; AllExitsList = sort(vertcat(AllExits{:}));
                % legth of time in desicion zone
                AllExitsList = cat(1,AllEntriesList(1,1), AllExitsList(:,1));%add imaging start as first exit
                AllEntriesList(1) = [];% remove imaging start from entries
                AllEntriesList(end+1,1) =  AllExitsList(end,1);% add block end as last entry into something
                AllExitsList(end)= [];% remove block end from exits


                clear DecFrames DecTypes 
                for et=1:size(AllExitsList,1)%event type
                    DecFrames(et,1) =  (new_event_frames(AllEntriesList(et)) - new_event_frames(AllExitsList(et)))/10;% check frames from exit to next entry
                    DecTypes(et,1) = event_type(AllExitsList(et)); % save what came before and what after
                    DecTypes(et,2) = event_type(AllEntriesList(et));
                end
                %check if there is less than 1.8s between desicions(entrz to exit and exit to entry) = error
                err = find(DecFrames < 1.8);
                if err ~= 0 
                   skippedCounter = skippedCounter + 1;
                   dec_skip = {fn(ii), fn1(j), s};
                   DecSkipped = [DecSkipped; dec_skip];
                   s = s+1;% if there is an error skip the session
                else
                    Desicion = {DecFrames, DecTypes, s};% seconds in desicion time and identidy of previous and next entry
                    DecSpecs = [DecSpecs; Desicion];% Decision zone: time, exit to entry types
                end
                %check how long the entries are
                clear PodFrames PodTimes
                for et=1:size(AllEntries,2)%event type (ef er ee es ed imaging starts)
                   NrEt = AllEntries{1,et};
                   NrEx = AllExits{1,et};
                   PodFrames = [];
                   for ev = 1:size(NrEt,1) %nr of events
                       PodFrames(ev,1) =  (new_event_frames(NrEx(ev)) - new_event_frames(NrEt(ev)))/10; %frames in the pod divided by framerate
                   end
                   PodTimes{1,et} = PodFrames;% time in pod (ef er ee es ed session length)
                end
                TimeSpecs = [TimeSpecs; PodTimes];%pod times: food, run, explore, social,drink, whole block

                %%
                % check if lick sensor is ok, otherwise skip drink snips
                ed_event_frames = new_event_frames(find(event_type=='enter_drink'));
                d_event_frames = new_event_frames(find(event_type=='drink'));
                ed_d_frames = cat(1,ed_event_frames, d_event_frames);   
                if size(ed_d_frames,1) == size(unique(ed_d_frames),1)
                    pp = 1:6;
                elseif size(ed_d_frames,1) > size(unique(ed_d_frames),1)
                    pp = [2 4:6];
                    skipped = {fn(ii), fn1(j), s};
                    skippedDrinkTrials = [skippedDrinkTrials; skipped];
                end   

                %% collect adlib sessions
                exd_event_frames = new_event_frames(find(event_type=='exit_drink'));
                ef_event_frames = new_event_frames(find(event_type=='enter_feed'));
                f_event_frames = new_event_frames(find(event_type=='retrieve_pellet'));
                exf_event_frames = new_event_frames(find(event_type=='exit_feed'));
                win = 30 ;%window for after exit used for baselining
                for start=1:numel(exf_event_frames)
                    adlib_f_snip = raw_f(ef_event_frames(start)-25:exf_event_frames(start)+win,:);%cuts all frames from enty -25 frames to exit
                    adlib_f_events = event_type(ef_events(start):exf_events(start),:);%cuts the events in that pod 
                    adlib_f_frames = new_event_frames(ef_events(start):exf_events(start),:);%gives the timestamps of the pellet retrievals
                    % correct timestamps so it matches the snip
                    st =adlib_f_frames(1)-1; %start at 1 
                    for fr =1:numel(adlib_f_frames)
                        adlib_f_frames(fr) = (adlib_f_frames(fr)-st)+25;
                    end
                    NrPel = size(find(adlib_f_events == 'retrieve_pellet'),1);%nr of pellets eaten in that pod
                    PelStats= {adlib_f_snip,adlib_f_events,adlib_f_frames,NrPel,fn{ii}};
                    Adlib_pellets = vertcat(Adlib_pellets, PelStats);
                end    
                
                for start=1:numel(exd_event_frames)
                    adlib_d_snip = raw_f(ed_event_frames(start)-25:exd_event_frames(start)+win,:);%cuts all frames from enty -25 frames to exit
                    adlib_d_events = event_type(ed_events(start):exd_events(start),:);%cuts the events in that pod 
                    adlib_d_frames = new_event_frames(ed_events(start):exd_events(start),:);%gives the timestamps of the pellet retrievals
                    % correct timestamps so it matches the snip
                    st =adlib_d_frames(1)-1; %start at 1 
                    for fr =1:numel(adlib_d_frames)
                        adlib_d_frames(fr) = (adlib_d_frames(fr)-st)+25;
                    end
                    NrDri = size(find(adlib_d_events == 'drink'),1);%nr of pellets eaten in that pod
                    DrinkStats= {adlib_d_snip,adlib_d_events,adlib_d_frames,NrDri,fn{ii}};
                    Adlib_drinks = vertcat(Adlib_drinks, DrinkStats);
                end 
                end
            end
        end
    end
end


%%      %% if sessions ended early skip block end and imaging_stops
              %  if ii == 4 & j ==1 & k == 1 & s== 18
               %     pp = pp(pp~=8);pp = pp(pp~=14);
                %else
                 %   pp=pp;
                %end    
% % events snips
%                 for p=pp
%                     clear snips snips_bslined snips2 snips_bslined2
%                     if p==1
%                         t='enter_drink';
%                         tt='enter drink';
%                     elseif p==2
%                         t='enter_feed';
%                         tt='enter feed';
%                     elseif p==3
%                         t='drink';
%                         tt='drink';
%                     elseif p==4
%                         t='retrieve_pellet';   
%                         tt='retrieve pellet';      
%                     elseif p==5
%                         t='exit_drink';
%                         tt='exit drink';
%                     elseif p==6
%                         t='exit_feed';
%                         tt='exit feed';
%                     
%                     end
%                     
%                     trigs=new_event_frames(find(event_type==t));%trigs are the indexes of events of event type t
%                     eventNr = find(event_type == t);
%                     if size(trigs,1)==0 % if that event does not exist in the session, move on
%                         p=p+1;
%                     else
%                         win = 30;% event window
%                         fte = 25; %frames to entry/exit; how long it takes an animal to go from door to entry beam
% 
%                         for i=1:size(trigs,1)%loop through the events of type t
%                             trig=trigs(i);
% %%%%%%%%%%%%%%%%% for chance triggers
%                             randomNumber = 0.1 + (0.9-0.1) * rand;
%                             trig = round(trig*randomNumber);
%                             evnt = eventNr(i);
%                             
%                                 if trig>win + fte + 1 & p==1 & event_type(evnt+1) == "drink" %% only if there are enough frames before trig
%                                     snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); %entry drink, entries take 3s/30 frames bf door entry (fte before entry beam break)
%                                     entry_d_bsl = mean(snips(1:10,:,i),1);% baselined to 30-20s before door entry
%                                     snips_bslined(:,:,i)=snips(:,:,i)-entry_d_bsl;
%                                     ed_averg=cat(3,ed_averg,snips(:,:,i));
%                                     ed_averg_bsl=cat(3,ed_averg_bsl,snips_bslined(:,:,i));
%                                     ed_all = {fn(ii), fn1(j), s, i};
%                                     DrinkEntryAll = [DrinkEntryAll; ed_all];
%                                 elseif trig>win + fte + 1 & p==1 & event_type(evnt+1) ~= "drink"%% entry drink NON Consumption
%                                     snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); 
%                                     entry_d_bsl = mean(snips(1:10,:,i),1);
%                                     snips_bslined(:,:,i)=snips(:,:,i)-entry_d_bsl;
%                                     ed_averg_nc=cat(3,ed_averg_nc,snips(:,:,i));
%                                     ed_averg_nc_bsl=cat(3,ed_averg_nc_bsl,snips_bslined(:,:,i));
%                                     ed_all = {fn(ii), fn1(j), s, i};
%                                     DrinkEntryAll = [DrinkEntryAll; ed_all];
%                                     ed_nc = {fn(ii), fn1(j), s, i};
%                                     DrinkEntryNC = [DrinkEntryNC; ed_nc];
%                                 elseif trig>win+ fte +1 & p==2 & event_type(evnt+1) == "retrieve_pellet"
%                                     snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); % entry food CONSUMPTION
%                                     entry_f_bsl = mean(snips(1:10,:,i),1);
%                                     snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
%                                     ef_averg=cat(3,ef_averg,snips(:,:,i));
%                                     ef_averg_bsl=cat(3,ef_averg_bsl,snips_bslined(:,:,i));
%                                     ef_all = {fn(ii), fn1(j), s, i};
%                                     FoodEntryAll = [FoodEntryAll; ef_all];
%                                 elseif trig>win+ fte +1 & p==2 & event_type(evnt+1) ~= "retrieve_pellet"
%                                     snips(:,:,i)=raw_f(trig-(fte + win):trig - fte,:); % entry food NON CONSUMPTION
%                                     entry_f_bsl = mean(snips(1:10,:,i),1);
%                                     snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
%                                     ef_averg_nc=cat(3,ef_averg_nc,snips(:,:,i));
%                                     ef_averg_nc_bsl=cat(3,ef_averg_nc_bsl,snips_bslined(:,:,i));
%                                     ef_all = {fn(ii), fn1(j), s, i};
%                                     FoodEntryAll = [FoodEntryAll; ef_all];
%                                     ef_nc = {fn(ii), fn1(j), s, i};
%                                     FoodEntryNC = [FoodEntryNC; ef_nc];
%                                        
%                                 consumptions with beam breaks    
%                                 elseif trig>win+1 & p==3 %drink, takes 3s before and 3s after drink 
%                                     snips(:,:,i)=raw_f(trig-win:trig+win,:);
%                                     d_bsl = mean(snips(1:10,:,i),1); %baselined to first first 10 frames of 60 frame event window (2s before drinking)
%                                     snips_bslined(:,:,i)=snips(:,:,i)-d_bsl;
%                                     drink_averg=cat(3,drink_averg,snips(:,:,i));
%                                     drink_averg_bsl=cat(3,drink_averg_bsl,snips_bslined(:,:,i));
%                                     d_all = {fn(ii), fn1(j), s, i};
%                                     DrinkAll = [DrinkAll; d_all];
%                                 elseif trig>win+1 & p==4 %pellet, takes 3s before and 3s after beam break 
%                                     snips(:,:,i)=raw_f(trig-win:trig+win,:);
%                                     f_bsl = mean(snips(1:10,:,i),1);%takes 1s right after entry before pellet
%                                     snips_bslined(:,:,i)=snips(:,:,i)-f_bsl;
%                                     eat_averg=cat(3,eat_averg,snips(:,:,i));
%                                     eat_averg_bsl=cat(3,eat_averg_bsl,snips_bslined(:,:,i));
%                                     e_all = {fn(ii), fn1(j), s, i};
%                                     EatAll = [EatAll; e_all];
% % exits
%                                    
%                                 elseif trig>win+1 & p==5%exit drink
%                                     snips(:,:,i)=raw_f(trig-win:trig+win,:); %maybe its more interesting what happens after the door is crossed (3s before beam ) and when animal is in center(+win)
%                                     exd_bsl = mean(snips(1:10,:,i),1);
%                                     snips_bslined(:,:,i)=snips(:,:,i)-exd_bsl;
%                                     exd_averg=cat(3,exd_averg,snips(:,:,i));
%                                     exd_averg_bsl=cat(3,exd_averg_bsl,snips_bslined(:,:,i));
%                                     exd_all = {fn(ii), fn1(j), s, i};
%                                     DrinkExitAll = [DrinkExitAll; exd_all];
%                                 elseif trig>win+1 & p==6
%                                     snips(:,:,i)=raw_f(trig-win:trig+win,:); %exit food
%                                     exf_bsl = mean(snips(1:10,:,i),1);
%                                     snips_bslined(:,:,i)=snips(:,:,i)-exf_bsl;
%                                     exf_averg=cat(3,exf_averg,snips(:,:,i));
%                                     exf_averg_bsl=cat(3,exf_averg_bsl,snips_bslined(:,:,i));
%                                     exf_all = {fn(ii), fn1(j), s, i};
%                                     FoodExitAll = [FoodExitAll; exf_all];
%                                 
%                                            
%                                 end                        
%                             end
%                     end 
%                 end
%                 end
%                 end
%             end    
%         end 
%     end
%     %saving into events_averg and Specs structs
%     Specs.(fn{ii}).PodTimes = TimeSpecs;
%     Specs.(fn{ii}).DecTimes = DecSpecs;
%     
%     events_averg.(fn{ii}).ed_averg = ed_averg;
%     events_averg.(fn{ii}).ef_averg = ef_averg;
%     
%     events_averg.(fn{ii}).drink_averg = drink_averg;
%     events_averg.(fn{ii}).eat_averg = eat_averg;
%     
%     events_averg.(fn{ii}).exd_averg = exd_averg;
%     events_averg.(fn{ii}).exf_averg = exf_averg;
%     
%     events_averg.(fn{ii}).ed_averg_nc = ed_averg_nc;
%     events_averg.(fn{ii}).ef_averg_nc = ef_averg_nc;
% 
%     events_averg.(fn{ii}).ed_averg_bsl = ed_averg_bsl;
%     events_averg.(fn{ii}).ef_averg_bsl = ef_averg_bsl;
%     
%     events_averg.(fn{ii}).drink_averg_bsl = drink_averg_bsl;
%     events_averg.(fn{ii}).eat_averg_bsl = eat_averg_bsl;
%     
%     events_averg.(fn{ii}).exd_averg_bsl = exd_averg_bsl;
%     events_averg.(fn{ii}).exf_averg_bsl = exf_averg_bsl;
%     
%     events_averg.(fn{ii}).ed_averg_nc_bsl = ed_averg_nc_bsl;
%     events_averg.(fn{ii}).ef_averg_nc_bsl = ef_averg_nc_bsl;
%     
