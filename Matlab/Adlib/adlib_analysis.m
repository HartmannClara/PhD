%%analzsis of adlib data 
%1.11.23
%takes in ci_data.adlib_corr with corrected missed frames
clearvars -except Adlib_pellets Adlib_drinks; close all;
load('ci_data_adlib_2604.mat');
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
%                 missed_frames = ci_data.adlib.(fn{ii}).(fn1{j})(s).(fn2{k}).missed_frames;
%                 if numel(missed_frames) >= skipSessionFrames
%                     skippedCounter = skippedCounter + 1;
%                     fr_skip = {fn(ii), fn1(j), s};
%                     FrameSkipped = [FrameSkipped; fr_skip];
%                     s = s+1;
%                 else
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
%                 clear PodFrames PodTimes
%                 for et=1:size(AllEntries,2)%event type (ef er ee es ed imaging starts)
%                    NrEt = AllEntries{1,et};
%                    NrEx = AllExits{1,et};
%                    PodFrames = [];
%                    for ev = 1:size(NrEt,1) %nr of events
%                        PodFrames(ev,1) =  (new_event_frames(NrEx(ev)) - new_event_frames(NrEt(ev)))/10; %frames in the pod divided by framerate
%                    end
%                    PodTimes{1,et} = PodFrames;% time in pod (ef er ee es ed session length)
%                 end
%                 TimeSpecs = [TimeSpecs; PodTimes];%pod times: food, run, explore, social,drink, whole block

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
                %end
            end
        end
    end
end

