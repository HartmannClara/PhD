%Drivemaze etho and ca imaging plotter 
%2022March session_starts
close all, clear all
%animals=table(['94331472'; '328340178184'; '335490249236'; '328340226232']);%
animals='3354903011';
%'94331472 grin1 female';
f1=figure;

% f4=figure;
% f5=figure;
% for 
a=1;%:3
clearvars -except a f1 f2 f3 f4 f5 animals binned b_trials b_dur binned_exp
% animal=cellstr(table2cell(animals(a,1)));
animal=animals;
% filename=['C:\Data\Drivemaze\Drivemaze_imaging_grin1' animal{1, 1} '_events.csv'];
filename=['D:\Exp_2_DM_exploration\PFC-LH\g5\' animal '_events.csv'];
%% Import csv
 opts = delimitedTextImportOptions("NumVariables", 7);
 opts.DataLines = [2, Inf];
 opts.Delimiter = ",";
 opts.VariableNames = ["Date_Time", "amount_consumed", "latency_to_consumption", "Type", "frame", "hardware_time", "experiment"];
 opts.VariableTypes = ["string", "double", "double", "categorical", "double", "double", "categorical"];
 opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
 opts = setvaropts(opts, [1, 4], "EmptyFieldRule", "auto");
 opts.ExtraColumnsRule = "ignore";
 opts.EmptyLineRule = "read";
 events = readtable(filename, opts);
 clear opts

x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
%% chop
chop_from=datetime('2022-09-23 14:37:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%chop the day
chop_to=datetime('2022-09-23 20:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];


x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);
%% cleanup

% find imaging starts
%u1=find(event_type=='frame');
u1=find(event_type=='imaging_stop');
%f_num=events.(3);
frames =events.(5);
abc = u1 - 1; 
imaging_ends=u1(find(frames(u1) > 1 ));%frames were aquired

u2=find(event_type=='imaging_start');
imaging_starts=[];%accumulator
for i=1:size(imaging_ends,1)
    u3=u2(find(u2 < imaging_ends(i)));%finds the previous imaging_starts
    imaging_starts=[imaging_starts; u3(end)];%adds the closest imaging start to list
end


%% plot
x=datetime(events.(1),'Format','yyyy-MM-dd HH:mm:ss.SSS');
y=NaN(size(x));
%frameplot=NaN(size(x));
%frameplot(find(event_type=='frame'))=f_num(find(event_type=='frame'));
%enter
y(find(event_type=='leave_nest'))=2;
y(find(event_type=='block_start'))=3.5; %change this to 2 before 28.1.22
y(find(event_type=='enter_feed'))=4; %food entry
enter(3).e=x(find(event_type=='enter_feed'));
y(find(event_type=='enter_run'))=3; %wheel entry
enter(2).e=x(find(event_type=='enter_run'));
y(find(event_type=='enter_social'))=4.5; %
enter(4).e=x(find(event_type=='enter_social'));
y(find(event_type=='enter_drink'))=2.5; %
enter(5).e=x(find(event_type=='enter_drink'));
y(find(event_type=='enter_explore'))=5; 
enter(1).e=x(find(event_type=='enter_explore'));
%consume
y(find(event_type=='retrieve_pellet'))=4;%pellet retrieval
y(find(event_type=='drink'))=2.5; %
y(find(event_type=='run'))=3; %wheel entry
pel=x(find(event_type=='retrieve_pellet'));
drink=x(find(event_type=='drink'));
run=x(find(event_type=='run'));
%exit
y(find(event_type=='block_end'))=2;
y(find(event_type=='block_available'))=2;
y(find(event_type=='exit_feed'))=3.5;
exit(3).e=x(find(event_type=='exit_feed'));
y(find(event_type=='exit_drink'))=3.5;
exit(5).e=x(find(event_type=='exit_drink'));
y(find(event_type=='exit_run'))=3.5;
exit(2).e=x(find(event_type=='exit_run'));
y(find(event_type=='exit_social'))=3.5;
exit(4).e=x(find(event_type=='exit_social'));
y(find(event_type=='exit_explore'))=3.5;
exit(1).e=x(find(event_type=='exit_explore'));

session_starts=find(event_type=='initialize');
%session_ends=find(event_type=='end session');
yy=fillmissing(y,'linear');

figure(f1);
subplot(1,1,a),plot([x(1) x(end)],[3.5 3.5],'b');hold on%decisionpoint
subplot(1,1,a),plot([x(1) x(end)],[2 2],'b');%nest
subplot(1,1,a),plot(x,yy,'k');hold on
subplot(1,1,a),plot(x,y,'ko');hold on
%subplot(1,1,a),plot(x,fed_pos,'b');
%subplot(1,1,a),plot(x,frameplot/1000,'g');hold on

title(animal);
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

%%
explore = numel(enter(1).e);
running = numel(enter(2).e);
food = numel(enter(3).e);
social = numel(enter(4).e);
water = numel(enter(5).e);
total_entries = explore+water+social+running+food;

%% 
starts = find(event_type == "block_start");
ends=find(event_type=='block_end');

%block duration and trial number
for i=1:size(starts,1)-1 %clean up
    if size(starts,1)>i
        if ends(i)>starts(i+1)
            starts(i)=[];
        end
    end
end
for i=1:size(ends,1)-1 %clean up
    if starts(i)>ends(i)
        ends(i)=[];
    end
end
%%
b_food_dur=[];

%session_starts
clear i
for i=1:size(ends,1) %blocks
    b_dur.dur(i)=x(ends(i))-x(starts(i));%duration of block
    events_block=events(starts(i):ends(i),:);%what happens in the block
    event_type_block=events_block.(4);
    b_trials.e(i)=length(find(event_type_block=='enter_explore'));
    b_trials.r(i)=length(find(event_type_block=='enter_run'));
    b_trials.f(i)=length(find(event_type_block=='enter_feed'));
    b_trials.s(i)=length(find(event_type_block=='enter_social'));
    b_trials.w(i)=length(find(event_type_block=='enter_drink'));
    b_trials.total(i)=b_trials.e(i)+b_trials.r(i)+b_trials.f(i)+b_trials.s(i)+b_trials.w(i);
    b_trials.entryrate(i) = (b_trials.total(i)/seconds(b_dur.dur(i)))*60;
    
    
    
    %high res analyses

    b_trials_explore=find(event_type_block=='enter_explore'); % checking if every pod entry has a pod exit
    b_trials_explore_end=find(event_type_block=='exit_explore');
    if length(b_trials_explore) ~= length(b_trials_explore_end)
        'problem inconsistent explore trials in this block'
        'correcting by adding next read'
        b_trials_explore_end=[b_trials_explore_end;b_trials_explore(end)+2];
    end

    b_trials_run=find(event_type_block=='enter_run');
    b_trials_run_end=find(event_type_block=='exit_run');
    if length(b_trials_run) ~= length(b_trials_run_end)
        'problem inconsistent run trials in this block'
        'correcting by adding next read'
        b_trials_run_end=[b_trials_run_end; b_trials_run(end)+2];
    end
    
    b_trials_food=find(event_type_block=='enter_feed');
    b_trials_food_end=find(event_type_block=='exit_feed');
    if length(b_trials_food) ~= length(b_trials_food_end)
        'problem inconsistent food trials in this block'
        'correcting by adding next read'
        b_trials_food_end=[b_trials_food_end;b_trials_food(end)+2];
    end

    b_trials_social=find(event_type_block=='enter_social');
    b_trials_social_end=find(event_type_block=='exit_social');
    if length(b_trials_social) ~= length(b_trials_social_end)
        'problem inconsistent social trials in this block'
        'correcting by adding next read'
        b_trials_social_end=[b_trials_social_end;b_trials_social(end)+2];
    end
   
    b_trials_water=find(event_type_block=='enter_drink');
    b_trials_water_end=find(event_type_block=='exit_drink');
    if length(b_trials_water) ~= length(b_trials_water_end)
        'problem inconsistent water trials in this block'
        'correcting by adding next read'
        b_trials_water_end=[b_trials_water_end; b_trials_water(end)+2];
    end
    
    b_trials_total=sort([b_trials_explore;b_trials_run; b_trials_food; b_trials_social; b_trials_water]);
    b_trials_total_end=sort([b_trials_explore_end;b_trials_run_end; b_trials_food_end; b_trials_social_end; b_trials_water_end]);
    if size(b_trials_total,1)>1

        j=0;%switch counter
        clear exploit
        ex_ind=1; 
        exploit(ex_ind)=0; %exploit counter
        b_trial_dur=[];%trial duration accumulator
        b_choice_dur=[];%choice duration accumulator
        for e=2:size(b_trials_total,1) % evaluating whether next choice was a switch: is it like the one before
            if find(b_trials_explore==b_trials_total(e)) & find(b_trials_run==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_explore==b_trials_total(e)) & find(b_trials_food==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_explore==b_trials_total(e)) & find(b_trials_social==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_explore==b_trials_total(e)) & find(b_trials_water==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;    

            elseif find(b_trials_run==b_trials_total(e)) & find(b_trials_explore==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_run==b_trials_total(e)) & find(b_trials_food==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_run==b_trials_total(e)) & find(b_trials_social==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_run==b_trials_total(e)) & find(b_trials_water==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;


            elseif find(b_trials_food==b_trials_total(e)) & find(b_trials_explore==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_food==b_trials_total(e)) & find(b_trials_run==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_food==b_trials_total(e)) & find(b_trials_social==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_food==b_trials_total(e)) & find(b_trials_water==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;


            elseif find(b_trials_social==b_trials_total(e)) & find(b_trials_explore==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_social==b_trials_total(e)) & find(b_trials_run==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_social==b_trials_total(e)) & find(b_trials_food==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_social==b_trials_total(e)) & find(b_trials_water==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;


            elseif find(b_trials_water==b_trials_total(e)) & find(b_trials_explore==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_water==b_trials_total(e)) & find(b_trials_run==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_water==b_trials_total(e)) & find(b_trials_food==b_trials_total(e-1))
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_water==b_trials_total(e)) & find(b_trials_social==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;

            else %it was the same then so exploitation index +1
              
                exploit(ex_ind)=exploit(ex_ind)+1;
                

            end
            b_trial_dur=[b_trial_dur; milliseconds(x(b_trials_total_end(e))-x(b_trials_total(e)))];
            b_choice_dur=[b_choice_dur; milliseconds(x(b_trials_total(e))-x(b_trials_total_end(e-1)))];  
            if find(b_trials_food==b_trials_total(e))%plots duration of food trials
                b_food_dur=[b_food_dur; milliseconds(x(b_trials_total_end(e))-x(b_trials_total(e)))];
            end    


        end
        b_trials.mean_choice_dur(i)=nanmean(b_choice_dur);
        b_trials.mean_trial_dur(i)=nanmean(b_trial_dur);

    else
        j=NaN;
        exploit=NaN;
        b_trials.mean_trial_dur(i)=NaN;
        b_trials.mean_choice_dur(i)=NaN;
    end
    
%first entry
    if event_type_block(b_trials_total(1)) == 'enter_feed'
            b_trials.first_entry(i) = 1;
    elseif event_type_block(b_trials_total(1)) == 'enter_run'
            b_trials.first_entry(i) = 2;
    elseif event_type_block(b_trials_total(1)) == 'enter_explore'
            b_trials.first_entry(i) = 3  ;  
    elseif event_type_block(b_trials_total(1)) == 'enter_social'
            b_trials.first_entry(i) = 4;
    elseif event_type_block(b_trials_total(1)) == 'enter_drink'
            b_trials.first_entry(i) = 5   ; 
    end  
   
    b_trials(a).switch_index(i)=j/(size(b_trials_total,1)-1);%total number of trials - the first trial
    b_trials(a).exploit_index_max(i)=max(exploit);
    b_trials(a).exploit_index_mean(i)=mean(exploit);
    b_trials(a).exploit_index_min(i)=min(exploit);
    b_trials(a).switches(i)=j;
    b_trials(a).start_time(i)=x(starts(i));
       
    %b_trials(a).revolutions(i)=sum(revolutions(starts(i):ends(i)));

  
end

%% check for faulty blocks
f3=figure;%plot total entries
figure(f3);
xx=x(ends);
subplot(size(animals,1),1,a),plot(xx,b_trials.total,'ko','MarkerFaceColor','k');hold on
%title(animal);
%yticks([0.5 1.5 2 2.5 3 3.5 4 4.5 5]);
%yticklabels({'frames','frames','nest','drink','run','decision point','food','social','explore'});
ylabel('total entries');

%get faulty block index
remove = find(b_trials.total > 50)



%remove faulty blocks
clear i; clear ii;
fn = fieldnames(b_trials);
for i=1:15(fn)
    for ii = 1:size(remove,2)
        b_trials.(fn{i})(remove(ii)) = NaN;
        %b_trials.(fn{i})(81) = NaN;
        %b_trials.(fn{i})(69) = NaN;
    end    
end    



%% figures

pod_entries(1, :) = b_trials.e;
pod_entries(2, :) = b_trials.r;
pod_entries(3, :) = b_trials.f;
pod_entries(4, :) = b_trials.s;
pod_entries(5, :) = b_trials.w;
pod_entries(6, :) = ends;


f2 = figure;
leg = {'food'};%{'explore', 'run', 'food', 'social', 'water', 'total'};
figure(f2);
xx=x(ends);
%plot(xx,pod_entries(1,:),'MarkerFaceColor','g');hold on;
%plot(xx,pod_entries(2,:),'MarkerFaceColor','m');
plot(xx,pod_entries(3,:),'MarkerFaceColor','y');
%plot(xx,pod_entries(4,:),'MarkerFaceColor','r');
%plot(xx,pod_entries(5,:),'MarkerFaceColor','w');
%plot(xx,b_trials.total,'ko','MarkerFaceColor','k');
%set(barb,{'FaceColor'},{'g';'m'; 'y';'r'; 'b'});
legend(leg)
ylabel('Trials')


%subplot(size(animals,1),1,a),plot(xx,b_trials.total,'ko','MarkerFaceColor','k');hold on
%ylabel('Trials');
%plot(xx,b_trials(a).e, '.');
%plot(xx,b_trials(a).r, '.');
%plot(xx,b_trials(a).f, '.');
%plot(xx,b_trials(a).s, '.');
%plot(xx,b_trials(a).w, '.');

f3=figure;
figure(f3);
xx=x(ends);
subplot(size(animals,1),1,a),plot(xx,b_trials.switch_index,'ko','MarkerFaceColor','k');hold on

%title(animal);
%yticks([0.5 1.5 2 2.5 3 3.5 4 4.5 5]);
%yticklabels({'frames','frames','nest','drink','run','decision point','food','social','explore'});
ylabel('switch index');
%ylim([-40 5.5]);



f4= figure;
figure(f4);
plot(xx,b_trials.entryrate,'ko','MarkerFaceColor','k');
%title(animal);
%yticks([0.5 1.5 2 2.5 3 3.5 4 4.5 5]);
%yticklabels({'frames','frames','nest','drink','run','decision point','food','social','explore'});
ylabel('entryrate per block (e/min)');
%ylim([-40 5.5]);
%yyaxis right;
%plot(xx,b_trials.total,'^','MarkerFaceColor','r');
%ylabel('total entries');

%%
% mean choice duration
f4=figure;
figure(f4);
xx=x(ends);
subplot(size(animals,1),1,a),plot(xx,b_trials.mean_choice_dur,'ko','MarkerFaceColor','k');hold on
ylabel('mean choice duration');

%mean trial duration
f5=figure;
figure(f5);
xx=x(ends);
xxx=1:39;
%subplot(size(animals,1),1,a),plot(xx,b_trials.mean_trial_dur,'ko','MarkerFaceColor','k');hold on
%ylabel('mean trial duration');
subplot(size(animals,1),1,a),plot(b_food_dur,xxx,'ko','MarkerFaceColor','k');hold on% food trial duration
xlabel('food trial duration');


%% first entry in block

f6=figure;
figure(f6);
xx=x(ends);
subplot(size(animals,1),1,a),plot(xx,b_trials.first_entry,'ko','MarkerFaceColor','k');hold on
ylabel('first entry');
yticks([1  2  3  4  5]);

%% separate
% bsl1_from=datetime('2022-07-06 9:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% switch1_from=datetime('2022-07-07 15:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% bsl2_from=datetime('2022-07-11 9:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
% switch2_from=datetime('2022-07-12 9:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');


%maybe just sorting the trials?

%bsl_trials = find(b_trials.start_time<)
%switch_trials = 





%% consumption



%%
%plot trial nr into
%first entry

%%
% switch index per block: nr of switches/total nr of entries in that block
% duration and nr of exploits (repeated entering and exploiting of one pod)
% latencies to consume
% time spent in pods
%


