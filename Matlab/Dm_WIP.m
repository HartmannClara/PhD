%Drivemaze etho and ca imaging plotter
%2022March 
close all, clear all
%animals=table(['94331472'; '328340178184'; '335490249236'; '328340226232']);%
animals='328340226232';
'94331472 grin1 female';
f1=figure;
f2=figure;
f3=figure;
% f4=figure;
% f5=figure;
% for 
a=1%:3
clearvars -except a f1 f2 f3 f4 f5 animals binned b_trials b_dur binned_exp
% animal=cellstr(table2cell(animals(a,1)));
animal=animals;
% filename=['C:\Data\Drivemaze\Drivemaze_imaging_grin1' animal{1, 1} '_events.csv'];
filename=['C:\Users\cha206\Data\Exp_2_DM_exploration\PFC-LH\g2\1_Raw\' animal '_events.csv'];
%% Import csv
opts = delimitedTextImportOptions("NumVariables", 4);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Date_Time", "amount_consumed", "latency_to_consumption", "Type"];
opts.VariableTypes = ["string", "double", "double", "categorical"];
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
chop_from=datetime('2022-04-03 17:05:40.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
chop_to=datetime('2022-04-13 17:31:04.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

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
u1=find(event_type=='frame');
f_num=events.(3);
imaging_starts=u1(find(f_num(u1)==1));
imaging_ends=[];%accumulator
for i=2:size(imaging_starts,1)
u2=u1(find(u1<imaging_starts(i)));
imaging_ends=[imaging_ends; u2(end)];
end
imaging_ends=[imaging_ends; u1(end)];

%% plot
x=datetime(events.(1),'Format','yyyy-MM-dd HH:mm:ss.SSS');
y=NaN(size(x));
frameplot=NaN(size(x));
frameplot(find(event_type=='frame'))=f_num(find(event_type=='frame'));
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
subplot(1,1,a),plot(x,frameplot/1000,'g');hold on

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
total_entries = explore+water+social+running+food

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


for i=1:size(ends,1) %blocks
    b_dur.dur(i)=x(ends(i))-x(starts(i));%duration
    events_block=events(starts(i):ends(i),:);%what happens in the block
    event_type_block=events_block.(4);
    b_trials.e(i)=length(find(event_type_block=='enter_explore'));   
    b_trials.r(i)=length(find(event_type_block=='enter_run'));
    b_trials.f(i)=length(find(event_type_block=='enter_feed'));
    b_trials.s(i)=length(find(event_type_block=='enter_social'));
    b_trials.w(i)=length(find(event_type_block=='enter_drink'));
    b_trials.total(i)=b_trials.e(i)+b_trials.r(i)+b_trials.f(i)+b_trials.s(i)+b_trials.w(i);
    
    %high res analyses

    b_trials_explore=find(event_type_block=='enter_explore');
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
        for e=2:size(b_trials_total,1)
            if find(b_trials_explore==b_trials_total(e)) & find(b_trials_run==b_trials_total(e-1)) | find(b_trials_food==b_trials_total(e-1))| find(b_trials_social==b_trials_total(e-1))| find(b_trials_water==b_trials_total(e-1)) %if other trial follows explore
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_run==b_trials_total(e)) & find(b_trials_explore==b_trials_total(e-1)) | find(b_trials_food==b_trials_total(e-1))| find(b_trials_social==b_trials_total(e-1))| find(b_trials_water==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;    
            elseif find(b_trials_food==b_trials_total(e)) & find(b_trials_run==b_trials_total(e-1)) | find(b_trials_explore==b_trials_total(e-1))| find(b_trials_social==b_trials_total(e-1))| find(b_trials_water==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;
            elseif find(b_trials_social==b_trials_total(e)) & find(b_trials_run==b_trials_total(e-1)) | find(b_trials_food==b_trials_total(e-1))| find(b_trials_explore==b_trials_total(e-1))| find(b_trials_water==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;    
            elseif find(b_trials_water==b_trials_total(e)) & find(b_trials_run==b_trials_total(e-1)) | find(b_trials_food==b_trials_total(e-1))| find(b_trials_social==b_trials_total(e-1))| find(b_trials_explore==b_trials_total(e-1)) 
                j=j+1;
                ex_ind=ex_ind+1;
                exploit(ex_ind)=0;    
            else %it was the same then
                exploit(ex_ind)=exploit(ex_ind)+1;
            end
            b_trial_dur=[b_trial_dur; milliseconds(x(b_trials_total_end(e))-x(b_trials_total(e)))];
            b_choice_dur=[b_choice_dur; milliseconds(x(b_trials_total(e))-x(b_trials_total_end(e-1)))];
        end
        b_trials.mean_choice_dur(i)=nanmean(b_choice_dur);
        b_trials.mean_trial_dur(i)=nanmean(b_trial_dur);
    else
        j=NaN;
        exploit=NaN;
        b_trials.mean_trial_dur(i)=NaN;
        b_trials.mean_choice_dur(i)=NaN;
    end
    b_trials(a).switch_index(i)=j/size(b_trials_total,1);
    b_trials(a).exploit_index_max(i)=max(exploit);
    b_trials(a).exploit_index_mean(i)=mean(exploit);
    b_trials(a).exploit_index_min(i)=min(exploit);
    b_trials(a).start_time(i)=x(starts(i));
    %b_trials(a).revolutions(i)=sum(revolutions(starts(i):ends(i)));
end


f2 = figure;
figure(f2);
xx=x(ends);
subplot(size(animals,1),1,a),plot(xx,b_trials.total,'ko','MarkerFaceColor','k');hold on
plot(xx,b_trials(a).e,'co','MarkerFaceColor','c');
plot(xx,b_trials(a).r,'yo','MarkerFaceColor','y');
plot(xx,b_trials(a).f,'co','MarkerFaceColor','c');
plot(xx,b_trials(a).s,'yo','MarkerFaceColor','y');
plot(xx,b_trials(a).w,'co','MarkerFaceColor','c');


%% switch index


Entries = event_type;
A = find(event_type == 'leave_nest');
B = find(event_type == 'block_start');
C= find(event_type == 'retrieve_pellet');
D= find(event_type == 'drink');
E= find(event_type == 'exit_drink');
F= find(event_type == 'exit_explore');
G= find(event_type == 'exit_feed');
H= find(event_type == 'exit_run');
I= find(event_type == 'exit_social');
J= find(event_type == 'frame');
K= find(event_type == 'block_end');
L= find(event_type == 'run');
remove_rows = cat(1, A, B, C, D, E, F, G, H , I, J, K, L);
Entries (remove_rows,:) = [];


% switch index whole 
clear i
switch_counter = 0;
for i = 1:size(Entries,1)-1
    AA = Entries(i);
    BB = Entries(i+1);
    if AA ~= BB
        switch_counter = switch_counter+1;
    end     
end
switch_index = switch_counter/size(Entries,1);





%%
% switch index per block: nr of switches/total nr of entries in that block
% duration and nr of exploits (repeated entering and exploiting of one pod)
% latencies to consume
% time spent in pods
%


