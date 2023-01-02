%Drivemaze etho and ca imaging plotter
%2022March 
clearvars -except events_averg ci_data events;close all;
cellevents=events;
%animals=table(['94331472'; '328340178184'; '335490249236']);%
animals='3354903011';
'3354903011 grin5 female';
f1=figure;
%f2=figure;
%f3=figure;
% f4=figure;
% f5=figure;
% for 
a=1;%:3;
clearvars -except a f1 f2 f3 f4 f5 animals binned b_trials b_dur binned_exp
% animal=cellstr(table2cell(animals(a,1)));
animal=animals;
% filename=['C:\Data\Drivemaze\Drivemaze_imaging_grin1' animal{1, 1} '_events.csv'];
filename=['E:\Exp_2_DM_exploration\PFC-LH\g5\' animal '_events.csv'];

%% Import csv
opts = delimitedTextImportOptions("NumVariables", 5);
opts.DataLines = [1, Inf];
opts.Delimiter = ",";
opts.VariableNames = ["Date_Time", "amount_consumed", "latency_to_consumption", "Type", "Frames"];
opts.VariableTypes = ["string", "double", "double", "categorical", "double"];
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
chop_from=datetime('2022-09-08 17:28:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
chop_to=datetime('2022-09-08 19:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');

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
f_num=events.(3);
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


%% imaging data
%select sessions with analysed raw data
topLevelFolder = "W:\BETA-NeuroSciences-Hypo\imaging\Data\miniscope_PFC_LH\g5_female_2022\g5_\2022_09_08";
% Get a list of all files and folders in this folder.
files = dir(topLevelFolder);
% Get a logical vector that tells which is a directory.
dirFlags = [files.isdir];
% Extract only those that are directories.
subFolders = files(dirFlags); % A structure with extra info.
% Get only the folder names into a cell array.
subFolderNames = {subFolders(3:end).name};
subFolderNames =subFolderNames';
folders = cell2table(subFolderNames);
imaging_datasets=folders(1:17,1);
imaging_sessions=1:17;
x_imaging_cum=[];%cumulators for averaging
raw_f_cum=[];
for s=1:length(imaging_sessions)
    stack=cellstr(table2cell(imaging_datasets(s,1)));
    stack_path=['W:\BETA-NeuroSciences-Hypo\imaging\Data\miniscope_PFC_LH\g5_female_2022\g5_\2022_09_08\' stack{1, 1} '\My_V4_Miniscope\Results.csv'];
    opts = delimitedTextImportOptions("NumVariables", 93);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    cells=23;
    VariableNames = ["VarName1"];
    VariableTypes = ["double"];
    for i = 1:cells
        newVariableNames= ["Area"+ i , "Mean"+ i, "Min"+ i, "Max"+ i];
        B = convertStringsToChars(newVariableNames);
        VariableNames = [VariableNames, B];
    
        newVariableTypes = ["double", "double", "double", "double"];
        C = convertStringsToChars(newVariableTypes);
        VariableTypes = [VariableTypes, C];
    end    
    opts.VariableNames = VariableNames;
    opts.VariableTypes = VariableTypes;
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    raw_f = readtable(stack_path, opts);
    raw_f = table2array(raw_f);
    clear opts
    raw_f=raw_f(:,3:4:end);
    cell_n=size(raw_f,2)   
    %good_cells = [1:23];%[1 2 3 4 5 7 8 9 10 11 12 13 15 16 17 18 19 20 21 22];
    good_cells = [1 4 14 9 15 23];
    %%make an x-axis
    f_dur_adhoc=milliseconds(x(imaging_ends(imaging_sessions(s)))-x(imaging_starts(imaging_sessions(s))))/size(raw_f,1)
    x_imaging=[x(imaging_starts(imaging_sessions(s)))];
    %for individual sessions
%     f_dur_adhoc=(milliseconds(x(imaging_ends(1,1))-x(imaging_starts(1,1)))/size(raw_f,1))
%     x_imaging=[x(imaging_starts(1,1))];
% 
%     for i=1:size(raw_f,1)-2
%         x_imaging=[x_imaging ; x(imaging_starts(1,1))+duration(0,0,0,floor(i*f_dur_adhoc))];
%     end
%     x_imaging=[x_imaging ; x(imaging_ends(1,1))];


    for i=1:size(raw_f,1)-2
        x_imaging=[x_imaging ; x(imaging_starts(imaging_sessions(s)))+duration(0,0,0,floor(i*f_dur_adhoc))];
    end
    x_imaging=[x_imaging ; x(imaging_ends(imaging_sessions(s)))];
    
    %semiautomatic cleaning
%     cut_s=1;cut_e=5;replacement=ones(cut_e-cut_s+1,cell_n).*raw_f(cut_e+1,:);raw_f(cut_s:cut_e,:)=replacement;
%     if s==1 
%         cut_s=4018;cut_e=4031;replacement=ones(cut_e-cut_s+1,cell_n).*raw_f(cut_s-1,:);raw_f(cut_s:cut_e,:)=replacement;
%         cut_s=12434;cut_e=12451;replacement=ones(cut_e-cut_s+1,cell_n).*raw_f(cut_s-1,:);raw_f(cut_s:cut_e,:)=replacement;
%     elseif s==2
%         cut_s=8049;cut_e=8205;replacement=ones(cut_e-cut_s+1,cell_n).*raw_f(cut_s-1,:);raw_f(cut_s:cut_e,:)=replacement;
%     end
    medians = [];
    for ROI = 1:size(raw_f,2)
        cell_n = [];
        cell_n = raw_f(:,ROI); 
        medians_cell = movmedian(cell_n, 1200);%window size in frames %120s
        medians = [medians, medians_cell];     
    end
    df_f_trace = [];
    for cell_median = 1: size(medians,2) % Loop df/f
       v = [];
       d = [];
       v = medians(:,cell_median);%select column from medians
       d = raw_f(:,cell_median);%select column from raw_f_res
       for f = 1: size(d,1)
           a = [];
           a =  v(f) ;
           b = (d(f) - 0.7*a)/0.7*a;
           df_f_trace(f,cell_median) = b;
       end
    end
    z_scored_raw= zscore(df_f_trace,0,1);% zscored using sample sd


    %z_scored_raw=zscore(raw_f,1,1);% z-scoring
    a=1;
    %plotting only good cells
    j=0;
    for i=good_cells
        subplot(1,1,a),plot(x_imaging,smooth(z_scored_raw(:,i)/7-j,3),'Color',[i/size(raw_f,2) 1-i/size(raw_f,2) i/size(raw_f,2)],'LineWidth',1);hold on;
        j=j+1;
    end

    

    %plotting all cells
    %for i=1:size(raw_f,2)
    %    subplot(1,1,a),plot(x_imaging,smooth(z_scored_raw(:,i)/7-i,30),'Color',[i/size(raw_f,2) 1-i/size(raw_f,2) i/size(raw_f,2)]);hold on;
    %end
    % a=2;
    % subplot(1,1,a),imagesc(x_imaging,[1:size(raw_f,2)],raw_f');hold on;
    % figure;
    % for i=1:size(raw_f,2)
    %     plot(x_imaging,raw_f(:,i)-mean(raw_f(:,i))-i*10);hold on;
    % end
    x_imaging_cum=[x_imaging_cum;x_imaging];
    raw_f_cum=[raw_f_cum;raw_f];
end

f1;hold on; p1 = plot([(x_imaging_cum(10603,1)) (x_imaging_cum(10603+600,1))], [1.0006 1.0006],"k");
p2=plot([(x_imaging_cum(10603,1)) (x_imaging_cum(10603,1))],[1.0006 1.2863],"k");p1.LineWidth = 2; p2.LineWidth = 2; 
%axes 2s.d. 1 min 2sd= 0.2857 (2/7)
%text((x_imaging_cum(10603,1)),1.0006,'60s');text((x_imaging_cum(10604,1)),1.15,'2 s.d.');


