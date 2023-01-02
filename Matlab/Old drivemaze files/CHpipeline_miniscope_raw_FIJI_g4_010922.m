clearvars -except ci_data; close all; 
set(0,'defaultAxesFontSize',18);
%% import data
topLevelFolder = "E:\Exp_2_DM_exploration\PFC-LH\g4\2022_09_01"; % or whatever, such as 'C:\Users\John\Documents\MATLAB\work'
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
ko=0;
imaging_datasets = folders(ko+1:7,1);
 %to add to linear imaging starts index, if imaging files start at 20, ko=19
bad=0; %bad imaging dataset to be excluded later
%imaging_datasets=table(['11_25_35'; '11_34_12']);%;'14_48_41';'14_56_33']);%;'';'';''
ed_averg=[];ef_averg=[];es_averg=[];ee_averg=[];er_averg=[];drink_averg=[];eat_averg=[];blockend_averg=[];soc_averg=[];
exd_averg=[];exf_averg=[];exs_averg=[];exe_averg=[];exr_averg=[];imstop_averg=[];run_averg=[];soc_averg_m=[];
ed_averg_bsl=[];ef_averg_bsl=[];es_averg_bsl=[];ee_averg_bsl=[];er_averg_bsl=[];drink_averg_bsl=[];eat_averg_bsl=[];blockend_averg_bsl=[];soc_averg_bsl=[];
exd_averg_bsl=[];exf_averg_bsl=[];exs_averg_bsl=[];exe_averg_bsl=[];exr_averg_bsl=[];imstop_averg_bsl=[];run_averg_bsl=[];soc_averg_bsl_m=[];

for s=1:size(imaging_datasets,1)%:5%:3%change which imaging sets here
%name paths and import
stack=cellstr(table2cell(imaging_datasets(s,1)));
stack_path=['E:\Exp_2_DM_exploration\PFC-LH\g4\2022_09_01\' stack{1, 1} '\My_V4_Miniscope\Results.csv'];
list_path=['E:\Exp_2_DM_exploration\PFC-LH\g4\2022_09_01\' stack{1, 1} '\My_V4_Miniscope\timeStamps.csv'];
events_path=['E:\Exp_2_DM_exploration\PFC-LH\g4\335490167178_events.csv'];
frame_list=import_framelist(list_path);
events=import_events(events_path);
raw_f_res=import_raw(stack_path); %raw data from results file
%df/f
medians = [];
    for ROI = 1:size(raw_f_res,2)
        cell_n = [];
        cell_n = raw_f_res(:,ROI); 
        medians_cell = movmedian(cell_n, 1200);%window size in frames %120s
        medians = [medians, medians_cell];     
    end
    df_f_trace = [];
    for cell_median = 1: size(medians,2) % Loop df/f
       v = [];
       d = [];
       v = medians(:,cell_median);%select column from medians
       d = raw_f_res(:,cell_median);%select column from raw_f_res
       for f = 1: size(d,1)
           a = [];
           a =  v(f) ;
           b = (d(f) - 0.7*a)/0.7*a;
           df_f_trace(f,cell_median) = b;
       end
    end
raw_f= zscore(df_f_trace,0,1);% zscored using sample sd


%% chop
chop_from=datetime('2022-09-01 10:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%chop the day
chop_to=datetime('2022-09-01 20:00:00.000','InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];


events=sortrows(events,1); %spliced in some events from accidental save in different subject
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS'); %make the x axis
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


%% refine events
%chop sessions %maybe do this automatic?
chop_from = x(imaging_starts(s+ko));
chop_to=x(imaging_ends(s+ko));


x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(isnat(x)==1),:)=[];
x(find(isnat(x)==1))=[];
event_type=events.(4);

%% find dropped frames
frame=frame_list.(1);
time_stamp=frame_list.(2);
%figure;subplot(2,1,1),plot(frame(2:end),diff(time_stamp),'ko-');
%title('before');

fd=median(diff(time_stamp));
fake_time_stamps=time_stamp;
missed_frames=[];

miss=find(diff(fake_time_stamps)>1.5*fd,1);%checks if there are missing frames
while isnan(miss)==0
    f=miss+1;
    missed_frames=[missed_frames f];
    new_stamp=fake_time_stamps(f-1)+fd;
    fake_time_stamps=[fake_time_stamps(1:f-1); new_stamp; fake_time_stamps(f:end)];
    miss=find(diff(fake_time_stamps)>1.5*fd,1);
end
%subplot(2,1,2),plot(diff(fake_time_stamps),'ko-'); %plots original frames
%vs corrected ones
%title('after');
missed_frames

%% correct event list for dropped frames
event_frames=events.(5);
new_event_frames=events.(5);
for i=1:length(missed_frames)
    to_shift=find(event_frames>missed_frames(i));
    new_event_frames(to_shift)=new_event_frames(to_shift)-1;
end
events.(5)=new_event_frames;

%% imaging data
%figure;stackedplot(raw_f);
size(raw_f);
win=30;
m=size(imaging_datasets,1);
n=8;

for p=1:15;
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

    trigs=new_event_frames(find(event_type==t));
    if size(trigs,1)==0
        p=p+1;
    else
        for i=1:size(trigs,1)
            trig=trigs(i);
            if trig>win+1 & p==1 
                snips(:,:,i)=raw_f(trig-win:trig,:); %entry drink
                entry_d_bsl = mean(snips(1:10,:,i),1);
                snips_bslined(:,:,i)=snips(:,:,i)-entry_d_bsl;
                ed_averg=cat(3,ed_averg,snips(:,:,i));
                ed_averg_bsl=cat(3,ed_averg_bsl,snips_bslined(:,:,i));
            elseif trig>win+1 & p==2
                snips(:,:,i)=raw_f(trig-win:trig,:); % entry food
                entry_f_bsl = mean(snips(1:10,:,i),1);
                snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
                ef_averg=cat(3,ef_averg,snips(:,:,i));
                ef_averg_bsl=cat(3,ef_averg_bsl,snips_bslined(:,:,i));
            elseif trig>win+1 & p==3
                snips(:,:,i)=raw_f(trig-win:trig,:); % entry social
                entry_s_bsl = mean(snips(1:10,:,i),1);
                snips_bslined(:,:,i)=snips(:,:,i)-entry_s_bsl;
                es_averg=cat(3,es_averg,snips(:,:,i));
                es_averg_bsl=cat(3,es_averg_bsl,snips_bslined(:,:,i));

                snips2(:,:,i)=raw_f(trig:trig+win,:); %social interact
                snips_bslined2(:,:,i)=snips(:,:,i)-entry_s_bsl;
                soc_averg=cat(3,soc_averg,snips2(:,:,i));
                soc_averg_bsl=cat(3,soc_averg_bsl,snips_bslined2(:,:,i));

            elseif trig>win+1 & p==4
                snips(:,:,i)=raw_f(trig-win:trig,:); %entry explore
                entry_e_bsl = mean(snips(1:10,:,i),1);
                snips_bslined(:,:,i)=snips(:,:,i)-entry_e_bsl;
                ee_averg=cat(3,ee_averg,snips(:,:,i));
                ee_averg_bsl=cat(3,ee_averg_bsl,snips_bslined(:,:,i));
            elseif trig>win+1 & p==5
                snips(:,:,i)=raw_f(trig-win:trig,:); %entry run
                entry_r_bsl = mean(snips(1:10,:,i),1);
                snips_bslined(:,:,i)=snips(:,:,i)-entry_r_bsl;
                er_averg=cat(3,er_averg,snips(:,:,i));
                er_averg_bsl=cat(3,er_averg_bsl,snips_bslined(:,:,i));
            elseif trig>win+1 & p==6 %drink
                snips(:,:,i)=raw_f(trig-win:trig+win,:);
                snips_bslined(:,:,i)=snips(:,:,i)-entry_d_bsl;
                drink_averg=cat(3,drink_averg,snips(:,:,i));
                drink_averg_bsl=cat(3,drink_averg_bsl,snips_bslined(:,:,i));
            elseif trig>win+1 & p==7 %pellet
                snips(:,:,i)=raw_f(trig-win:trig+win,:);
                snips_bslined(:,:,i)=snips(:,:,i)-entry_f_bsl;
                eat_averg=cat(3,eat_averg,snips(:,:,i));
                eat_averg_bsl=cat(3,eat_averg_bsl,snips_bslined(:,:,i));
            elseif trig>win+1 & p==8 %blockend
                snips(:,:,i)=raw_f(trig-win:trig+win,:);
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
            ed_averg_m= nanmean(ed_averg,3);
            ef_averg_m= nanmean(ef_averg,3);
            es_averg_m= nanmean(es_averg,3);
            soc_averg_m= nanmean(soc_averg,3);
            ee_averg_m= nanmean(ee_averg,3);
            er_averg_m= nanmean(er_averg,3);
            drink_averg_m= nanmean(drink_averg,3);
            eat_averg_m= nanmean(eat_averg,3);
            blockend_averg_m= nanmean(blockend_averg,3);
            exd_averg_m= nanmean(exd_averg,3);
            exf_asverg_m= nanmean(exf_averg,3);
            exs_averg_m= nanmean(exs_averg,3);
            exe_averg_m= nanmean(exe_averg,3);
            exr_averg_m= nanmean(exr_averg,3);
            imstop_averg_m= nanmean(imstop_averg,3);
            run_averg_m= nanmean(run_averg,3);

            ed_averg_bsl_m= nanmean(ed_averg_bsl,3);
            ef_averg_bsl_m= nanmean(ef_averg_bsl,3);
            es_averg_bsl_m= nanmean(es_averg_bsl,3);
            soc_averg_bsl_m= nanmean(soc_averg_bsl,3);
            ee_averg_bsl_m= nanmean(ee_averg_bsl,3);
            er_averg_bsl_m= nanmean(er_averg_bsl,3);
            drink_averg_bsl_m= nanmean(drink_averg_bsl,3);
            eat_averg_bsl_m= nanmean(eat_averg_bsl,3);
            blockend_averg_bsl_m= nanmean(blockend_averg_bsl,3);
            exd_averg_bsl_m= nanmean(exd_averg_bsl,3);
            exf_asverg_m= nanmean(exf_averg_bsl,3);
            exs_averg_bsl_m= nanmean(exs_averg_bsl,3);
            exe_averg_bsl_m= nanmean(exe_averg_bsl,3);
            exr_averg_bsl_m= nanmean(exr_averg_bsl,3);
            imstop_averg_bsl_m= nanmean(imstop_averg_bsl,3);
            run_averg_bsl_m= nanmean(run_averg_bsl,3);
        end
%         response_mean=nanmean(snips_bslined,3);
%        response_mean=nanmean(snips,3); %makes individual plots per session according to nr in imaging sessions
%        figure(f3);hold on;
%        subplot(m,n,p+s*7-7),imagesc([-win:win],[1:size(response_mean,2)],response_mean');
%        title(['\fontsize{10}',tt,' trials:',num2str(size(snips,3))]);
%        session(s).event(p).snips =snips; %something like this to save and append over
% %         sessions
        
    end
ci_data.g4.Sep01(s).session.raw_f_res = raw_f_res;
ci_data.g4.Sep01(s).session.events = events;
ci_data.g4.Sep01(s).session.missed_frames = missed_frames;    
end




end

%% averg plots
lim = [-3 3];
f4=figure;% averages across all sessions
figure(f4);hold on;
        subplot(1,5,1),imagesc([-win:0],[1:size(ed_averg_m,2)],ed_averg_m');
        title(['\fontsize{10}enter drink',' trials:',num2str(size(ed_averg,3))]);clim(lim)

        subplot(1,5,2),imagesc([-win:0],[1:size(ef_averg_m,2)],ef_averg_m');
        title(['\fontsize{10}enter feed',' trials:',num2str(size(ef_averg,3))]);clim(lim)

        subplot(1,5,3),imagesc([-win:0],[1:size(es_averg_m,2)],es_averg_m');
        title(['\fontsize{10}enter social',' trials:',num2str(size(es_averg,3))]);clim(lim)

        subplot(1,5,4),imagesc([-win:0],[1:size(ee_averg_m,2)],ee_averg_m');
        title(['\fontsize{10}enter explore',' trials:',num2str(size(ee_averg,3))]);clim(lim)

        subplot(1,5,5),imagesc([-win:0],[1:size(er_averg_m,2)],er_averg_m');
        title(['\fontsize{10}enter run',' trials:',num2str(size(er_averg,3))]);clim(lim)
        h = axes(f4,'visible','off');
        c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);
        caxis(h,lim);

f6=figure;% averages across all sessions
figure(f6);hold on;
        subplot(1,5,1),imagesc([-win:win],[1:size(exd_averg_m,2)],exd_averg_m');
        title(['\fontsize{10}exit drink',' trials:',num2str(size(exd_averg,3))]);clim(lim)

        subplot(1,5,2),imagesc([-win:win],[1:size(exf_asverg_m,2)],exf_asverg_m');
        title(['\fontsize{10}exit feed',' trials:',num2str(size(exf_averg,3))]);clim(lim)

        subplot(1,5,3),imagesc([-win:win],[1:size(exs_averg_m,2)],exs_averg_m');
        title(['\fontsize{10}exit social',' trials:',num2str(size(exs_averg,3))]);clim(lim)

        subplot(1,5,4),imagesc([-win:win],[1:size(exe_averg_m,2)],exe_averg_m');
        title(['\fontsize{10}exit explore',' trials:',num2str(size(exe_averg,3))]);clim(lim)

        subplot(1,5,5),imagesc([-win:win],[1:size(exr_averg_m,2)],exr_averg_m');
        title(['\fontsize{10}exit run',' trials:',num2str(size(exr_averg,3))]);clim(lim)
        h = axes(f6,'visible','off');
        c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);
        caxis(h,lim);        

lim = [-3 3];
f5=figure;% averages across all sessions
figure(f5);hold on;
        subplot(1,6,1),imagesc([-win:win],[1:size(drink_averg_m,2)],drink_averg_m');
        title(['\fontsize{10}drink',' trials:',num2str(size(drink_averg,3))]);clim(lim)

        subplot(1,6,2),imagesc([-win:win],[1:size(eat_averg_m,2)],eat_averg_m');
        title(['\fontsize{10}retrieve pellet',' trials:',num2str(size(eat_averg,3))]);clim(lim)

        subplot(1,6,3),imagesc([0:win],[1:size(soc_averg_m,2)],soc_averg_m');
        title(['\fontsize{10}social interaction',' trials:',num2str(size(soc_averg,3))]);clim(lim)
       
        subplot(1,6,4),imagesc([-win:win],[1:size(run_averg_m,2)],run_averg_m');
        title(['\fontsize{10}run',' trials:',num2str(size(run_averg,3))]);clim(lim)

        subplot(1,6,5),imagesc([-win:win],[1:size(blockend_averg_m,2)],blockend_averg_m');
        title(['\fontsize{10}block end',' trials:',num2str(size(blockend_averg,3))]);clim(lim)

        subplot(1,6,6),imagesc([-win:0],[1:size(imstop_averg_m,2)],imstop_averg_m');
        title(['\fontsize{10}imaging end',' trials:',num2str(size(imstop_averg,3))]);clim(lim)
        h = axes(f5,'visible','off');
        c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);
        caxis(h,lim);        

 %% averg plots baselined
lim = [-3 3];
f7=figure;% averages across all sessions
figure(f7);hold on;
        subplot(1,5,1),imagesc([-win:0],[1:size(ed_averg_bsl_m,2)],ed_averg_bsl_m');
        title(['\fontsize{10}enter drink',' trials:',num2str(size(ed_averg_bsl,3))]);clim(lim)

        subplot(1,5,2),imagesc([-win:0],[1:size(ef_averg_bsl_m,2)],ef_averg_bsl_m');
        title(['\fontsize{10}enter feed',' trials:',num2str(size(ef_averg_bsl,3))]);clim(lim)

        subplot(1,5,3),imagesc([-win:0],[1:size(es_averg_bsl_m,2)],es_averg_bsl_m');
        title(['\fontsize{10}enter social',' trials:',num2str(size(es_averg_bsl,3))]);clim(lim)

        subplot(1,5,4),imagesc([-win:0],[1:size(ee_averg_bsl_m,2)],ee_averg_bsl_m');
        title(['\fontsize{10}enter explore',' trials:',num2str(size(ee_averg_bsl,3))]);clim(lim)

        subplot(1,5,5),imagesc([-win:0],[1:size(er_averg_bsl_m,2)],er_averg_bsl_m');
        title(['\fontsize{10}enter run',' trials:',num2str(size(er_averg_bsl,3))]);clim(lim)
        h = axes(f7,'visible','off');
        c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);
        caxis(h,lim);

f8=figure;% averages across all sessions
figure(f8);hold on;
        subplot(1,5,1),imagesc([-win:win],[1:size(exd_averg_bsl_m,2)],exd_averg_bsl_m');
        title(['\fontsize{10}exit drink',' trials:',num2str(size(exd_averg_bsl,3))]);clim(lim)

        subplot(1,5,2),imagesc([-win:win],[1:size(exf_asverg_m,2)],exf_asverg_m');
        title(['\fontsize{10}exit feed',' trials:',num2str(size(exf_averg_bsl,3))]);clim(lim)

        subplot(1,5,3),imagesc([-win:win],[1:size(exs_averg_bsl_m,2)],exs_averg_bsl_m');
        title(['\fontsize{10}exit social',' trials:',num2str(size(exs_averg_bsl,3))]);clim(lim)

        subplot(1,5,4),imagesc([-win:win],[1:size(exe_averg_bsl_m,2)],exe_averg_bsl_m');
        title(['\fontsize{10}exit explore',' trials:',num2str(size(exe_averg_bsl,3))]);clim(lim)

        subplot(1,5,5),imagesc([-win:win],[1:size(exr_averg_bsl_m,2)],exr_averg_bsl_m');
        title(['\fontsize{10}exit run',' trials:',num2str(size(exr_averg_bsl,3))]);clim(lim)
        h = axes(f8,'visible','off');
        c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);
        caxis(h,lim);        

lim = [-3 3];
f9=figure;% averages across all sessions
figure(f9);hold on;
        subplot(1,6,1),imagesc([-win:win],[1:size(drink_averg_bsl_m,2)],drink_averg_bsl_m');
        title(['\fontsize{10}drink',' trials:',num2str(size(drink_averg_bsl,3))]);clim(lim)

        subplot(1,6,2),imagesc([-win:win],[1:size(eat_averg_bsl_m,2)],eat_averg_bsl_m');
        title(['\fontsize{10}retrieve pellet',' trials:',num2str(size(eat_averg_bsl,3))]);clim(lim)

        subplot(1,6,3),imagesc([0:win],[1:size(soc_averg_bsl_m,2)],soc_averg_bsl_m');
        title(['\fontsize{10}social interaction',' trials:',num2str(size(soc_averg_bsl,3))]);clim(lim)
       
        subplot(1,6,4),imagesc([-win:win],[1:size(run_averg_bsl_m,2)],run_averg_bsl_m');
        title(['\fontsize{10}run',' trials:',num2str(size(run_averg_bsl,3))]);clim(lim)

        subplot(1,6,5),imagesc([-win:win],[1:size(blockend_averg_bsl_m,2)],blockend_averg_bsl_m');
        title(['\fontsize{10}block end',' trials:',num2str(size(blockend_averg_bsl,3))]);clim(lim)

        subplot(1,6,6),imagesc([-win:0],[1:size(imstop_averg_bsl_m,2)],imstop_averg_bsl_m');
        title(['\fontsize{10}imaging end',' trials:',num2str(size(imstop_averg_bsl,3))]);clim(lim)
        h = axes(f9,'visible','off');
        c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);
        caxis(h,lim);        


function raw_f_res=import_raw(filename)
    %% Import csv
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
    Results = readtable(filename, opts);
    Results = table2array(Results);
    raw_f_res=Results(:,3:4:end);
%     medians = [];
%     for ROI = 1:size(raw_f_res,2)
%         cell_n = [];
%         cell_n = raw_f_res(:,ROI); 
%         medians_cell = movmedian(cell_n, 1200);%window size in frames %120s
%         medians = [medians, medians_cell];     
%     end
%     df_f_trace = [];
%     for cell_median = 1: size(medians,2) % Loop df/f
%        v = [];
%        d = [];
%        v = medians(:,cell_median);%select column from medians
%        d = raw_f_res(:,cell_median);%select column from raw_f
%        for f = 1: size(d,1)
%            a = [];
%            a =  v(f) ;
%            b = (d(f) - 0.7*a)/0.7*a;
%            df_f_trace(f,cell_median) = b;
%        end
%     end
%     
%     raw_f= zscore(df_f_trace,0,1);% zscored using sample sd
    clear opts
    %cells=unique(Results(:,2));
    %for i=1:length(cells)
      %  raw_f(:,i)=Results(find(Results(:,2)==cells(i)),4);
    %end
end
function events=import_events(filename)
    %% Import csv
    opts = delimitedTextImportOptions("NumVariables", 7);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["Date_Time", "amount_consumed", "latency_to_consumption", "Type", "frame", "hardware_time", "experiment"];
    opts.VariableTypes = ["string", "double", "double", "categorical", "double", "double", "double"];
    opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 4], "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    events = readtable(filename, opts);
    clear opts
end
function frame_list=import_framelist(list_path)
    %% Import csv
    opts = delimitedTextImportOptions("NumVariables", 3);
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";
    opts.VariableNames = ["FrameNumber", "TimeStampms", "BufferIndex"];
    opts.VariableTypes = ["double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    frame_list = readtable(list_path, opts);
    clear opts
end
