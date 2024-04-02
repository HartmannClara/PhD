clear all; 
close all;

%% import data
%animals=table({'g10_femP150'; 'g12_femP150'; 'g15_mP200'});%
%animal_tags=table({'3354916686'; '3354909574'; '1251667485196'});%

for animal=1
    if animal==1
        load('3354916686_g10_femP150.mat')
        good_cells=[1 2 4 5 8 10 11 13];
    elseif animal==2
        load('3354909574_g12_femP150.mat')
        good_cells=[1 2 3 4 5 6 7 8 9];
    elseif animal==3
        load('1251667485196_g15_mP200.mat')
        good_cells=[1:17];    
    end
end

blockfig=figure;
cols=0;
for i=1:size(meta,2)
    cols=cols+size(cat(3,meta(i).imaging.stack_id),3);
end
if cols>10
    blockfig2=figure;
    cols=ceil(cols/2);
end
col=0; %column to plot into

%% make psth
F_meals=[];
F_singles=[];
for exp=2
    
    if exp==1
        analyse=string('habituation ad lib')
    elseif exp==2
        analyse=string('bsl ad lib')
    elseif exp==3
        analyse=string('swapFoodDrink ad lib lockin')
    end
    data_in=meta(find(cat(1,meta.exp_note)==analyse)).imaging;

    win=100;
    meal_win=600; %2min post win for meals
    data_out(exp).all_snips(4).snips=[];
    F_meals(exp).all_snips(4).snips=[];
    F_singles(exp).all_snips(4).snips=[];
    for s=1:size(data_in,2)
        clearvars -except F_singles F_meals data_in meta analyse data_out s win meal_win data_out exp good_cells col cols blockfig blockfig2
        events=data_in(s).events;
        event_type=events.(4);
        event_frames=events.(5);    
        x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
        diff_f=data_in(s).cell_f-1*data_in(s).np_f;% %0.95 leaves sig component
        temp=[mean(data_in(s).np_f,2) mean(diff_f(:,good_cells),2) diff_f(:,good_cells)];
        diff_f_z=zscore(temp,1,1);
        for c=1:size(diff_f_z,2)
            diff_f_z(:,c)=smooth(diff_f_z(:,c),10);
        end
        
        col=col+1;
        if col>cols
            figure(blockfig2);
            col2=col-cols;
        else
            figure(blockfig);
            col2=col;
        end
        % plot 
        x=datetime(events.(1),'Format','yyyy-MM-dd HH:mm:ss.SSS');
        y=NaN(size(x));
        y(find(event_type=='imaging_start'))=2;
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
        enter_all=x(find(event_type=='enter_explore' | event_type=='enter_run' | event_type=='enter_feed' | event_type=='enter_social' | event_type=='enter_drink'));
        %consume
        y(find(event_type=='retrieve_pellet'))=4;%pellet retrieval
        y(find(event_type=='drink'))=2.5; %
        y(find(event_type=='run'))=3; %wheel entry
        pel=x(find(event_type=='retrieve_pellet'));
        drink=x(find(event_type=='drink'));
        run=x(find(event_type=='run'));
        %exit
        y(find(event_type=='block_end'))=2;
        y(find(event_type=='imaging_stop'))=2.3;
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
        exit_all=x(find(event_type=='exit_explore' | event_type=='exit_run' | event_type=='exit_feed' | event_type=='exit_social' | event_type=='exit_drink'));
        yy=fillmissing(y,'linear');
        ax(col)=subplot(2,cols,col2),plot([x(1) x(end)],[3.5 3.5],'b');hold on%decisionpoint
        ax(col)=subplot(2,cols,col2),plot([x(1) x(end)],[2 2],'b');%nest
        ax(col)=subplot(2,cols,col2),plot(x,yy,'k');hold on
        ax(col)=subplot(2,cols,col2),plot(x,y,'ko');hold on
        if col2==1
            yticks([0.5 1.5 2 2.5 3 3.5 4 4.5 5]);
            yticklabels({'frames','frames','nest','drink','run','decision point','food','social','explore'});
            ylabel('area');
        end
        ylim([-40 5.5]);
        for i=1:1:length(pel)
            ax(col)=subplot(2,cols,col2),plot([pel(i) pel(i)],[3.8 4.2],'y','LineWidth',3);
        end
        for i=1:1:length(drink)
            ax(col)=subplot(2,cols,col2),plot([drink(i) drink(i)],[2.4 2.8],'c','LineWidth',3);
        end
        for i=1:1:length(run)
            ax(col)=subplot(2,cols,col2),plot([run(i) run(i)],[2.8 3.2],'m','LineWidth',3);
        end
        axis 'tight'
        f_dur_adhoc=milliseconds(x(find(event_type=='imaging_stop'))-x(find(event_type=='imaging_start')))/size(diff_f_z,1)
        if f_dur_adhoc>110 | f_dur_adhoc<90
            stoptsop
        end
        x_imaging=[x(find(event_type=='imaging_start'))];
        for i=1:size(diff_f_z,1)-1
            x_imaging=[x_imaging ; x(find(event_type=='imaging_start'))+duration(0,0,0,floor(i*f_dur_adhoc))];
        end
        j=0;
        for i=1:2
            j=j+1;
            ax(col+cols)=subplot(2,cols,col2+cols),plot(x_imaging,smooth(diff_f_z(:,i)/7-j,10),'Color',[i/10 0 0]);hold on;
        end
        for i=3:size(diff_f_z,2)
            j=j+1;
            ax(col+cols)=subplot(2,cols,col2+cols),plot(x_imaging,smooth(diff_f_z(:,i)/7-j,10),'Color',[i/size(diff_f_z,2) 1-i/size(diff_f_z,2) i/size(diff_f_z,2)]);hold on;
        end
        axis 'tight'
        linkaxes([ax(col) ax(col+cols)], 'x');           
        %tsop
       %clip snips here
       for p=1:4 %four diff events
            trigs=[];
            if p<3 %all inclusive event snips
                if p==1
                    t(p).t='enter_feed';% enter_feed
                elseif p==2
                    t(p).t='enter_drink';% enter_drink
                end
                trigs=event_frames(find(event_type==t(p).t));
            end

            if p>2 %selective event snips
                if p==3
                    t(p).t='retrieve_pellet';
                    pre='enter_feed';
                    post='exit_feed';
                elseif p==4
%                     t(p).t='drink';
%                     pre='enter_drink';
                    t(p).t='exit_feed';
                    pre='retrieve_pellet';
                end
                j=1;
                all=find(event_type==t(p).t);
                for i=1:size(all,1)
                    if event_type(all(i)-1)==pre
%                         if seconds(x(all(i))-x(all(i)-1))>0.5*win/10 %non-overlapping
                            trigs(j)=event_frames(all(i));
                            j=j+1;
%                         end
                    end
                end 
            end

            if size(trigs,1)==0
                'no trials'
            else
                if p~=4
                    for i=1:length(trigs)
                        trig=trigs(i);
                        if trig>win+1 & trig<(size(diff_f_z,1)-win-1)
                            snips(:,:,i)=diff_f_z(trig-win:trig+win,:);
    %                         snips_bslined(:,:,i)=snips(:,:,i)-mean(snips(win-10:win,:,i),1);
                            data_out(exp).all_snips(p).snips=cat(3,data_out(exp).all_snips(p).snips, snips(:,:,i));%append over blocks
                        else
                            snips(:,:,i)=NaN(win+1+win,size(diff_f_z,2),1);
                        end
                    end
                end
                
                if p==4
%                     data_out(exp).all_snips(p).snips=[];
                    snips=[];
                    for i=1:length(trigs)
                        trig=trigs(i);
                        if trig>meal_win+1 & trig<(size(diff_f_z,1)-win-1)
                            snips(:,:,i)=diff_f_z(trig-meal_win:trig+win,:);
    %                         snips_bslined(:,:,i)=snips(:,:,i)-mean(snips(win-10:win,:,i),1);
                            data_out(exp).all_snips(p).snips=cat(3,data_out(exp).all_snips(p).snips, snips(:,:,i));%append over blocks
                        else
                            snips(:,:,i)=NaN(meal_win+1+win,size(diff_f_z,2),1);
                        end
                    end
                end
                
                if p==3
                    food_exits=event_frames(find(event_type==post))';
                    if size(trigs,2)<size(food_exits,2)
                        j=1;
                        while j<size(trigs,2)+1
                            if food_exits(j)<trigs(j)
                                food_exits(j)=[];
                            else
                                j=j+1;
                            end
                        end
                    end
                    while size(trigs,2)<size(food_exits,2)
                        food_exits(end)=[];
                    end
                    if size(trigs,2)~=size(food_exits,2)
                        'feeding cycle mismatch in event list'
                        stop
                    end
                    for i=1:length(trigs)%widthwidth
                        trig=trigs(i);
                        if trig>win+1 & trig<(size(diff_f_z,1)-win-1)
                            if trig+meal_win>food_exits(i)
                                meal_snips(:,:,i)=[diff_f_z(trig-win:food_exits(i),:); nan(trig+meal_win-food_exits(i),size(diff_f_z,2))] ;
                            else
                                meal_snips(:,:,i)=diff_f_z(trig-win:trig+meal_win,:);
                            end
                            if find(event_frames==food_exits(i))-find(event_frames==trig)<2
                                F_singles(exp).all_snips(p).snips=cat(3,F_singles(exp).all_snips(p).snips, meal_snips(:,:,i));
                            else
                                F_meals(exp).all_snips(p).snips=cat(3,F_meals(exp).all_snips(p).snips, meal_snips(:,:,i));%append over blocks
                            end
                        else
                            meal_snips(:,:,i)=NaN(win+1+meal_win,size(diff_f_z,2),1);
                        end
                    end
                end
            end        
       end
    end
end

%% plot psth
f1=figure;
figure(f1);
n=4;
for exp=1:2%3
    for p=1:4
        response_mean=nanmean(data_out(exp).all_snips(p).snips,3);
        subplot(n,3,((p*3)-3+exp)),imagesc([-win:win],[1:size(response_mean,2)],response_mean');hold on;
        title([t(p).t ' trials:',num2str(size(data_out(exp).all_snips(p).snips,3))]);
        plot([0 0],[0 size(response_mean,2)],'k--');
    end
    if exp==1
        analyse=string('habituation ad lib')
    elseif exp==2
        analyse=string('bsl ad lib')
    elseif exp==3
        analyse=string('swapFoodDrink ad lib lockin')
    end
    subplot(n,3,((1*3)-3+exp)),text(-20,-5,analyse);
end
f2=figure;
f3=figure;
figure(f2);
n=2
for exp=1%:3%:2
    p=3
    figure(f2);
    response_mean=nanmean(F_meals(exp).all_snips(p).snips,3);
    subplot(n,3,exp),imagesc([-win:meal_win],[1:size(response_mean,2)],response_mean',[-2 2]);hold on;
    title(['long feeding cycles:',num2str(size(F_meals(exp).all_snips(p).snips,3))]);
    plot([0 0],[0 size(response_mean,2)],'k--');
    figure(f3);
    plot(mean(response_mean(:,3:end),2),'k','LineWidth',2);hold on;
    
    figure(f2);
    response_mean=nanmean(F_singles(exp).all_snips(p).snips,3);
    subplot(n,3,exp+3),imagesc([-win:meal_win],[1:size(response_mean,2)],response_mean',[-2 2]);hold on;
    title(['food singles:',num2str(size(F_singles(exp).all_snips(p).snips,3))]);
    plot([0 0],[0 size(response_mean,2)],'k--');
    figure(f3);
    plot(mean(response_mean(:,3:end),2),'r','LineWidth',1);
    legend('long cycles','singles');
    plot([win win],[-0.5 0.3],'k--');
    
end

