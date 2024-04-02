clear all; close all;


%% get cycle length and F
R(2).F=[];
R(2).L=[];

% import data
for animal=1:2
    if animal==1
        load('3354916686_g10_femP150.mat')
        good_cells=[1 2 4 5 8 10 11 13];
        k=1;
    elseif animal==2
        load('3354909574_g12_femP150.mat')
        good_cells=[1 2 3 4 5 6 7 8 9];
        k=1;
    elseif animal==3
        load('1251667485196_g15_mP200.mat')
        good_cells=[1:17]; 
        k=1;
    end
    for exp=2
        if exp==1
            analyse=string('habituation ad lib')
        elseif exp==2
            analyse=string('bsl ad lib')
        elseif exp==3
            analyse=string('swapFoodDrink ad lib lockin')
        end
        data_in=meta(find(cat(1,meta.exp_note)==analyse)).imaging;

        for s=1:size(data_in,2)
            clearvars -except data_in meta analyse R k s exp good_cells animal
            events=data_in(s).events;
            event_type=events.(4);
            event_frames=events.(5);    
            x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
            diff_f=data_in(s).np_f;%data_in(s).cell_f-
            diff_f=diff_f(:,good_cells);
            diff_f_z=zscore(diff_f,1,1);

           %find L and F in this imaging block
            t='exit_feed';
            consum='retrieve_pellet';
            pre='enter_feed';        
            all_t=find(event_type==t);
            all_consum=find(event_type==consum);
            all_pre=find(event_type==pre);
            if size(all_t)==size(all_pre)
                'check ok'
            else
                'event list error'
                stop
            end
            amount_consumed=events.(2);
            del=find(amount_consumed(all_t)<3);
            all_t(del)=[];
            all_pre(del)=[];
            R(animal).L=[R(animal).L;amount_consumed(all_t)];  
            if isempty(all_t)
                'no trials'
            else
                for i=1:size(all_t,1)
                    t1=all_consum(find(all_consum>all_pre(i),1,'first'));
                    t2=all_consum(find(all_consum<all_t(i),1,'last'));
                    R(animal).F(k,:)=mean(diff_f_z(event_frames(t1):event_frames(t2),:));
                    k=k+1;
                end
            end
        end
    end
end

%% plot feeding F v cycle length
f1=figure;
f2=figure;
f3=figure;
for animal=1:2
    L=R(animal).L;
    F=R(animal).F;
    [L inds]=sort(L,'ascend');
    F=F(inds,:);
    figure(f1);
    for cell=1:size(F,2)
        plot(L,F(:,cell),'Color',[cell/size(F,2) 1-cell/size(F,2) cell/size(F,2)],'LineWidth',3/sqrt(cell));hold on;
    end

figure;histogram(L);

figure(f2);
for cell=1:size(F,2)
    c=[cell/size(F,2) 1-cell/size(F,2) cell/size(F,2)];
    scatter(L,F(:,cell),[],c,'filled');hold on;
end


figure(f3);
binF=[];
binL=[];
bins=[2 3;
      4 8;
      9 12;
      13 max(L)];
for bin=1:size(bins,1)
    clear inds
    inds=find(L>bins(bin,1)-1 & L<bins(bin,2)+1);
    binL(bin)=mean(L(inds));
    binF(bin,:)=mean(F(inds,:));   
end
for cell=1:size(F,2)
    plot(binL,binF(:,cell),'Color',[cell/size(F,2) 1-cell/size(F,2) cell/size(F,2)],'LineWidth',3/sqrt(cell));hold on;
end
    plot(binL,mean(binF(:,cell),2),'ko-','LineWidth',4);
end