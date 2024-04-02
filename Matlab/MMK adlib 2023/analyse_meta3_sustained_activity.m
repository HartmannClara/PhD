clear all;
close all;
load('collated_meal_single_snips_bothanimals.mat');
f2=figure;
f3=figure;
figure(f2);
n=2
anwin=[win+100:win+meal_win];

%% clean up
% del_blocks=
del_cells1=[3 6 7 8]+2
for exp=1:3
    for p=3
        F_meals1(exp).all_snips(p).snips(:,del_cells1,:)=[];
        p=4;
        data_out1(exp).all_snips(p).snips(:,del_cells1,:)=[];
    end
end

%%
f4=figure;
for exp=1%:2
    p=3
    on=2;
    off=on;
    figure(f2);
    on_examples=cat(2,F_meals1(exp).all_snips(p).snips(:,3:end,on),F_meals2(exp).all_snips(p).snips(:,3:end,on));
    response_mean1=nanmean(F_meals1(exp).all_snips(p).snips(:,3:end,:),3);
    response_mean2=nanmean(F_meals2(exp).all_snips(p).snips(:,3:end,:),3);
    response_mean=cat(2,response_mean1,response_mean2);
    subplot(n,3,exp),imagesc([-win:meal_win],[1:size(response_mean,2)],response_mean',[-2 2]);hold on;
    title(['long feeding cycles:',num2str(size(F_meals1(exp).all_snips(p).snips,3)),'/',num2str(size(F_meals2(exp).all_snips(p).snips,3))]);
    plot([0 0],[0 size(response_mean,2)],'k--');

    figure(f3);
%     plot(smooth(mean(response_mean(:,find(mean(response_mean(anwin,:),1)<0.0 & mean(response_mean(anwin,:),1)>-0.1)),2)),'b','LineWidth',1);hold on;
    up_cells=find(mean(response_mean(anwin,:),1)>-0.2);
    down_cells=find(mean(response_mean(anwin,:),1)<-0.2);
    subplot(1,2,1),plot(smooth(mean(response_mean(:,up_cells),2),1),'r','LineWidth',2);hold on;
    subplot(1,2,1),plot(smooth(mean(response_mean(:,down_cells),2),1),'b','LineWidth',1);hold on;
    
    legend('sustained','inhibited');
    subplot(1,2,1),plot([win win],[-0.5 0.3],'k--');
    ylim([-0.6 0.9]);
    figure(f4);
    subplot(1,2,1),stackedplot(on_examples);
    
    
end
f5=figure;
for exp=1%:2
    p=4
    figure(f2);
    
    off_examples=cat(2,data_out1(exp).all_snips(p).snips(:,3:end,off),data_out2(exp).all_snips(p).snips(:,3:end,off));
    response_mean1=nanmean(data_out1(exp).all_snips(p).snips(:,3:end,:),3);
    response_mean2=nanmean(data_out2(exp).all_snips(p).snips(:,3:end,:),3);
    response_mean=cat(2,response_mean1,response_mean2);
    subplot(n,3,2),imagesc([-win:meal_win],[1:size(response_mean,2)],response_mean',[-2 2]);hold on;
    title(['end of feeding cycle:',num2str(size(F_meals1(exp).all_snips(p).snips,3)),'/',num2str(size(F_meals2(exp).all_snips(p).snips,3))]);
    plot([0 0],[0 size(response_mean,2)],'k--');

    figure(f3);
%     plot(smooth(mean(response_mean(:,find(mean(response_mean(anwin,:),1)<0.0 & mean(response_mean(anwin,:),1)>-0.1)),2)),'b','LineWidth',1);hold on;
    subplot(1,2,2),plot(smooth(mean(response_mean(:,up_cells),2),1),'r','LineWidth',2);hold on;
    subplot(1,2,2),plot(smooth(mean(response_mean(:,down_cells),2),1),'b','LineWidth',1);hold on;
    
    legend('sustained','inhibited');
    subplot(1,2,2),plot([meal_win meal_win],[-0.5 0.3],'k--');
    ylim([-0.6 0.9]);
    
    figure(f4);
    subplot(1,2,2),stackedplot(off_examples);
    
    
end

