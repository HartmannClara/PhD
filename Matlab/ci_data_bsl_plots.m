clearvars -except events_averg ci_data; close all;
set(0,'defaultAxesFontSize',18);

%%
fn=fieldnames(events_averg);
%loop through the fields
for ii=1: numel(fn)
    fn1=fieldnames(events_averg.(fn{ii}));
        %for j=16:numel(fn1)
        events_averg.(fn{ii}).ed_averg_bsl_m = nanmean(events_averg.(fn{ii}).(fn1{16}),3);
        events_averg.(fn{ii}).ef_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{17}),3);
        events_averg.(fn{ii}).es_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{18}),3);
        events_averg.(fn{ii}).soc_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{19}),3);
        events_averg.(fn{ii}).ee_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{20}),3);
        events_averg.(fn{ii}).er_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{21}),3);
        events_averg.(fn{ii}).drink_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{22}),3);
        events_averg.(fn{ii}).eat_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{23}),3);
        events_averg.(fn{ii}).blockend_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{24}),3);
        events_averg.(fn{ii}).exd_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{25}),3);
        events_averg.(fn{ii}).exf_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{26}),3);
        events_averg.(fn{ii}).exs_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{27}),3);
        events_averg.(fn{ii}).exe_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{28}),3);
        events_averg.(fn{ii}).exr_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{29}),3);
        events_averg.(fn{ii}).run_averg_bsl_m= nanmean(events_averg.(fn{ii}).(fn1{30}),3); 

end    

%% concatenating over all cells

events_averg.all.ed = cat(2,events_averg.g5.ed_averg_bsl_m,events_averg.g2.ed_averg_bsl_m,events_averg.g4.ed_averg_bsl_m);
events_averg.all.ef = cat(2,events_averg.g5.ef_averg_bsl_m,events_averg.g2.ef_averg_bsl_m,events_averg.g4.ef_averg_bsl_m);
events_averg.all.es = cat(2,events_averg.g5.es_averg_bsl_m,events_averg.g2.es_averg_bsl_m,events_averg.g4.es_averg_bsl_m);
events_averg.all.soc = cat(2,events_averg.g5.soc_averg_bsl_m,events_averg.g2.soc_averg_bsl_m,events_averg.g4.soc_averg_bsl_m);
events_averg.all.ee = cat(2,events_averg.g5.ee_averg_bsl_m,events_averg.g2.ee_averg_bsl_m,events_averg.g4.ee_averg_bsl_m);
events_averg.all.er = cat(2,events_averg.g5.er_averg_bsl_m,events_averg.g2.er_averg_bsl_m,events_averg.g4.er_averg_bsl_m);
events_averg.all.drink = cat(2,events_averg.g5.drink_averg_bsl_m,events_averg.g2.drink_averg_bsl_m,events_averg.g4.drink_averg_bsl_m);
events_averg.all.eat = cat(2,events_averg.g5.eat_averg_bsl_m,events_averg.g2.eat_averg_bsl_m,events_averg.g4.eat_averg_bsl_m);
events_averg.all.be = cat(2,events_averg.g5.blockend_averg_bsl_m,events_averg.g2.blockend_averg_bsl_m,events_averg.g4.blockend_averg_bsl_m);
events_averg.all.exd = cat(2,events_averg.g5.exd_averg_bsl_m,events_averg.g2.exd_averg_bsl_m,events_averg.g4.exd_averg_bsl_m);
events_averg.all.exf = cat(2,events_averg.g5.exf_averg_bsl_m,events_averg.g2.exf_averg_bsl_m,events_averg.g4.exf_averg_bsl_m);
events_averg.all.exs = cat(2,events_averg.g5.exs_averg_bsl_m,events_averg.g2.exs_averg_bsl_m,events_averg.g4.exs_averg_bsl_m);
events_averg.all.exe = cat(2,events_averg.g5.exe_averg_bsl_m,events_averg.g2.exe_averg_bsl_m,events_averg.g4.exe_averg_bsl_m);
events_averg.all.exr = cat(2,events_averg.g5.exr_averg_bsl_m,events_averg.g2.exr_averg_bsl_m,events_averg.g4.exr_averg_bsl_m);
events_averg.all.run = cat(2,events_averg.g5.run_averg_bsl_m,events_averg.g2.run_averg_bsl_m,events_averg.g4.run_averg_bsl_m);

%% averages around events

events_averg.all.ed_averg = mean(events_averg.all.ed(15:30,:),1);
events_averg.all.ef_averg = mean(events_averg.all.ef(15:30,:),1);
events_averg.all.es_averg = mean(events_averg.all.es(15:30,:),1);
events_averg.all.soc_averg = mean(events_averg.all.soc(15:30,:),1);
events_averg.all.ee_averg = mean(events_averg.all.ee(15:30,:),1);
events_averg.all.er_averg = mean(events_averg.all.er(15:30,:),1);
events_averg.all.drink_averg = mean(events_averg.all.drink(30:45,:),1);
events_averg.all.drink_antic_averg = mean(events_averg.all.drink(15:30,:),1);
events_averg.all.eat_averg = mean(events_averg.all.eat(30:45,:),1);
events_averg.all.eat_antic_averg = mean(events_averg.all.eat(15:30,:),1);
events_averg.all.be_averg = mean(events_averg.all.be(15:30,:),1);
events_averg.all.home_averg = mean(events_averg.all.be(30:90,:),1);
events_averg.all.exd_averg = mean(events_averg.all.exd(15:30,:),1);
events_averg.all.exf_averg = mean(events_averg.all.exf(15:30,:),1);
events_averg.all.exs_averg = mean(events_averg.all.exs(15:30,:),1);
events_averg.all.exe_averg = mean(events_averg.all.exe(15:30,:),1);
events_averg.all.exr_averg = mean(events_averg.all.exr(15:30,:),1);
events_averg.all.run_averg = mean(events_averg.all.run(15:30,:),1);


%% scatterplots
 f1= figure; 
 figure(f1);hold on,
    subplot(2,2,1);scatter(events_averg.all.eat_averg,events_averg.all.drink_averg,	"g", LineWidth=1);
    xlabel('\fontsize{12}Ca+-response eat');ylabel('\fontsize{12}drink ');hold on;plot([0 2],[0 0],"k");plot([0 0],[0 2],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-2, 0.5],	"m");
    xlim([-1 2]);  %xlim([-1 1])
    scatter(([events_averg.all.eat_averg(1,(eat_cells(:,1)))]),([events_averg.all.drink_averg(1,(eat_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.eat_averg(1,(eat_cells_inh(:,1)))]),([events_averg.all.drink_averg(1,(eat_cells_inh(:,1)))]),'b', LineWidth=1);
    scatter(([events_averg.all.eat_averg(1,(drink_cells(:,1)))]),([events_averg.all.drink_averg(1,(drink_cells(:,1)))]),"*",'r');scatter(([events_averg.all.eat_averg(1,(drink_cells_inh(:,1)))]),([events_averg.all.drink_averg(1,(drink_cells_inh(:,1)))]),"*",'b');
    
    subplot(2,2,2);scatter(events_averg.all.eat_antic_averg,events_averg.all.drink_antic_averg,	"g", LineWidth=1);hold on;plot([0 2],[0 0],"k");plot([0 0],[0 2],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response eat anticipatory');ylabel('\fontsize{12}drink anticipatory ');
    xlim([-1 2]);
    scatter(([events_averg.all.eat_antic_averg(1,(eat_antic_cells(:,1)))]),([events_averg.all.drink_antic_averg(1,(eat_antic_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.eat_antic_averg(1,(eat_antic_cells_inh(:,1)))]),([events_averg.all.drink_antic_averg(1,(eat_antic_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.eat_antic_averg(1,(drink_antic_cells(:,1)))]),([events_averg.all.drink_antic_averg(1,(drink_antic_cells(:,1)))]),"*",'r');scatter(([events_averg.all.eat_antic_averg(1,(drink_antic_cells_inh(:,1)))]),([events_averg.all.drink_antic_averg(1,(drink_antic_cells_inh(:,1)))]),"*",'b');
    
    subplot(2,2,3);scatter(events_averg.all.eat_averg,events_averg.all.soc_averg,"g", LineWidth=1);hold on;plot([0 2],[0 0],"k");plot([0 0],[0 2],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-2, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response eat');ylabel('\fontsize{12}social contact');
    xlim([-1 2]);
    scatter(([events_averg.all.eat_averg(1,(eat_cells(:,1)))]),([events_averg.all.soc_averg(1,(eat_cells(:,1)))]),'r',  LineWidth=1);scatter(([events_averg.all.eat_averg(1,(eat_cells_inh(:,1)))]),([events_averg.all.soc_averg(1,(eat_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.eat_averg(1,(soc_cells(:,1)))]),([events_averg.all.soc_averg(1,(soc_cells(:,1)))]),"*",'r');scatter(([events_averg.all.eat_averg(1,(soc_cells_inh(:,1)))]),([events_averg.all.soc_averg(1,(soc_cells_inh(:,1)))]),"*",'b');
    
    subplot(2,2,4);scatter(events_averg.all.be_averg,events_averg.all.home_averg,"g", LineWidth=1);hold on;plot([0 2],[0 0],"k");plot([0 0],[0 3],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response block end');ylabel('\fontsize{12}homecage');
    xlim([-1 2]);
    scatter(([events_averg.all.be_averg(1,(be_cells(:,1)))]),([events_averg.all.home_averg(1,(be_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.be_averg(1,(be_cells_inh(:,1)))]),([events_averg.all.home_averg(1,(be_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.be_averg(1,(home_cells(:,1)))]),([events_averg.all.home_averg(1,(home_cells(:,1)))]),"*",'r');scatter(([events_averg.all.be_averg(1,(home_cells_inh(:,1)))]),([events_averg.all.home_averg(1,(home_cells_inh(:,1)))]),"*",'b');


 f2= figure; 
 figure(f2);hold on;
    subplot(2,2,1);scatter(events_averg.all.ed_averg,events_averg.all.ef_averg,	"g", LineWidth=1);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");%plot(X,Y)
    xlabel('\fontsize{12}Ca+-response enter drink');ylabel('\fontsize{12}enter food ');
    xlim([-1 1.5]); 
    scatter(([events_averg.all.ed_averg(1,(ed_cells(:,1)))]),([events_averg.all.ef_averg(1,(ed_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.ed_averg(1,(ed_cells_inh(:,1)))]),([events_averg.all.ef_averg(1,(ed_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.ed_averg(1,(ef_cells(:,1)))]),([events_averg.all.ef_averg(1,(ef_cells(:,1)))]),"*",'r');scatter(([events_averg.all.ed_averg(1,(ef_cells_inh(:,1)))]),([events_averg.all.ef_averg(1,(ef_cells_inh(:,1)))]),"*",'b');
    legend({'ed cells ex','ed cells inh','ef cells ex' , 'ef cells inh'});
    subplot(2,2,2);scatter(events_averg.all.ed_averg,events_averg.all.es_averg,	"g", LineWidth=1);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter drink');ylabel('\fontsize{12}enter social ');
    xlim([-1 1.5]);
    scatter(([events_averg.all.ed_averg(1,(ed_cells(:,1)))]),([events_averg.all.es_averg(1,(ed_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.ed_averg(1,(ed_cells_inh(:,1)))]),([events_averg.all.er_averg(1,(ed_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.ed_averg(1,(es_cells(:,1)))]),([events_averg.all.es_averg(1,(es_cells(:,1)))]),"*",'r');scatter(([events_averg.all.ed_averg(1,(es_cells_inh(:,1)))]),([events_averg.all.es_averg(1,(es_cells_inh(:,1)))]),"*",'b');


    subplot(2,2,3);scatter(events_averg.all.ed_averg,events_averg.all.er_averg,	"g", LineWidth=1);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter drink');ylabel('\fontsize{12}enter run');
    xlim([-1 1.5]);
    scatter(([events_averg.all.ed_averg(1,(ed_cells(:,1)))]),([events_averg.all.er_averg(1,(ed_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.ed_averg(1,(ed_cells_inh(:,1)))]),([events_averg.all.er_averg(1,(ed_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.ed_averg(1,(er_cells(:,1)))]),([events_averg.all.er_averg(1,(er_cells(:,1)))]),"*",'r');scatter(([events_averg.all.ed_averg(1,(er_cells_inh(:,1)))]),([events_averg.all.er_averg(1,(er_cells_inh(:,1)))]),"*",'b');


    subplot(2,2,4);scatter(events_averg.all.ed_averg,events_averg.all.be_averg,	"g", LineWidth=1);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter drink');ylabel('\fontsize{12}block end');
    xlim([-1 1.5]);   
    scatter(([events_averg.all.ed_averg(1,(ed_cells(:,1)))]),([events_averg.all.be_averg(1,(ed_cells(:,1)))]),'r', LineWidth=1);scatter(([events_averg.all.ed_averg(1,(ed_cells_inh(:,1)))]),([events_averg.all.ef_averg(1,(ed_cells_inh(:,1)))]),'b',LineWidth=1);
    scatter(([events_averg.all.ed_averg(1,(be_cells(:,1)))]),([events_averg.all.be_averg(1,(be_cells(:,1)))]),"*",'r');scatter(([events_averg.all.be_averg(1,(be_cells_inh(:,1)))]),([events_averg.all.be_averg(1,(be_cells_inh(:,1)))]),"*",'b');



f3= figure; 
 figure(f3);hold on;
    subplot(2,2,1);scatter(events_averg.all.ef_averg,events_averg.all.es_averg);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");%plot(X,Y)
    xlabel('\fontsize{12}Ca+-response enter food');ylabel('\fontsize{12}Ca+-response enter social ');
    xlim([-1 1.5]);  %xlim([-1 1])
    %plot([0 2],[0 0]);plot([0 0],[0 2])
    subplot(2,2,2);scatter(events_averg.all.ef_averg,events_averg.all.ee_averg);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter food');ylabel('\fontsize{12}Ca+-response enter explore ');
    xlim([-1 1.5]);
    subplot(2,2,3);scatter(events_averg.all.ef_averg,events_averg.all.er_averg);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter food');ylabel('\fontsize{12}Ca+-response enter run');
    xlim([-1 1.5]);
    subplot(2,2,4);scatter(events_averg.all.ef_averg,events_averg.all.be_averg);hold on;plot([0 1.5],[0 0],"k");plot([0 0],[0 1.5],"k");plot([-1, 0],[-1,0],"k");%plot([-1, 0.5],[0.5,0.5],	"m");plot([0.5,0.5],[-1, 0.5],	"m");
    xlabel('\fontsize{12}Ca+-response enter food');ylabel('\fontsize{12}Ca+-response block end');
    xlim([-1 1.5]);



%% histogram



%% binary

%% std of the population

std_ed = std(events_averg.all.ed(1:10,:), 1);
std_ef = std(events_averg.all.ef(1:10,:), 1);
std_es = std(events_averg.all.es(1:10,:), 1);
std_ee = std(events_averg.all.ee(1:10,:), 1);
std_er = std(events_averg.all.er(1:10,:), 1);
std_be = std(events_averg.all.be(1:10,:), 1);
%figure; plot(std_ed);hold on;plot(std_ef);plot(std_es);plot(std_es);plot(std_ee);plot(std_er);plot(std_be);

%cll =30;
%figure; plot(events_averg.all.ed(:,cll), LineWidth=2);hold on;plot(events_averg.all.ef(:,cll), LineWidth=2);plot(events_averg.all.es(:,cll), LineWidth=2);plot(events_averg.all.ee(:,cll), LineWidth=2);plot(events_averg.all.er(:,cll), LineWidth=2);plot(events_averg.all.be(1:31,cll), LineWidth=2);

%thresh= 0.5;
%lowthresh = -0.5;
fact= 1.5;
events=[];
column= 0;
fn3 = fieldnames(events_averg.all);
for hh=16:33
    column= column+1;
    for kk=1:59
        if ((hh == 16) || (hh ==22) || (hh ==23) || (hh ==18)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ed)
            events(column,kk) = 1;
            %kk = kk + 1;
        elseif  ((hh == 16) || (hh ==22) ||(hh ==23)||(hh ==18)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ed))
               events(column,kk) = 2;
  
        elseif  ((hh == 17)||(hh == 24)||(hh == 25)||(hh == 29)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ef)
               events(column,kk) = 1;
        elseif  ((hh == 17)||(hh == 24)||(hh == 25)||(hh == 29)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ef))
               events(column,kk) = 2;    

        elseif ( (hh == 18)||(hh == 19)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_es)
               events(column,kk) = 1;
        elseif  ((hh == 18)||(hh == 19)) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_es))
               events(column,kk) = 2; 


        elseif  (hh == 20) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_ee)
               events(column,kk) = 1;
        elseif  (hh == 20) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_ee))
               events(column,kk) = 2; 

        elseif ( (hh == 21)||(hh == 33)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_er)
               events(column,kk) = 1;
        elseif ( (hh == 21)||(hh == 33) ) & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_er))
               events(column,kk) = 2; 

        elseif ( (hh == 26)||(hh == 27)) & (events_averg.all.(fn3{hh})(1,kk) >= fact*std_be)
               events(column,kk) = 1;
        elseif ( (hh == 26)||(hh == 27))  & (events_averg.all.(fn3{hh})(1,kk) <= -(fact*std_be))
               events(column,kk) = 2;        

        else    
            events(column,kk) = 0;
        end    
    end
    
end

[r,c]=find(events==1);%finds ones
[r2,c2]=find(events ==2);%finds 2

ed_cells =c(find(r==1));ef_cells =c(find(r==2));es_cells =c(find(r==3));ee_cells =c(find(r==5));er_cells =c(find(r==6));
drink_antic_cells =c(find(r==8));drink_cells =c(find(r==7));eat_antic_cells =c(find(r==10));eat_cells =c(find(r==9));
be_cells =c(find(r==11));home_cells =c(find(r==12));soc_cells=c(find(r==4));run_cells=c(find(r==18));
exd_cells =c(find(r==13));exf_cells =c(find(r==14));exs_cells =c(find(r==15));exe_cells =c(find(r==16));exr_cells =c(find(r==17));
drink_antic_cells =drink_antic_cells( [3 4 5 7:15], :); eat_antic_cells= eat_antic_cells([1:7 10:13],:); %without entry cells

ed_cells_inh =c2(find(r2==1));ef_cells_inh =c2(find(r2==2));es_cells_inh =c2(find(r2==3));ee_cells_inh =c2(find(r2==5));er_cells_inh =c2(find(r2==6));
drink_antic_cells_inh =c2(find(r2==8));drink_cells_inh =c2(find(r2==7));eat_antic_cells_inh =c2(find(r2==10));eat_cells_inh =c2(find(r2==9));
be_cells_inh =c2(find(r2==11));home_cells_inh =c2(find(r2==12));soc_cells_inh=c2(find(r2==4));run_cells_inh=c2(find(r2==18));
exd_cells_inh =c2(find(r2==13));exf_cells_inh =c2(find(r2==14));exs_cells_inh =c2(find(r2==15));exe_cells_inh =c2(find(r2==16));exr_cells_inh =c2(find(r2==17));
drink_antic_cells_inh =drink_antic_cells_inh( [2 4:7 9:14 16:17], :); eat_antic_cells_inh= eat_antic_cells_inh; %without entry cells

%% percentages
all=59;
%ev = er_cells;
%perc = (size(ev,1)/all)*100);
%size(ev,1)
events_averg.perc.ed_cells = [(size(ed_cells,1)/all)*100 ; size(ed_cells,1)];events_averg.perc.exd_cells = [(size(exd_cells,1)/all)*100 ; size(exd_cells,1)];
events_averg.perc.ef_cells = [(size(ef_cells,1)/all)*100 ; size(ef_cells,1)];events_averg.perc.exf_cells = [(size(exf_cells,1)/all)*100 ; size(exf_cells,1)];
events_averg.perc.es_cells = [(size(es_cells,1)/all)*100 ; size(es_cells,1)];events_averg.perc.exs_cells = [(size(exs_cells,1)/all)*100 ; size(exs_cells,1)];
events_averg.perc.ee_cells = [(size(ee_cells,1)/all)*100 ; size(ee_cells,1)];events_averg.perc.exe_cells = [(size(exe_cells,1)/all)*100 ; size(exe_cells,1)];
events_averg.perc.er_cells = [(size(er_cells,1)/all)*100 ; size(er_cells,1)];events_averg.perc.exr_cells = [(size(exr_cells,1)/all)*100 ; size(exr_cells,1)];
events_averg.perc.be_cells = [(size(be_cells,1)/all)*100 ; size(be_cells,1)];events_averg.perc.home_cells = [(size(home_cells,1)/all)*100 ; size(home_cells,1)];
events_averg.perc.drink_cells = [(size(drink_cells,1)/all)*100 ; size(drink_cells,1)];events_averg.perc.drink_antic_cells = [(size(drink_antic_cells,1)/all)*100 ; size(drink_antic_cells,1)];
events_averg.perc.eat_cells = [(size(eat_cells,1)/all)*100 ; size(eat_cells,1)];events_averg.perc.eat_antic_cells = [(size(eat_antic_cells,1)/all)*100 ; size(eat_antic_cells,1)];
events_averg.perc.soc_cells = [(size(soc_cells,1)/all)*100 ; size(soc_cells,1)];
events_averg.perc.run_cells = [(size(run_cells,1)/all)*100 ; size(run_cells,1)];
%inhib cells
events_averg.perc.ed_cells_inh = [(size(ed_cells_inh,1)/all)*100 ; size(ed_cells_inh,1)];events_averg.perc.exd_cells_inh = [(size(exd_cells_inh,1)/all)*100 ; size(exd_cells_inh,1)];
events_averg.perc.ef_cells_inh = [(size(ef_cells_inh,1)/all)*100 ; size(ef_cells_inh,1)];events_averg.perc.exf_cells_inh = [(size(exf_cells_inh,1)/all)*100 ; size(exf_cells_inh,1)];
events_averg.perc.es_cells_inh = [(size(es_cells_inh,1)/all)*100 ; size(es_cells_inh,1)];events_averg.perc.exs_cells_inh = [(size(exs_cells_inh,1)/all)*100 ; size(exs_cells_inh,1)];
events_averg.perc.ee_cells_inh = [(size(ee_cells_inh,1)/all)*100 ; size(ee_cells_inh,1)];events_averg.perc.exe_cells_inh = [(size(exe_cells_inh,1)/all)*100 ; size(exe_cells_inh,1)];
events_averg.perc.er_cells_inh = [(size(er_cells_inh,1)/all)*100 ; size(er_cells_inh,1)];events_averg.perc.exr_cells_inh = [(size(exr_cells_inh,1)/all)*100 ; size(exr_cells_inh,1)];
events_averg.perc.be_cells_inh = [(size(be_cells_inh,1)/all)*100 ; size(be_cells_inh,1)];events_averg.perc.home_cells_inh = [(size(home_cells_inh,1)/all)*100 ; size(home_cells_inh,1)];
events_averg.perc.drink_cells_inh = [(size(drink_cells_inh,1)/all)*100 ; size(drink_cells_inh,1)];events_averg.perc.drink_antic_cells_inh = [(size(drink_antic_cells_inh,1)/all)*100 ; size(drink_antic_cells_inh,1)];
events_averg.perc.eat_cells_inh = [(size(eat_cells_inh,1)/all)*100 ; size(eat_cells_inh,1)];events_averg.perc.eat_antic_cells_inh = [(size(eat_antic_cells_inh,1)/all)*100 ; size(eat_antic_cells_inh,1)];
events_averg.perc.soc_cells_inh = [(size(soc_cells_inh,1)/all)*100 ; size(soc_cells_inh,1)];
events_averg.perc.run_cells_inh = [(size(run_cells_inh,1)/all)*100 ; size(run_cells_inh,1)];
%% piecharts
noEventCells = find((sum(events,1)) == 0);

figure;label={'2% no event cells', '98% event cells'};pie([size(noEventCells,1)/59,(59-(size(noEventCells,1)))/59], label);


%figure;label={'drink cells', 'non-drink cells'},pie([size(drink_cells,1)/59,(59-(size(drink_cells,1)))/59],label); 
%figure;label={'food cells', 'non-food cells'},pie([size(eat_cells,1)/59,(59-(size(eat_cells,1)))/59],label);

%% venn diagrams
A = 1:59;
B = drink_cells;
C = drink_antic_cells;
D = eat_cells;
E = eat_antic_cells;
F = soc_cells;
G = se_cells;
H = home_cells;
I = de_cells;
setListData = {A B C D E F G H};
setLabels = ["all"; "drink"; "drink antic"; "eat"; "eat antic"; "social"; "social entry"; "home"];
h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);
%h.ShowIntersectionCounts = true;

%entries
setListData2={A be_cells ef_cells ed_cells es_cells er_cells};
setLabels2=[" "; "block end 8.5% (5/59)" ;"enter food 6.8% (4/59) ";"enter drink 6.8% (4/59)";"enter social 25.4% (15/59)";"enter run 3.9% (2/59)"];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);
%hh.ShowIntersectionCounts = true;
setListData2={A es_cells_inh };
setLabels2=[" " ; "es "];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);
%
setListData2={ef_cells_inh ed_cells_inh es_cells_inh ee_cells_inh};
setLabels2=["\fontsize{17}enter food 3.9% (2/59) ";"\fontsize{20}enter drink 8.5% (5/59)";"\fontsize{20}enter social 27.1% (16/59)";"\fontsize{20}enter explore 5.1% (3/59)"];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);
%ex cons
setListData2={A drink_cells drink_antic_cells eat_cells eat_antic_cells soc_cells home_cells};
setLabels2=["" ;"\fontsize{17}drink 17.0% (10/59)"; " \fontsize{17}drink ant 20.2% (12/59) "; "\fontsize{17}food 35.6% (21/59)"; "\fontsize{17}food ant 18.6% (11/59)"; "\fontsize{17}social contact 28.8% (17/59)"; "\fontsize{17}home 17.0% (10/59) "];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);
%inh cons
setListData2={A drink_cells_inh drink_antic_cells_inh eat_cells_inh eat_antic_cells_inh soc_cells_inh home_cells_inh};
setLabels2=[ " "; "\fontsize{17}drink 30.5% (18/59)"; " \fontsize{17}drink ant 22.0% (13/59) "; "\fontsize{17}food 18.6% (11/59)"; "\fontsize{17}food ant 15.3% (9/59)"; "\fontsize{17}social contact 39.0% (23/59)"; "\fontsize{17}home 3.4% (2/59) "];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);

setListData2={A drink_cells drink_antic_cells ed_cells};
setLabels2=[ " "; "\fontsize{15}drink 17.0% (10/59)"; " \fontsize{15}drink ant 20.2% (12/59) "; "\fontsize{15}enter drink 6.8% (4/59)"];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);

setListData2={A drink_cells_inh drink_antic_cells_inh ed_cells_inh};
setLabels2=[ " ";"\fontsize{15}drink inh 30.5% (18/59)";" \fontsize{15}drink ant inh 22.0% (13/59) ";"\fontsize{15}enter drink inh 8.5% (5/59)"];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);


setListData2={A home_cells drink_cells_inh};
setLabels2=[ " ";"\fontsize{17}home ex 17.0% (10/59) ";"\fontsize{17}food ex 35.6% (21/59)";"\fontsize{15}drink inh 30.5% (18/59)"];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);

setListData2={A soc_cells soc_cells_inh es_cells es_cells_inh};
setLabels2=[ " ";"\fontsize{15}social contact ex 28.8% (17/59)" ; "\fontsize{15}social contact inh 39.0% (23/59)"; "\fontsize{15}enter social ex 25.4% (15/59)";"\fontsize{15}enter social inh 27.1% (16/59)"];
hh=vennEulerDiagram(setListData2, setLabels2, 'drawProportional', true);

%% traces
enter = 1606;
c40_ed = ci_data.g4.Sep02(4).session.raw_f_res(enter-200:enter+100,4);z_c40 = (c40_ed-nanmean(c40_ed(1:50,1)))/nanstd(c40_ed(1:50,1));
c52_ed = ci_data.g4.Sep02(4).session.raw_f_res(enter-200:enter+100,14);z_c52 = (c52_ed-nanmean(c52_ed(1:50,1)))/nanstd(c52_ed(1:50,1));
c56_ed = ci_data.g4.Sep02(4).session.raw_f_res(enter-200:enter+100,16);z_c56 = (c56_ed-nanmean(c56_ed(1:50,1)))/nanstd(c56_ed(1:50,1));
c54_ed = ci_data.g4.Sep02(4).session.raw_f_res(enter-200:enter+100,20);z_c54 = (c54_ed-nanmean(c54_ed(1:50,1)))/nanstd(c54_ed(1:50,1));
%x-nanmean(b)/nanstd(b)
exmpl_ed = table(smooth(zscore(c40_ed)), smooth(zscore(c52_ed)), smooth(zscore(c56_ed)), smooth(zscore(c54_ed)));
%figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
%plot(smooth(zscore(c40_ed)));hold on;plot(smooth(zscore(c52_ed)));plot(smooth(zscore(c56_ed)));plot(smooth(zscore(c54_ed)));
p1=plot(smooth(z_c40));hold on;p2=plot(smooth(z_c52));p3=plot(smooth(z_c56));p4=plot(smooth(z_c54));
p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;p4.LineWidth=1;
xl= xline(199,'-','entry drink');xl.LineWidth = 1; 
xxl= xline(199+20,'-','drink');xxl.LineWidth = 1;
xxxl= xline(200-185,'-','eat');xxxl.LineWidth = 1;
xxxl= xline(200-192,'-','');xxxl.LineWidth = 1;
xxxxl= xline(200+90,'-','exit drink');xxxxl.LineWidth = 1;
xticks([0 50 100 150 200 250 300]);xticklabels({'-20','-15','-10','-5','0','5','10'});xlabel('\fontsize{12}time [s]')
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');xlim([0 300]);
legend({'\fontsize{12}c40 (ed ex)','\fontsize{12}c50', '\fontsize{12}c52', '\fontsize{12}c56 (ef ex'}, Location="northwest");



cell_5_ed = events_averg.g4.ed_averg_bsl(:,5,:);
c5 = reshape(cell_5_ed ,[],70);

figure; plot(c5, LineWidth=1);hold on; plot(events_averg.all.ed(:,5),'k' ,LineWidth=4);
xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time [s]');xlim([1 31]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');title(['\fontsize{10}entry drink inhibited cell']);



cell_5_ed = events_averg.g5.soc_averg_bsl(:,14,:);
c5 = reshape(cell_5_ed ,[],145);c5 = c5(:,25:75);

figure; plot(c5, LineWidth=1);hold on; plot(events_averg.all.ed(:,5),'k' ,LineWidth=4);
xticks([1 11 21 31]);xticklabels({'-3','-2','-1', '0'});xlabel('\fontsize{12}time [s]');xlim([1 31]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');title(['\fontsize{10}entry drink inhibited cell']);



%% enter food g5
enter = 2337;
c40_ed = ci_data.g5.Aug08(17).session.raw_f_res(enter-200:enter+100,4);z_c40 = (c40_ed-nanmean(c40_ed(1:50,1)))/nanstd(c40_ed(1:50,1));
c52_ed = ci_data.g5.Aug08(17).session.raw_f_res(enter-200:enter+100,8);z_c52 = (c52_ed-nanmean(c52_ed(1:50,1)))/nanstd(c52_ed(1:50,1));
c56_ed = ci_data.g5.Aug08(17).session.raw_f_res(enter-200:enter+100,6);z_c56 = (c56_ed-nanmean(c56_ed(1:50,1)))/nanstd(c56_ed(1:50,1));
c54_ed = ci_data.g5.Aug08(17).session.raw_f_res(enter-200:enter+100,2);z_c54 = (c54_ed-nanmean(c54_ed(1:50,1)))/nanstd(c54_ed(1:50,1));
%x-nanmean(b)/nanstd(b)
exmpl_ed = table(smooth(zscore(c40_ed)), smooth(zscore(c52_ed)), smooth(zscore(c56_ed)), smooth(zscore(c54_ed)));
%figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
%plot(smooth(zscore(c40_ed)));hold on;plot(smooth(zscore(c52_ed)));plot(smooth(zscore(c56_ed)));plot(smooth(zscore(c54_ed)));
p1=plot(smooth(z_c40));hold on;p2=plot(smooth(z_c52));p3=plot(smooth(z_c56));%p4=plot(smooth(z_c54));
p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;%p4.LineWidth=1;
xl= xline(200,'-','entry food');xl.LineWidth = 1; 
xxl = xline(200+24,'-','cons food');xxl.LineWidth = 1;
xxxl = xline(200-34,'-','exit social');xxxl.LineWidth = 1;
xticks([0 50 100 150 200 250 300]);xticklabels({'-20','-15','-10','-5','0','5','10'});xlabel('\fontsize{12}time [s]');xlim([0 300]);
ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');xlim([0 300]);
legend({'\fontsize{12}c4 (eat ex)','\fontsize{12}c8 (eat ex)', '\fontsize{12}c6 (ef/eat inh)'}, Location="northwest");

%%
enter = 9012;
c40_ed = ci_data.g5.Aug27(10).session.raw_f_res(enter-3000:enter,18);%z_c40 = (c40_ed-nanmean(c40_ed(1:50,1)))/nanstd(c40_ed(1:50,1));
c52_ed = ci_data.g5.Aug27(10).session.raw_f_res(enter-3000:enter,9);%z_c52 = (c52_ed-nanmean(c52_ed(1:50,1)))/nanstd(c52_ed(1:50,1));
c56_ed = ci_data.g5.Aug27(10).session.raw_f_res(enter-3000:enter,10);%z_c56 = (c56_ed-nanmean(c56_ed(1:50,1)))/nanstd(c56_ed(1:50,1));
c54_ed = ci_data.g5.Aug27(10).session.raw_f_res(enter-3000:enter,12);%z_c54 = (c54_ed-nanmean(c54_ed(1:50,1)))/nanstd(c54_ed(1:50,1));
%x-nanmean(b)/nanstd(b)
exmpl_ed = table(smooth(zscore(c40_ed)), smooth(zscore(c52_ed)), smooth(zscore(c56_ed)), smooth(zscore(c54_ed)));
figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
plot(smooth(zscore(c40_ed)+1));hold on;plot(smooth(zscore(c52_ed)+5));plot(smooth(zscore(c56_ed)+9));plot(smooth(zscore(c54_ed)+13));
p1 = plot([550 250],[-2.5 -2.5],"k");p2=plot([250 250],[-2.5 -0.5],"k");p1.LineWidth = 2; p2.LineWidth = 2; 
text(320,-3,'0.5s');text(30,-1,'2 s.d.');
xlim([0 3000]);
%p1=plot(smooth(z_c40));hold on;p2=plot(smooth(z_c52));p3=plot(smooth(z_c56));p4=plot(smooth(z_c54));
%p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;p4.LineWidth=1;
%xl= xline(200,'-','entry drink');xl.LineWidth = 1; 
xticks([0 600 1200 1800 2400 3000]);xticklabels({'0','1','2','3','4','5'});xlabel('\fontsize{12}time [s]')
%ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');xlim([0 3000]);


enter = 9012;
c40_ed = ci_data.g5.Aug27(10).session.raw_f_res(enter-3000:enter,16);%z_c40 = (c40_ed-nanmean(c40_ed(1:50,1)))/nanstd(c40_ed(1:50,1));
c52_ed = ci_data.g5.Aug27(10).session.raw_f_res(enter-3000:enter,4);%z_c52 = (c52_ed-nanmean(c52_ed(1:50,1)))/nanstd(c52_ed(1:50,1));
%c56_ed = ci_data.g5.Aug27(10).session.raw_f_res(enter-3000:enter,10);%z_c56 = (c56_ed-nanmean(c56_ed(1:50,1)))/nanstd(c56_ed(1:50,1));
%c54_ed = ci_data.g5.Aug27(10).session.raw_f_res(enter-3000:enter,12);%z_c54 = (c54_ed-nanmean(c54_ed(1:50,1)))/nanstd(c54_ed(1:50,1));
%x-nanmean(b)/nanstd(b)
exmpl_ed = table(smooth(zscore(c40_ed)), smooth(zscore(c52_ed)), smooth(zscore(c56_ed)), smooth(zscore(c54_ed)));
figure; stackedplot(exmpl_ed);
figure;set(0,'defaultAxesFontSize',13);
plot(smooth(zscore(c40_ed)+1));hold on;plot(smooth(zscore(c52_ed)+5));%plot(smooth(zscore(c56_ed)+9));plot(smooth(zscore(c54_ed)+13));
p1 = plot([550 250],[-2.5 -2.5],"k");p2=plot([250 250],[-2.5 -0.5],"k");p1.LineWidth = 2; p2.LineWidth = 2; 
text(320,-3,'0.5s');text(30,-1,'2 s.d.');
xlim([0 3000]);
%p1=plot(smooth(z_c40));hold on;p2=plot(smooth(z_c52));p3=plot(smooth(z_c56));p4=plot(smooth(z_c54));
%p1.LineWidth=1;p2.LineWidth=1;p3.LineWidth=1;p4.LineWidth=1;
%xl= xline(200,'-','entry drink');xl.LineWidth = 1; 
xticks([0 600 1200 1800 2400 3000]);xticklabels({'0','1','2','3','4','5'});xlabel('\fontsize{12}time [s]')
%ylabel('\fontsize{12}Ca-activity (z-scored, s.d.)');xlim([0 3000]);
%% baseline and thresh fig


figure; plot(events_averg.all.ed (:, [14 16 5 10 57]), LineWidth=2);hold on; xlim([1 30]); 
xticks([1 5 10 15 20 25 30]); xticklabels({'-3','-2.5' ,'-2','-1.5' ,'-1' , '-0.5', '0'});xlabel('\fontsize{12}time from drink entry [s]'); ylabel('\fontsize{12}averg Ca-activity (z-scored, s.d.)');
plot([1 10], [-0.9 -0.9],"k", LineWidth=3);plot([15 30], [-0.9 -0.9],"k", LineWidth=3);
text(2.5,-0.8,'baseline');text(16,-0.8,'event window');
legend({'\fontsize{12}c14(enter drink ex)','\fontsize{12}c16(enter drink ex)', '\fontsize{12}c5(enter drink inh)', '\fontsize{12}c10 (enter drink inh)','\fontsize{12}c57 (enter food ex)'}, Location="northwest");


%%
drink_ant = events(8,:);
food_ant = events(10,:);
x=table(drink_ant',food_ant');

x= table([2;6],[8;41]);
[hv,p,stats] = fishertest(x)


%%
cell_31_drink = events_averg.g2.drink_averg(:,9,:);
c31 = reshape(cell_31_drink,[],46);

%lim = [-1 9];
fA = figure;
figure(fA);hold on;
        s1= subplot(2,1,1);imagesc([0:3],[1:size(c31,2)],c31');xlim([0 3]);
        title(['\fontsize{14}Drink cell']);xticklabels({'','','','','','',''});ylabel('\fontsize{12}Trials');
        c = colorbar(s1,'Position',[0.93 0.583837209302326 0.022 0.341162790697675]);
        s2= subplot(2,1,2);stdshade(c31',0.3,[0 0 1]);
        xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]')
        ylabel('\fontsize{12}Ca-activity (z-score)');
        set(s2, 'Position', [0.13,0.38,0.775,0.20]);


cell_42_eat = events_averg.g4.eat_averg(:,7,:);
c42 = reshape(cell_42_eat,[],110);
c42 = c42(:,25:100);

%lim = [-1 9];
fA = figure;
figure(fA);hold on;
        s1= subplot(2,1,1);imagesc([0:3],[1:size(c42,2)],c42');xlim([0 3]);
        title(['\fontsize{14}Food anticipatory cell']);xticklabels({'','','','','','',''});ylabel('\fontsize{12}Trials');
        c = colorbar(s1,'Position',[0.93 0.583837209302326 0.022 0.341162790697675]);
        s2= subplot(2,1,2);stdshade(c42',0.3,[0 0 1]);
        xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]')
        ylabel('\fontsize{12}Ca-activity (z-score)');
        set(s2, 'Position', [0.13,0.38,0.775,0.20]); 






cell_43_eat = events_averg.g4.eat_averg(:,7,:);
c43 = reshape(cell_43_eat,[],119);
c43 = c43(:,25:100);

%lim = [-1 9];
fA = figure;
figure(fA);hold on;
        s1= subplot(2,1,1);imagesc([0:3],[1:size(c43,2)],c43');xlim([0 3]);
        title(['\fontsize{14}enter food']);xticklabels({'','','','','','',''});ylabel('\fontsize{12}Trials');
        c = colorbar(s1,'Position',[0.93 0.583837209302326 0.022 0.341162790697675]);
        s2= subplot(2,1,2);stdshade(c43',0.3,[0 0 1]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xlabel('\fontsize{12}time [s]')
        ylabel('\fontsize{12}Ca-activity (z-score)');
        set(s2, 'Position', [0.13,0.38,0.775,0.20]);


cell_40_ed = events_averg.g4.ed_averg(:,4,:);
c40 = reshape(cell_40_ed,[],70);
c40 = c40(:,1:70);

%lim = [-1 9];
fA = figure;
figure(fA);hold on;
        s1= subplot(2,1,1);imagesc([0:3],[1:size(c40,2)],c40');xlim([0 3]);
        title(['\fontsize{14} entry drink']);xticklabels({'','','','','','',''});ylabel('\fontsize{12}Trials');
        c = colorbar(s1,'Position',[0.93 0.583837209302326 0.022 0.341162790697675]);
        s2= subplot(2,1,2);stdshade(c40',0.3,[0 0 1]);
        xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]');
        ylabel('\fontsize{12}Ca-activity (z-score)');
        set(s2, 'Position', [0.13,0.38,0.775,0.20]);


%
cell_14_ed = events_averg.g5.ed_averg(:,14,:);
c14 = reshape(cell_14_ed,[],82);
c14 = c14(:,1:70);

%lim = [-1 9];
fA = figure;
figure(fA);hold on;
        s1= subplot(2,1,1);imagesc([0:3],[1:size(c14,2)],c14');xlim([0 3]);
        title(['\fontsize{14} entry drink']);xticklabels({'','','','','','',''});ylabel('\fontsize{12}Trials');
        c = colorbar(s1,'Position',[0.93 0.583837209302326 0.022 0.341162790697675]);
        s2= subplot(2,1,2);stdshade(c40',0.3,[0 0 1]);
        xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]');
        ylabel('\fontsize{12}Ca-activity (z-score)');
        set(s2, 'Position', [0.13,0.38,0.775,0.20]);



       

%[xLeft, yBottom, width, height]

%stdshade(amatrix,alpha,acolor,F,smth)


%%
cell_31_drink = events_averg.g2.drink_averg(:,9,:);
c31 = reshape(cell_31_drink,[],46);

%lim = [-1 9];
fA = figure;
figure(fA);hold on;
        s1= subplot(2,1,1);imagesc([0:3],[1:size(c31,2)],c31');xlim([0 3]);
        title(['\fontsize{14}Drink cell']);xticklabels({'','','','','','',''});ylabel('\fontsize{12}Trials');
        c = colorbar(s1,'Position',[0.93 0.583837209302326 0.022 0.341162790697675]);
        s2= subplot(2,1,2);stdshade(c31',0.3,[0 0 1]);
        xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]')
        ylabel('\fontsize{12}Ca-activity (z-score)');
        set(s2, 'Position', [0.13,0.38,0.775,0.20]);


win= 30;
lim = [-2 2];
limaverg = [-0.3 0.5];
f7=figure;% averages across 2all sessions
figure(f7);hold on;
        subplot(2,6,1),imagesc([-win:0],[1:size(events_averg.all.ed,2)],events_averg.all.ed');xline(-20,'--b');
        title(['\fontsize{10}enter drink',' trials: 70-88']);clim(lim)
        xticklabels({' ',' ',' ',''});ylabel('\fontsize{15}Cells');
        sp1= subplot(2,6,7), stdshade(events_averg.all.ed',0.3,[0 0 1]);set(sp1, 'Position', [0.13,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-scored)');
        ylim(limaverg);
        subplot(2,6,2),imagesc([-win:0],[1:size(events_averg.all.ef,2)],events_averg.all.ef');xline(-20,'--b');
        title(['\fontsize{10}enter feed',' trials: 107-160']);clim(lim)
        xticklabels({' ',' ',' ',''});
        sp2= subplot(2,6,8), stdshade(events_averg.all.ef',0.3,[0 0 1]);set(sp2, 'Position', [0.2645,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
        subplot(2,6,3),imagesc([-win:0],[1:size(events_averg.all.es,2)],events_averg.all.es');xline(-20,'--b');
        title(['\fontsize{10}enter social',' trials: 85-149']);clim(lim)
         xticklabels({' ',' ',' ',''});xlabel('\fontsize{15}time [s]');
        sp3= subplot(2,6,9), stdshade(events_averg.all.es',0.3,[0 0 1]);set(sp3, 'Position', [0.3991,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
        subplot(2,6,4),imagesc([-win:0],[1:size(events_averg.all.ee,2)],events_averg.all.ee');xline(-20,'--b');
        title(['\fontsize{10}enter explore',' trials: 20-33']);clim(lim)
         xticklabels({' ',' ',' ',''});
        sp4= subplot(2,6,10), stdshade(events_averg.all.ee',0.3,[0 0 1]);set(sp4, 'Position', [0.5336,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
        subplot(2,6,5),imagesc([-win:0],[1:size(events_averg.all.er,2)],events_averg.all.er');xline(-20,'--b');
        title(['\fontsize{10}enter run',' trials: 8-34']);clim(lim)
        xticklabels({'-3','-2','-1',' '});
        sp5= subplot(2,6,11), stdshade(events_averg.all.er',0.3,[0 0 1]);set(sp5, 'Position', [0.6682,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
        subplot(2,6,6),imagesc([-win:0],[1:size(events_averg.all.be(1:31,:),2)],events_averg.all.be(1:31,:)');xline(-20,'--b');
        title(['\fontsize{10}block end',' trials: 11-72']);clim(lim)
        xticks([-30 -20 -10 0]); xticklabels({' ',' ',' ',' '});
        sp6= subplot(2,6,12), stdshade((events_averg.all.be(1:31,:))',0.3,[0 0 1]);set(sp6, 'Position', [0.8027,0.43,0.099,0.15]);
        xticks([1 11 21 31]);xticklabels({'-3','-2','-1','0'});
        ylim(limaverg);
        h = axes(f7,'visible','off');
        c = colorbar(h,'Position',[0.93 0.425 0.022 0.5]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 3.2;
        caxis(h,lim);

f8=figure;% averages across all sessions
figure(f8);hold on;
        subplot(1,5,1),imagesc([-win:win],[1:size(events_averg.all.exd,2)],events_averg.all.exd');xline(0,'--w');
        title(['\fontsize{10}exit drink',' trials: 70-88']);clim(lim)
        xticks([-30 -20 -10 0 10 20 30]);xticklabels({'-3','-2','-1','0','1','2','3'});ylabel('\fontsize{15}Cells');
        subplot(1,5,2),imagesc([-win:win],[1:size(events_averg.all.exf,2)],events_averg.all.exf');xline(0,'--w');
        title(['\fontsize{10}exit feed',' trials: 109-160']);clim(lim)
        xticks([-30 -20 -10 0 10 20 30]);xticklabels({'-3','-2','-1','0','1','2','3'})
        subplot(1,5,3),imagesc([-win:win],[1:size(events_averg.all.exs,2)],events_averg.all.exs');xline(0,'--w');
        title(['\fontsize{10}exit social',' trials: 85-150';clim(lim)
        xticks([-30 -20 -10 0 10 20 30]);xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{15}time [s]');
        subplot(1,5,4),imagesc([-win:win],[1:size(events_averg.all.exe,2)],events_averg.all.exe');xline(0,'--w');
        title(['\fontsize{10}exit explore',' trials: 25-33']);clim(lim)
        xticks([-30 -20 -10 0 10 20 30]);xticklabels({'-3','-2','-1','0','1','2','3'})
        subplot(1,5,5),imagesc([-win:win],[1:size(events_averg.all.exr,2)],events_averg.all.exr');xline(0,'--w');
        title(['\fontsize{10}exit run',' trials: 9-34']);clim(lim)
        xticks([-30 -20 -10 0 10 20 30]);xticklabels({'-3','-2','-1','0','1','2','3'})
        h = axes(f8,'visible','off');
        c = colorbar(h,'Position',[0.93 0.168 0.022 0.7]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',16,'Rotation',270);
        c.Label.Position(1) = 2.5;
        caxis(h,lim);        

lim = [-3 3];limaverg2= [-1 1];
f9=figure;% averages across all sessions
figure(f9);hold on;
        subplot(2,4,1),imagesc([-win:win],[1:size(events_averg.all.drink,2)],events_averg.all.drink');xline(0,'--w');
        title(['\fontsize{10}drink',' trials: 46-82']);clim(lim)
        xticks([1 11 21 31 41 51 61]); xticklabels({'-3','-2','-1','0','1','2','3'});;ylabel('\fontsize{15}Cells');xlabel('\fontsize{15}time [s]');
        sp7= subplot(2,4,5), stdshade(events_averg.all.drink',0.3,[0 0 1]);set(sp7, 'Position', [0.13,0.43,0.1566,0.15]);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});xlabel('\fontsize{12}time [s]'); ylabel('\fontsize{12}Ca-activity (z-score)');
        ylim(limaverg2);

        subplot(2,4,2),imagesc([-win:win],[1:size(events_averg.all.eat,2)],events_averg.all.eat');xline(0,'--w');
        title(['\fontsize{10}retrieve pellet',' trials: 35-110']);clim(lim)
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});
        sp8= subplot(2,4,6), stdshade(events_averg.all.eat',0.3,[0 0 1]);set(sp8, 'Position', [0.3361,0.43,0.1566,0.15]);
        xticks([1 11 21 31 41 51 61]);xticklabels({'-3','-2','-1','0','1','2','3'});
        ylim(limaverg2);

        subplot(2,4,3),imagesc([0:win],[1:size(events_averg.all.soc,2)],events_averg.all.soc');
        title(['\fontsize{10}social interaction',' trials: 85-149']);clim(lim)
        xticks([1 11 21 31]); xticklabels({' ','1','2',' '});
        sp9= subplot(2,4,7), stdshade(events_averg.all.soc',0.3,[0 0 1]);set(sp9, 'Position', [0.5422,0.43,0.1566,0.15]);
        xticks([1 11 21 31]);xticklabels({'0','-1','-2','-3'}); 
        ylim(limaverg2);
%         subplot(1,4,4),imagesc([-win:win],[1:size(events_averg.all.run,2)],events_averg.all.run');
%         title(['\fontsize{10}run',' trials:',num2str(size(events_averg.g5.run_averg_bsl,3)+size(events_averg.g2.run_averg_bsl,3)+ size(events_averg.g4.run_averg_bsl,3))]);clim(lim)

        subplot(2,4,4),imagesc([0:100],[1:size(events_averg.all.be(31:131,:),2)],events_averg.all.be(31:131,:)');xline(0,'--w');
        title(['\fontsize{10}home',' trials: 11-72']);clim(lim)
        xticks([1 21 41 61 81 101]);xticklabels({' ','2','4','6','8',' '});
        sp10= subplot(2,4,8), stdshade((events_averg.all.be(31:131,:))',0.3,[0 0 1]);set(sp10, 'Position', [0.7484,0.43,0.1566,0.15]);
        xticks([1 21 41 61 81 101]);xticklabels({'0','2','4','6','8','10'});
        ylim(limaverg2);
%         subplot(1,6,6),imagesc([-win:0],[1:size(imstop_averg_bsl_m,2)],imstop_averg_bsl_m');
%         title(['\fontsize{10}imaging end',' trials:',num2str(size(events_averg.g5.ed_averg_bsl,3)+size(events_averg.g2.ed_averg_bsl,3)+ size(events_averg.g4.ed_averg_bsl,3))]);clim(lim)
        h = axes(f9,'visible','off');
        c = colorbar(h,'Position',[0.93 0.45 0.02 0.45]);
        ylabel(c,'Ca+-response (z-scored, s.d.)','FontSize',10,'Rotation',270);
        c.Label.Position(1) = 2.5;
        caxis(h,lim);        


%% nr of cells histogram
nr_before_ex = 
nr_before_inh = 
nr_cons_ex = 
nr_cons_inh = 

bar(nr_before_ex)
y = [nr_before_ex nr_cons_ex ; nr_before_inh nr_cons_inh] 
b = bar(y,"stacked");

