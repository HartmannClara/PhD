clearvars -except ci_data; close all;
set(0,'defaultAxesFontSize',18);
%%var accumulators
% ed_averg=[];ef_averg=[];es_averg=[];ee_averg=[];er_averg=[];drink_averg=[];eat_averg=[];blockend_averg=[];soc_averg=[];
% exd_averg=[];exf_averg=[];exs_averg=[];exe_averg=[];exr_averg=[];imstop_averg=[];run_averg=[];soc_averg_m=[];
% ed_averg_bsl=[];ef_averg_bsl=[];es_averg_bsl=[];ee_averg_bsl=[];er_averg_bsl=[];drink_averg_bsl=[];eat_averg_bsl=[];blockend_averg_bsl=[];soc_averg_bsl=[];
% exd_averg_bsl=[];exf_averg_bsl=[];exs_averg_bsl=[];exe_averg_bsl=[];exr_averg_bsl=[];imstop_averg_bsl=[];run_averg_bsl=[];soc_averg_bsl_m=[];
%%
%cells
%g5=23;
%g2=13;
%g4=23;



fn=fieldnames(ci_data);
%loop through the fields
for ii=1: numel(fn)
    fn1=fieldnames(ci_data.(fn{ii}));
    ed_averg=[];ef_averg=[];es_averg=[];ee_averg=[];er_averg=[];drink_averg=[];eat_averg=[];blockend_averg=[];soc_averg=[];
    exd_averg=[];exf_averg=[];exs_averg=[];exe_averg=[];exr_averg=[];imstop_averg=[];run_averg=[];soc_averg_m=[];
    ed_averg_bsl=[];ef_averg_bsl=[];es_averg_bsl=[];ee_averg_bsl=[];er_averg_bsl=[];drink_averg_bsl=[];eat_averg_bsl=[];blockend_averg_bsl=[];soc_averg_bsl=[];
    exd_averg_bsl=[];exf_averg_bsl=[];exs_averg_bsl=[];exe_averg_bsl=[];exr_averg_bsl=[];imstop_averg_bsl=[];run_averg_bsl=[];soc_averg_bsl_m=[];
    for j=1: numel(fn1)
        fn2=fieldnames(ci_data.(fn{ii}).(fn1{j}));
        for k=1: numel(fn2)
            for s=1:size(ci_data.(fn{ii}).(fn1{j}),2)
                %access the data
                fn3= fieldnames(ci_data.(fn{ii}).(fn1{j})(s).(fn2{k}));
                event_type = ci_data.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2}).(4);
                new_event_frames= ci_data.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{2}).(5);
                raw_f_res=ci_data.(fn{ii}).(fn1{j})(s).(fn2{k}).(fn3{1});
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

                %%
                win =30;
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
                                    snips_bslined2(:,:,i)=snips2(:,:,i)-entry_s_bsl;
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
                                    snips(:,:,i)=raw_f(trig-win:trig+100,:);%blockend and 10 seconds
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
    %                             elseif trig>win+1 & p==14 & trig<= size(raw_f,1)%imaging stop 
    %                                 snips(:,:,i)=raw_f(trig-win:trig,:);
    %                                 snips_bslined(:,:,i)=snips(:,:,i)-be_bsl;
    %                                 imstop_averg=cat(3,imstop_averg,snips(:,:,i));
    %                                 imstop_averg_bsl=cat(3,imstop_averg_bsl,snips_bslined(:,:,i));
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
                            exf_averg_m= nanmean(exf_averg,3);
                            exs_averg_m= nanmean(exs_averg,3);
                            exe_averg_m= nanmean(exe_averg,3);
                            exr_averg_m= nanmean(exr_averg,3);
                            %imstop_averg_m= nanmean(imstop_averg,3);
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
                            exf_averg_bsl_m= nanmean(exf_averg_bsl,3);
                            exs_averg_bsl_m= nanmean(exs_averg_bsl,3);
                            exe_averg_bsl_m= nanmean(exe_averg_bsl,3);
                            exr_averg_bsl_m= nanmean(exr_averg_bsl,3);
                            %imstop_averg_bsl_m= nanmean(imstop_averg_bsl,3);
                            run_averg_bsl_m= nanmean(run_averg_bsl,3);

                            events_averg.(fn{ii}).ed_averg = ed_averg;
                            events_averg.(fn{ii}).ef_averg = ef_averg;
                            events_averg.(fn{ii}).es_averg = es_averg;
                            events_averg.(fn{ii}).soc_averg = soc_averg;
                            events_averg.(fn{ii}).ee_averg = ee_averg;
                            events_averg.(fn{ii}).er_averg = er_averg;
                            events_averg.(fn{ii}).drink_averg = drink_averg;
                            events_averg.(fn{ii}).eat_averg = eat_averg;
                            events_averg.(fn{ii}).blockend_averg = blockend_averg;
                            events_averg.(fn{ii}).exd_averg = exd_averg;
                            events_averg.(fn{ii}).exf_averg = exf_averg;
                            events_averg.(fn{ii}).exs_averg = exs_averg;
                            events_averg.(fn{ii}).exe_averg = exe_averg;
                            events_averg.(fn{ii}).exr_averg = exr_averg;
                            events_averg.(fn{ii}).run_averg = run_averg;

                            events_averg.(fn{ii}).ed_averg_bsl = ed_averg_bsl;
                            events_averg.(fn{ii}).ef_averg_bsl = ef_averg_bsl;
                            events_averg.(fn{ii}).es_averg_bsl = es_averg_bsl;
                            events_averg.(fn{ii}).soc_averg_bsl = soc_averg_bsl;
                            events_averg.(fn{ii}).ee_averg_bsl = ee_averg_bsl;
                            events_averg.(fn{ii}).er_averg_bsl = er_averg_bsl;
                            events_averg.(fn{ii}).drink_averg_bsl = drink_averg_bsl;
                            events_averg.(fn{ii}).eat_averg_bsl = eat_averg_bsl;
                            events_averg.(fn{ii}).blockend_averg_bsl = blockend_averg_bsl;
                            events_averg.(fn{ii}).exd_averg_bsl = exd_averg_bsl;
                            events_averg.(fn{ii}).exf_averg_bsl = exf_averg_bsl;
                            events_averg.(fn{ii}).exs_averg_bsl = exs_averg_bsl;
                            events_averg.(fn{ii}).exe_averg_bsl = exe_averg_bsl;
                            events_averg.(fn{ii}).exr_averg_bsl = exr_averg_bsl;
                            events_averg.(fn{ii}).run_averg_bsl = run_averg_bsl;
                           



                        end
                    end  
                end


            end    
        end 
    end
end

%%


