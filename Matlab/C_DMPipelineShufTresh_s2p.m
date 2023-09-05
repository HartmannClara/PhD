clearvars -except events_averg ci_data; close all;
set(0,'defaultAxesFontSize',10);
%% Takes in events_averg file of shuffled data and nr of iterations
%% creates 95 and 5 percentile tables for every cell for every event
%09.06.2023
%% for shuffled data
fn=fieldnames(events_averg); Iter = 100;%%% Change nur of iterations here!!!
for ii=1:size(fn,1)
    fn1=fieldnames(events_averg.(fn{ii})); 
        for j=20:numel(fn1) % only process baselined data
            % divide trials by then nr of iterations 
            start = 1;
            Eend = (size(events_averg.(fn{ii}).(fn1{j}),3)/Iter);
            for k =1:Iter
            events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr = events_averg.(fn{ii}).(fn1{j})(:,:,start:Eend);
            start = start + (size(events_averg.(fn{ii}).(fn1{j}),3)/Iter);
            Eend = Eend + (size(events_averg.(fn{ii}).(fn1{j}),3)/Iter);
            end
        end
end

% duplicate events which are analyzed in multiple ways
fn = fieldnames(events_averg.shuf);
for ii = 1:size(fn,1)
    events_averg.shuf.(fn{ii}).be_averg = events_averg.shuf.(fn{ii}).blockend_averg_bsl;
    events_averg.shuf.(fn{ii}).home_averg =events_averg.shuf.(fn{ii}).blockend_averg_bsl;
    events_averg.shuf.(fn{ii}).drink_antic = events_averg.shuf.(fn{ii}).drink_averg_bsl;
    events_averg.shuf.(fn{ii}).eat_antic = events_averg.shuf.(fn{ii}).eat_averg_bsl;
    events_averg.shuf.(fn{ii}) = rmfield(events_averg.shuf.(fn{ii}), 'blockend_averg_bsl');% remove blocked field
end    
events_averg1= events_averg;
%make averages
%% makes an average over all trials per cell
fn=fieldnames(events_averg.shuf);
%loop through the fields to make mean all bsl trials
for ii=1: numel(fn)
    fn1=fieldnames(events_averg.shuf.(fn{ii}));
        for j=1:numel(fn1)
            %different analysis windows for 
            arr = events_averg.shuf.(fn{ii}).(fn1{j}).Nr;empt = isempty(arr); % if array is empty move on
            if empt == true
                j=j+1;
            elseif j == 4 || j ==15%soc expl
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(10:25,:),1);% takes 1.5 s window 1s after beam break,most likelihood 
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            elseif j == 7 || j == 8 %drink eat
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(31:45,:),1);% 1.5s from beam break 
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            elseif   j ==20%home
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(30:89,:),1);% 6s window 3 s after exit beam
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            elseif  j == 21 || j == 22% antic
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(16:30,:),1);%1.5s right before beam break
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            else % all the other ones 
                for k = 1:Iter
                avergTrials = nanmean(events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr,3);
                avergTrialsMean = nanmean(avergTrials(16:30,:),1);
                events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean = avergTrialsMean;
                end
            end 
        end
end 
%manually delete run

fn=fieldnames(events_averg.shuf);
%combine shuffled means in one array
for ii=1:numel(fn1)
    fn1=fieldnames(events_averg.shuf.(fn{ii}));
        for j=1:numel(fn1)
            arr = events_averg.shuf.(fn{ii}).(fn1{j})(k).Nr;empt = isempty(arr); % if array is empty move on
            if empt == true
                j=j+1;
            else    
                shufAnimal = [];
                for k = 1:Iter
                    shufAnimal(k,:) = events_averg.shuf.(fn{ii}).(fn1{j})(k).Mean;
                end
                events_averg.shufCells.(fn{ii}).(fn1{j}) = shufAnimal;
            end
        end    
end

% combine all cells together
fn1=fieldnames(events_averg.shufCells.g5);
    for j=1:numel(fn1)
        events_averg.shufCells.all.(fn1{j}) = cat(2,events_averg.shufCells.g5.(fn1{j}),events_averg.shufCells.g4.(fn1{j}),events_averg.shufCells.g2.(fn1{j}),events_averg.shufCells.g12.(fn1{j}));
    end 
% get 5 and 95 percentile values
fn = fieldnames(events_averg.shufCells.all);
for ii = 1:numel(fn)
    for k = 1:115
       perc95(ii,k) = prctile(events_averg.shufCells.all.(fn{ii})(:,k), 95, 1);
       perc5(ii,k) =  prctile(events_averg.shufCells.all.(fn{ii})(:,k), 5, 1);
    end   
end    

% only animals with run  and er(cons) cells
fn1=fieldnames(events_averg.shufCells.g12);
    for j=1:numel(fn1)
        events_averg.shufCells.all.(fn1{j}) = cat(2,events_averg.shufCells.g2.(fn1{j}),events_averg.shufCells.g12.(fn1{j}));
    end 
% get 5 and 95 percentile values 
for ii = 20:numel(fn)
    for k = 1:66 % enter run
       perc95_run(ii-19,k) = prctile(events_averg.shufCells.all.(fn{ii})(:,k), 95, 1);
       perc5_run(ii-19,k) =  prctile(events_averg.shufCells.all.(fn{ii})(:,k), 5, 1);
    end
end  
ze=zeros(2,49);
perc95_run = cat(2,ze, perc95_run); perc95 = cat(1,perc95, perc95_run);% concatenate run and er for cells where that information exists
perc5_run = cat(2,ze, perc5_run);perc5 = cat(1,perc5, perc5_run);

fn = fieldnames(events_averg.shufCells.all);% make a table with correct row descriptions
Perc5 = array2table(perc5,"RowNames",fn); 
Perc95= array2table(perc95,"RowNames",fn);

save("Perc95_s2p_100", "Perc95");
save("Perc5_s2p_100", "Perc5");

events_averg_shuf = events_averg;
