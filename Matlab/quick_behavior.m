% code to quickly analyze behavior data only
% input:  events file from pi
% output: cleaned and organized metadata bh_data
% 05.09.2023

clearvars -except bh_data; 
%% import data
animal = 4;


if animal == 2 
    RFID = '328340226232';
elseif animal == 4
    RFID = '335490167178';
elseif animal == 5    
    RFID ='3354903011';
elseif animal == 6    
    RFID ='8501114126';
elseif animal == 7    
    RFID ='8502199200';
elseif animal == 9    
    RFID ='8501181185';
elseif animal == 10    
    RFID ='3354916686';
elseif animal == 11    
    RFID ='8500188177';
elseif animal == 12    
    RFID ='3354909574';
elseif animal == 14    
    RFID ='8500142131';  
elseif animal == 15    
    RFID ='1251667485196';    
end    
events_path=['E:\Exp_2_DM_exploration\PFC-LH\' RFID '_events.csv'];

level = wildcardPattern + "\";pat = asManyOfPattern(level);
day = extractAfter(topLevelFolder,pat);day2 = strrep(day,'_','-');% get the day
chop_from=datetime([day2 ' 09:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');%chop the day
chop_to=datetime([day2 ' 22:00:00.000'],'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x<chop_from),:)=[];%chops events to match the day
x=datetime(events.(1),'InputFormat','yyyy-MM-dd HH:mm:ss.SSS');
events(find(x>chop_to),:)=[];