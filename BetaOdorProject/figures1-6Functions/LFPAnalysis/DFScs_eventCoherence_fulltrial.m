
%----- Params -----%
%animals = {'CS41','CS42','CS44','CS31','CS33','CS34','CS35','CS39'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS
animals = {'CS34'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS

regions = {'CA1-PFC'};
%regions = {'PFC-OB','CA1-OB'};

topDir = cs_setPaths();
dataDir = [topDir,animals{1},'Expt\',animals{1},'_direct\'];
%---- Triggers -----%
trigtypes = {'allTriggers'};
freqband = 'low';
window = [1 5];
fpass = [0 40];

%Triggers = 'stemTimes';
Triggers = 'odorTriggers';
 if strcmp(Triggers,'odorTriggers')
     trigstr = '';
 else 
     trigstr = [Triggers,'_'];
 end
 
 %---- Filters -----%
 runepochfilter = 'isequal($environment, ''odorplace'')';
 odorstr = '';
%Novel odor: 
% odorstr = 'novel_';
% runepochfilter = 'isequal($environment, ''novelodor'') || isequal($environment, ''novelodor2'')';
 
 %----- Iterator -----% 
iterator = 'tetpaireeganal'; 
 

 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(regions)
    region = regions{r};
    
%----- Select Data -----%

switch region
    case 'CA1-PFC'
        tetfilter = {'(strcmp($area, ''CA1'' )) && $numcells >0','(strcmp($area, ''PFC'' )) && $numcells >0'};
        %tetfilter = {'(strcmp($area, ''CA1'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''hpcRef'' ) && strcmp($descrip2, ''betatet'' ))',...
           % '(strcmp($area, ''PFC'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''pfcRef'' ) && strcmp($descrip2, ''betatet'' ))'};
    case 'PFC-OB'
        tetfilter = {'(strcmp($area, ''PFC'')) && $numcells >0', '(strcmp($area, ''OB'' ))'};
        %tetfilter = {'(strcmp($area, ''PFC'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''pfcRef'' ) && strcmp($descrip2, ''betatet'' ))',...
        %    '(strcmp($area, ''OB'') && strcmp($descrip2,''betatet''))'};
    case 'CA1-OB'
        tetfilter = {'(strcmp($area, ''CA1'' )) && $numcells >0','(strcmp($area, ''OB''))'};
        %tetfilter = {'(strcmp($area, ''CA1'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''hpcRef'' ) && strcmp($descrip2, ''betatet'' ))',...
         %   '(strcmp($area, ''OB'') && strcmp($descrip2,''betatet''))'};
    case 'V1-OB'
        tetfilter = {'(strcmp($area, ''V1'' ))','(strcmp($area, ''OB''))'};
end


    
%----- Filter creation -----%
    
psf = createfilter('animal',animals,'epochs',runepochfilter,'eegtetrodepairs',tetfilter, 'iterator', iterator);
   
    
%----- Set analysis function -----%
    
out = setfilterfunction(psf, 'DFAcs_eventcoherence',{'eeg',Triggers},'trigtypes',trigtypes,'window',window, 'fpass',fpass);
    
%----- Run analysis -----%
out_all = runfilter(out);

    

disp('Done with all');

%---- Save file with params -----%
coherenceData.data = out_all;
coherenceData.animals = animals;
coherenceData.region = region;
coherenceData.trigtypes = trigtypes;
coherenceData.fpass = fpass;

filename = [animals{1},'coherence', region,'_',freqband,'_fulltrial.mat'];

save([dataDir,filename],'coherenceData');
%save(['E:\AnalysesAcrossAnimals\',filename],'coherenceData');

%%%% run cs_plotCoherence to plot %%%% 
%cs_plotCoherence(['E:\AnalysesAcrossAnimals\',filename])
cs_plotCoherence([dataDir,filename])
close all

end


