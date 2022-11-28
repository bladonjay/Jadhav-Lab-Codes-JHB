clear
%----- Params -----%
animals = {'CS41','CS42','CS44','CS31','CS33','CS34','CS35','CS39'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS
%animals = {'CS34'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS

regions = {'CA1','PFC','OB','TC'};
%regions = {'PFC-OB','CA1-OB'};

topDir = cs_setPaths();
dataDir = [topDir,'AnalysesAcrossAnimals\'];
%---- Triggers -----%
trigtypes = {'allTriggers'};
freqband = 'beta';
window = [0.5 1];

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
iterator = 'cs_eeganal'; 
 

 
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:length(regions)
    region = regions{r};
    
%----- Select Data -----%

switch region
    case 'CA1'
        tetfilter = '(strcmp($area, ''CA1'' ))';
        %tetfilter = {'(strcmp($area, ''CA1'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''hpcRef'' ) && strcmp($descrip2, ''betatet'' ))',...
           % '(strcmp($area, ''PFC'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''pfcRef'' ) && strcmp($descrip2, ''betatet'' ))'};
    case 'PFC'
        tetfilter = '(strcmp($area, ''PFC''))';
        %tetfilter = {'(strcmp($area, ''PFC'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''pfcRef'' ) && strcmp($descrip2, ''betatet'' ))',...
        %    '(strcmp($area, ''OB'') && strcmp($descrip2,''betatet''))'};
    case 'OB'
        tetfilter = '(strcmp($area, ''OB''))';
        %tetfilter = {'(strcmp($area, ''CA1'' ) && strcmp($descrip2, ''betatet'' )) || (strcmp($descrip, ''hpcRef'' ) && strcmp($descrip2, ''betatet'' ))',...
         %   '(strcmp($area, ''OB'') && strcmp($descrip2,''betatet''))'};
    case 'TC'
        tetfilter = '(strcmp($area, ''TC'' ))';
end


    
%----- Filter creation -----%
    
psf = createfilter('animal',animals,'epochs',runepochfilter,'eegtetrodes',tetfilter, 'iterator', iterator);
   

%saves tet pairs separately across epochs. only consider unique tet pairs
%for whole animal, not separate across days/epochs

totaltets = 0;
for n = 1:size(psf,2)
tetlists = psf(n).eegdata{1};
list = cell2mat(tetlists');
tets = unique(list,'rows');
numtets = size(tets,1);
totaltets = totaltets+numtets;
end

tetlist.(['r',num2str(r)]) = totaltets;
tetlist.regions = regions;

end

save([topDir,'AnalysesAcrossAnimals\totaltets'],'tetlist');