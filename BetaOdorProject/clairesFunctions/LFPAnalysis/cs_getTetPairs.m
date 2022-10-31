clear
%----- Params -----%
animals = {'CS41','CS42','CS44','CS31','CS33','CS34','CS35','CS39'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS
%animals = {'CS34'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS

regions = {'CA1-PFC','CA1-OB','PFC-OB'};
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
   

%saves tet pairs separately across epochs. only consider unique tet pairs
%for whole animal, not separate across days/epochs

totalpairs = 0;
for n = 1:size(psf,2)
pairlists = psf(n).eegdata{1};
pairlist = cell2mat(pairlists');
pairs = unique(pairlist,'rows');
numpairs = size(pairs,1);
totalpairs = totalpairs+numpairs;
end

tetpairs.(['r',num2str(r)]) = totalpairs;
tetpairs.regions = regions;

end

save([topDir,'AnalysesAcrossAnimals\totaltetpairs'],'tetpairs');