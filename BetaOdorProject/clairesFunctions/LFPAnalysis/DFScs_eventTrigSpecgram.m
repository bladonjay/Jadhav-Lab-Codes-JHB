%% CHECK ALL PARAMS EACH TIME SCRIPT IS RUN. 
%% ----- Params -----%
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS
%animals = {'CS41','CS42'};
%animals = {'CS39','CS41','CS42'};

regions = {'OB','CA1','PFC'};
%regions = {'CA1'};
%regions = {'riptet'};
[dataDir, figDir] = cs_setPaths(); 
dataDir = [dataDir,'AnalysesAcrossAnimals\SpecgramData\'];


savefile = 1;

%Triggers = 'odorTriggers';
Triggers = 'odorOffset';
%Triggers = 'laserTriggers';
%Triggers = 'stemTimes';
%Triggers = 'rewardTimes';
%Triggers = 'nperrors';
 
switch Triggers
    case 'odorTriggers'
        trigstr = '';
    case 'odorOffset'
        trigstr = '_odoroff';
    case 'stemTimes'
        trigstr = '_stem';
    case 'laserTriggers'
        trigstr = '_laser';
    case 'rewardTimes'
        trigstr = '_reward';
    case 'nperrors'
        trigstr = '_nperrors';
end

envfilt = 'odorplace';
runepochfilter = ['isequal($environment, ''',envfilt,''')'];
if strcmp(envfilt,'odorplace')
    envstr = '';
else
    envstr = envfilt;
end

novelodors = 0;
if novelodors == 1
    runepochfilter = 'isequal($environment, ''novel1'')';
    trigstr = '_novel';
end

%trigtypes = {'incorrectTriggers'}; trigstr = '_incorrect';
%trigtypes = {'correctTriggers'}; trigstr = '_correct';
trigtypes = {'allTriggers'};

tapers = [1 1];
reforgnd = 'gnd'; 
if strcmp(reforgnd,'gnd')
    do_wrtgnd = 1;
    refstr = 'ref_';
else
    do_wrtgnd = 0;
    refstr = '';
end
%reforgnd = 'ref'; refstr = 'ref_';
%freqband = 'high'; fpass = [0 300];  movingwin = [100 10]/1000; %mid = 0-100 Hz
%freqband = 'mid'; fpass = [50 100];  movingwin = [400 20]/1000; %mid = 0-100 Hz
freqband = 'low'; fpass = [0 40]; movingwin = [1000 20]/1000;  %low = 0-40 Hz
%freqband = 'floor'; fpass = [0 15]; movingwin = [2000 25]/1000; tapers = [1 1];%floor = 0-15 Hz
win = [1.5 0.5];

createbaselinespecs = 1; 
if createbaselinespecs == 1
    cs_createBaselineSpecs(animals, regions, do_wrtgnd, fpass, movingwin, tapers, envfilt) 
end

eegspecfile = ['eeg',reforgnd,'spec',freqband]; %(eg eegrefspeclow, eeggndspecmid, etc)

%create file strings
binsizestr = [num2str(movingwin(2)*1000),'msBins'];
winstring = [num2str(-win(1)*1000),'-',num2str(win(2)*1000),'ms'];
taperstring = ['tapers',num2str(tapers(1)),'-', num2str(tapers(2))];

filenametitle = ['eventTrigSpecgramData_'];
filenameparams = [trigstr,'_',freqband,'_', envstr,refstr, winstring,'_',taperstring,'_', binsizestr, '.mat'];

varargstr = {'trigtypes',trigtypes, 'gnd', 1, 'win', win, 'tapers',tapers, 'movingwin', movingwin,'fpass',fpass};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(regions)
    region = regions{r};
disp(['Doing ',region])
%----- Select Data -----%

%runepochfilter = 'isequal($environment, ''laseroff'')';

switch region
    case 'PFC'
        tetfilter = '(isequal($area, ''PFC''))';
    case 'CA1'
        tetfilter = '(isequal($area, ''CA1''))';
    case 'OB'
        tetfilter = '(isequal($area, ''OB''))';
    case 'V1'
        tetfilter = '(isequal($area, ''V1''))';
    case 'TC'
        tetfilter = '(isequal($area, ''TC''))';
    case 'riptet'
        tetfilter = '(strcmp($descrip, ''riptet''))';
        
end

%----- Iterator -----%
   
iterator = 'cs_eeganal'; 
    
%----- Filter creation -----%
    
psf = createfilter('animal',animals,'epochs',runepochfilter,'eegtetrodes', tetfilter, 'iterator', iterator);

% psf.epochs{1, 1}= [3,2; 3,5];  
    
%----- Set analysis function -----%
    
out = setfilterfunction(psf, 'DFAcs_eventTrigSpecgram',{'eeg', eegspecfile, Triggers},varargstr);
    
%----- Run analysis -----%
out_all = runfilter(out);

    

disp(['Done with ',region]);



for t = 1:length(trigtypes)
        trigtype = trigtypes{t};
        
        newdata = [];
        for a = 1:length(animals)
            data = [out_all(a).output{1,1}.(trigtype)];
            data = data(find(~cellfun(@isempty,{data.Smean})));
            data = cat(3, data.Smean);
            
            data = mean(data,3);
            newdata(:,:,a) = data;
        end
        
        MeanSpec = mean(newdata,3);
        eventTrigSpecgramData.(trigtype)= MeanSpec;
end

eventTrigSpecgramData.ReforGnd = reforgnd;
eventTrigSpecgramData.freqband = freqband;
eventTrigSpecgramData.trigtypes = trigtypes;
eventTrigSpecgramData.region = region;
eventTrigSpecgramData.win = win;
eventTrigSpecgramData.tapers = tapers;
eventTrigSpecgramData.movingwin = movingwin;

if savefile == 1
filename = [filenametitle, region, filenameparams];
    
save([dataDir,filename],'eventTrigSpecgramData');
end

[maxval, minval] = cs_plotSpecgram([dataDir,filename],'freqband',freqband);
maxvals(r) = maxval;
minvals(r) = minval;

end

% limits(1) = min(minvals);
% limits(2) = max(maxvals);
% 
% 
% for r = 1:length(regions)
%     figure(r)
%     caxis(limits) 
% end