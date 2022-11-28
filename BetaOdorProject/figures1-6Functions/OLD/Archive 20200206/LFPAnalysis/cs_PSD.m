% 
%----- Params -----%
animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};  %CHANGE _direct PATH IN 'animaldef.m' FOR DIFFERENT COMPUTERS
regions = {'OB','CA1','PFC'};
[dataDir, figDir] = cs_setPaths(); 
dataDir = [dataDir,'AnalysesAcrossAnimals\'];

savefile = 1;

Triggers = 'odorTriggers';
 
switch Triggers
    case 'odorTriggers'
        trigstr = '';
end

runepochfilter = 'isequal($environment, ''odorplace'')';

novelodors = 0;
if novelodors == 1
    runepochfilter = 'isequal($environment, ''novel1'')';
    trigstr = '_novel';
end

trigtypes = {'allTriggers'};

tapers = [1 1];
reforgnd = 'ref'; refstr = 'ref_';
if strcmp(reforgnd,'gnd')
    do_wrtgnd = 1;
else
    do_wrtgnd = 0;
end

win = [0 1.5];
fpass = [0 55];
movingwin = [1000 20]/1000;

varargstr = {'gnd', 1, 'win', win, 'tapers',tapers, 'movingwin', movingwin,'fpass',fpass};


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1:length(regions)
    region = regions{r};

%----- Select Data -----%

%runepochfilter = 'isequal($environment, ''laseroff'')';

switch region
    case 'PFC'
        tetfilter = 'isequal($area, ''PFC'')&& $numcells > 0';
    case 'CA1'
        tetfilter = 'isequal($area, ''CA1'') && $numcells > 0';
    case 'OB'
        tetfilter = '(isequal($area, ''OB''))';
        
end

%----- Iterator -----%
   
iterator = 'cs_eeganal'; 
    
%----- Filter creation -----%
    
psf = createfilter('animal',animals,'epochs',runepochfilter,'eegtetrodes', tetfilter, 'iterator', iterator);

% psf.epochs{1, 1}= [3,2; 3,5];  
    
%----- Set analysis function -----%
    
out = setfilterfunction(psf, 'cs_calculatePSD',{'eeg', Triggers},varargstr);
    
%----- Run analysis -----%
out_all = runfilter(out);

    

disp(['Done with ',region]);




newdata = [];
for a = 1:length(animals)
    data = [out_all(a).output{1,1}];
    data = data(find(~cellfun(@isempty,{data.psd})));
    data = cat(3, data.psd);

    data = mean(data,3)';
    newdata(:,a) = data;
end

PSDdata.data = newdata;
PSDdata.fpass = fpass;
PSDdata.region = region;
PSDdata.win = win;
PSDdata.tapers = tapers;
PSDdata.movingwin = movingwin;

filename = ['PSDdata_', region];
    
save([dataDir,filename],'PSDdata');

psd = mean(newdata,2);
freqs = fpass(2)/length(psd):fpass(2)/length(psd):fpass(2);

plot(freqs,psd);
ylabel('Power');
xlabel('Frequency');

 figfile = [figDir,'EEG\PSD_',region];
        
        print('-djpeg', figfile);
        print('-dpdf', figfile);

end