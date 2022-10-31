% 
%----- Params -----%
 %animals = {'CS31','CS33','CS34','CS35','CS39','CS41','CS42','CS44'};
 animals = {'CS41','CS42'};
%regions = {'CA1','PFC','OB'};
regions = {'TC'};
 %regions = {'TC'};
%animals = {'CS41','CS42'};
[dataDir, figDir] = cs_setPaths(); 
dataDir = [dataDir,'AnalysesAcrossAnimals\'];

win = [0 3]; savetag = 'run';
inflection = 0.8;

% %win = [0 1.5]; savetag = 'pre_';
fpass = [0 20]; freqband = 'floor';
movingwin = [1000 20]/1000;

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
reforgnd = 'gnd'; %refstr = 'ref_';

createbaselinespecs = 0; 
if createbaselinespecs == 1
    cs_createBaselineSpecs(animals, regions, 1, fpass, movingwin, tapers, 'odorplace') 
end

eegspecfile = ['eeg',reforgnd,'spec',freqband];


varargstr = {'gnd', 1, 'win', win, 'tapers',tapers, 'movingwin', movingwin,'fpass',fpass,'inflection',inflection};


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
    case 'TC'
        tetfilter = '(isequal($area,''TC''))';
        
end

%----- Iterator -----%
   
iterator = 'cs_eeganal'; 
    
%----- Filter creation -----%
    
psf = createfilter('animal',animals,'epochs',runepochfilter,'eegtetrodes', tetfilter, 'iterator', iterator);

% psf.epochs{1, 1}= [3,2; 3,5];  
    
%----- Set analysis function -----%
    
out = setfilterfunction(psf, 'cs_calculatePSD',{'eeg', eegspecfile, Triggers},varargstr);
    
%----- Run analysis -----%
out_all = runfilter(out);

    

disp(['Done with ',region]);




data_sessions_pre = [];
data_sessions_post = [];

data_trials_pre = [];
data_trials_post = [];

data_sessions_pre_norm = [];
data_sessions_post_norm = [];

data_trials_pre_norm = [];
data_trials_post_norm = [];
for a = 1:length(animals)
    data = [out_all(a).output{1,1}];
   
    inds = vertcat(data.index);
    uniqueinds = unique(inds(:,1:2),'rows');
    
    %separate data by session
    for i = 1:size(uniqueinds,1)
        ind = uniqueinds(i,:);
        datainds = find(ismember(inds(:,1:2),ind,'rows'));
        epdata = data(datainds);
        
        %% raw data
        data_pre = vertcat(epdata.psd_pre);
        data_post = vertcat(epdata.psd_post);
        
        %combine across tetrodes
        data_pre = mean(data_pre,1);
        data_post = mean(data_post,1);

        data_sessions_pre = [data_sessions_pre;data_pre];
        data_sessions_post = [data_sessions_post;data_post];
        
        %combine across tetrodes
        data_pre = cat(3,epdata.psd_all_pre);
        data_pre = mean(data_pre,3)';
        data_post = cat(3,epdata.psd_all_post);
        data_post = mean(data_post,3)';
    
    
        data_trials_pre = [data_trials_pre;data_pre];
        data_trials_post = [data_trials_post;data_post];
        
        %% normalized data
        data_pre = vertcat(epdata.psd_pre_norm);
        data_post = vertcat(epdata.psd_post_norm);
        
        %combine across tetrodes
        data_pre = mean(data_pre,1);
        data_post = mean(data_post,1);

        data_sessions_pre_norm = [data_sessions_pre_norm;data_pre];
        data_sessions_post_norm = [data_sessions_post_norm;data_post];
        
        %combine across tetrodes
        data_pre = cat(3,epdata.psd_all_pre_norm);
        data_pre = mean(data_pre,3)';
        data_post = cat(3,epdata.psd_all_post_norm);
        data_post = mean(data_post,3)';
    
    
        data_trials_pre_norm = [data_trials_pre_norm;data_pre];
        data_trials_post_norm = [data_trials_post_norm;data_post];
    end
        
end

%average for each day


PSDdata.sessions_pre = data_sessions_pre;
PSDdata.sessions_post = data_sessions_post;
PSDdata.trials_pre = data_trials_pre;
PSDdata.trials_post = data_trials_post;

PSDdata.sessions_pre_norm = data_sessions_pre_norm;
PSDdata.sessions_post_norm = data_sessions_post_norm;
PSDdata.trials_pre_norm = data_trials_pre_norm;
PSDdata.trials_post_norm = data_trials_post_norm;


PSDdata.fpass = fpass;
PSDdata.region = region;
PSDdata.win = win;
PSDdata.tapers = tapers;
PSDdata.movingwin = movingwin;

filename = ['PSDdata_',savetag, region];
    
save([dataDir,filename],'PSDdata');

end