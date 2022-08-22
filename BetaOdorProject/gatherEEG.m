%% GatherEEGdata (just gather the data, no analysis, newer)


eegDir='I:\BrandeisData\SymanskiData';


% set region parameters here!!
regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')];

%% what does RR look like anyways

%temp=load('E:\ClaireData\CS31_direct\EEG\CS31resp02-04-20');

load('respfilter');
load('betafilter');
load('ripplefilter');

%% Initial code block to index the file from which I'll pullRR data

% load eeg data,generate a contiuous datastream (raw)
%if ~isfield(SuperRat,'CA1beta')

% the gist is pull the eeg from whichever tetrode has the most cells in
% each region
% the alternative is to pull whichver tetrode has task responsive cells in
% each region
    
for i=1:length(SuperRat)
    
    % need to use a local beta metric
    for j=1:length(regions)
        % get units in that area that arent MUA
        myunits=find(contains({SuperRat(i).units.area},regions{j}));
        if ~isempty(myunits)
            [~,LFPtet]=max(accumarray([SuperRat(i).units(myunits).tet]',1));
        else % if no units, grab the first one in that region
            oktets=SuperRat(i).tetinfo(cellfun(@(a) ~isempty(a), {SuperRat(i).tetinfo.area}));
            LFPtet=oktets(find(contains({oktets.area},regions{j}),1,'first')).tetnum;
        end
        if isempty(LFPtet) || SuperRat(i).longTrack==0
            continue;
        end
        

        % generate filename here
        [~,name]=system('hostname');
        if strcmpi(strtrim(name),'desktop-t3q247t')
            ratname=sprintf('%sExpt',SuperRat(i).name);
            allLFPfiles=dir(fullfile(eegDir,ratname,'**/*.mat'));
        else
            ratname=sprintf('%s_direct',SuperRat(i).name);
            lfpdir='EEG';
            allLFPfiles=dir(fullfile(eegDir,ratname,lfpdir));
        end
        % get lfp files from this tetrode and today
        sessname=sprintf('%seeg%02d',SuperRat(i).name,SuperRat(i).daynum);
        todayfiles=contains({allLFPfiles.name},sessname);

        % skip middle number (epoch) and go straight to last number
        % (tetrode)
        mytet=sprintf('%02d.mat',LFPtet);
        tetfiles=contains({allLFPfiles.name},mytet);
        loadfiles=allLFPfiles(todayfiles & tetfiles);
        % sort the files in ascending order
        [~,index] = sortrows({loadfiles.name}.'); loadfiles = loadfiles(index); clear index
        respcont=[]; betacont=[]; ripplecont=[];
        clear lfpData
        for k=1:length(loadfiles)
            lfpBit=load(fullfile(loadfiles(k).folder,loadfiles(k).name));
            tempstruct=lfpBit.eeg{SuperRat(i).daynum}{k}{LFPtet};
            tempstruct.filename=loadfiles(k).name;
            lfpData(k)=tempstruct;
            if tempstruct.data_voltage_inverted==1
                tempstruct.data=tempstruct.data.*-1;
                tempstruct.data_voltage_inverted=0;
            end
            % now filter those data!
            respfiltered=filtereeg2(tempstruct,respfilter);
            betafiltered=filtereeg2(tempstruct,betafilter);
            ripplefiltered=filtereeg2(tempstruct,ripplefilter);
            % generate continuous data (ts, amp, instaphase, envelope)
            respcont=[respcont; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(respfiltered.data)];
            betacont=[betacont; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(betafiltered.data)];
            ripplecont=[ripplecont; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(ripplefiltered.data)];
        end
        
        clear tempstruct;
        
        
        % save struct out
        SuperRat(i).([regions{j} lfpdir])=lfpData;
        % save continuous phase data out
        SuperRat(i).([regions{j} 'resp'])=sortrows(respcont,1); % sort epochs by time!
        SuperRat(i).([regions{j} 'beta'])=sortrows(betacont,1); % sort epochs by time!
        SuperRat(i).([regions{j} 'ripple'])=sortrows(betacont,1); % sort epochs by time!


    end
    fprintf('session number %d, animal %s day %d done\n',i,SuperRat(i).name,SuperRat(i).daynum);
    
end
