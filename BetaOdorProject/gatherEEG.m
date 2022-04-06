%% GatherEEGdata


% set region parameters here!!
regions={'PFC','CA1','OB'};
colors=[rgbcolormap('DarkAquamarine'); rgbcolormap('LightCoral'); rgbcolormap('DarkOrange')];
rhythmcolors=[rgbcolormap('navy'); rgbcolormap('DeepPink')];

%% what does RR look like anyways

%temp=load('E:\ClaireData\CS31_direct\EEG\CS31resp02-04-20');

load('respfilter');
load('betafilter');
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
        myunits=find(contains({SuperRat(i).units.area},regions{j}) & ...
            ~contains({SuperRat(i).units.tag},'mua'));
        %myunits=find(contains({SuperRat(i).units.area},regions{j}));
        if ~isempty(myunits)
            [~,LFPtet]=max(accumarray([SuperRat(i).units(myunits).tet]',1));
        else % if no units, grab the first one in that region
            oktets=SuperRat(i).tetinfo(cellfun(@(a) ~isempty(a), {SuperRat(i).tetinfo.area}));
            LFPtet=oktets(find(contains({oktets.area},regions{j}),1,'first')).tetnum;
        end
        if isempty(LFPtet)
            continue;
        end
        
        % generate filename here
        ratname=sprintf('%s_direct',SuperRat(i).name);
        lfpdir='EEG';
        allLFPfiles=dir(fullfile('E:\ClaireData',ratname,lfpdir));
        
        % get lfp files from this tetrode and today
        sessname=sprintf('%seeg%02d',SuperRat(i).name,SuperRat(i).daynum);
        todayfiles=contains({allLFPfiles.name},sessname);
        mytet=sprintf('%02d.mat',LFPtet);
        tetfiles=contains({allLFPfiles.name},mytet);
        loadfiles=allLFPfiles(todayfiles & tetfiles);
        % sort the files in ascending order
        [~,index] = sortrows({loadfiles.name}.'); loadfiles = loadfiles(index); clear index
        contdata=[]; contdata2=[];
        clear lfpData
        for k=1:length(loadfiles)
            lfpBit=load(fullfile(loadfiles(k).folder,loadfiles(k).name));
            tempstruct=lfpBit.eeg{SuperRat(i).daynum}{k}{LFPtet};
            tempstruct.tet=LFPtet; tempstruct.filename=loadfiles(k).name;
            lfpData(k)=tempstruct;
            filtered=filtereeg2(tempstruct,respfilter);
            filtered2=filtereeg2(tempstruct,betafilter);

            % generate continuous data (ts, amp, instaphase, envelope)
            contdata=[contdata; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(filtered.data)];
            contdata2=[contdata2; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(filtered2.data)];
        end
        % now filter those data!
        clear tempstruct;
        
        
        % save struct out
        SuperRat(i).([regions{j} lfpdir])=lfpData;
        % save continuous phase data out
        SuperRat(i).([regions{j} 'resp'])=sortrows(contdata,1); % sort epochs by time!
        SuperRat(i).([regions{j} 'beta'])=sortrows(contdata,1); % sort epochs by time!

    end
    fprintf('session number %d, animal %s day %d done\n',i,SuperRat(i).name,SuperRat(i).daynum);
    
end
