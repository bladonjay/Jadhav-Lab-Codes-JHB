%
%
%
%%


% so we want to calculate the phase amp cfc during the odor period.
% we can either normalize it to random other periods (the same times before
% sampling) or just bootstrap it by shuffling.


%%
%{
 % probably need a different filter here...
        % here i'm butterworthing and then hilberting.... could fiir it too
        h=waitbar(0,'gathering wavelets');
        LECspectam=abs(cwt(LECdata,bases,'cmor16-.5')); % the spectrogram for ALLLLL
        HPCspectam=abs(cwt(HPCdata,bases,'cmor16-.5'));
        
        waitbar(0,h,sprintf('Finished Wavelets, now comodulogram'));
        % then is the fast freq
        for fl=1:size(LoFreqs,1)
            % grab the lo frequencies
            [LoLphase,LoLamp]=GetLFPBand(LECdata,LFPts,LoFreqs(fl,:),1);
            [LoHphase,LoHamp]=GetLFPBand(HPCdata,LFPts,LoFreqs(fl,:),1);
            
            for fh=1:length(HiFreqs)
                
                % here you can trialwise boot by shuffling the trials and
                % correlating- if they're all the same length the
                % concatenation should line up the same
                % so to input, we'll need two transforms on each cell: 1.
                % the lo amp, and 2 the high phase that way we can only
                % pull the high phases that are also lo amp good
                loamps=cell2mat(cellfun(@(a) LoLamp(a), evinds, 'UniformOutput', false))';
                lophases=cell2mat(cellfun(@(a) LoLphase(a), evinds, 'UniformOutput', false))';
                hiamps=cell2mat(cellfun(@(a) LECspectam(fh,a), evinds, 'UniformOutput', false)');
                [~,mu,sig]=zscore(LoLamp);
                [LLcomod{fl,fh}(ses),~,~,~,~,LLnull{ses}] = CalcPhaseAmpCFC(lophases,...
                    hiamps,'LoAmps',loamps,'Cutoff',mu+sig,'runstat',100);
                
                % maybe gather the mean and std from this null, or cat all the nulls to get
                % the true dist, and then patch the real from the means
                [HHcomod{fl,fh}(ses),~,~,~,~,HHnull{ses}] = CalcPhaseAmpCFC(LoHphase(intersect(allevinds,find(OKphaseh)))',...
                    HPCspectam(fh,intersect(allevinds,find(OKphaseh)))','runstat',100);
            end
        end
%}
        
%%

% so we're first doing cross frequency coupling, not phase amp coupling


%{
[C,phi,S12,S1,S2,f,confC,phistd,Cerr]=coherencyc(data1,data2,params)
data1=raw first lfp
data2=raw second lfp

params includes:
params=struct('tapers',[3 5],'fs',1500,'pad',1,'Fpass',[2 40]);


%}

% the biggest choice here is which HPC lfp to use.  I am inclined to use
% the tet with the most odor responsive units or the tet with the highest
% beta to gamma power (15-30) hz to 30-45 hz ratio
%% lets see if we can get some beta to gamma ratios in tetrodes
% if there is a good spread, i'll try to relate that to place fields and
% odor fields
% per igarashi et al extended data figure 5


for i=1:length(SuperRat)
    tic
    % get behavior epochs
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10;
    odorTimes=odorTimes(diff(odorTimes')>.5,:);
    trialCorrect=trialCorrect(diff(odorTimes')>.5);
    
    % find HPC tetrodes
    myunits=find(contains({SuperRat(i).units.area},'CA1'));
    LFPtet=find(accumarray([SuperRat(i).units(myunits).tet]',1)>0);
    % generate filename here
    ratname=sprintf('%s_direct',SuperRat(i).name);
    lfpdir='EEG';
    allLFPfiles=dir(fullfile('E:\ClaireData',ratname,lfpdir));
    
    % get lfp files from today
    sessname=sprintf('%seeg%02d',SuperRat(i).name,SuperRat(i).daynum);
    todayfiles=allLFPfiles(contains({allLFPfiles.name},sessname));
    for j=1:length(todayfiles)
        todayfiles(j).tetrode=str2double(todayfiles(j).name(end-5:end-4));
    end
    % for each tetrode, accumulate all the data, get timestamps,
    % just use odor period, calculate your plots, and get ratio
    for j=1:length(LFPtet)
        % gather all files for this tetrode
        
        loadfiles=todayfiles([todayfiles.tetrode]==LFPtet(j));
        contdata=[];
        for k=1:length(loadfiles)
            lfpBit=load(fullfile(loadfiles(k).folder,loadfiles(k).name));
            tempstruct=lfpBit.eeg{SuperRat(i).daynum}{k}{LFPtet(j)};
            tempstruct.tet=LFPtet(j); tempstruct.filename=loadfiles(k).name;
            %filtered=filtereeg2(tempstruct,respfilter);
            % generate continuous data (ts, amp, instaphase, envelope)
            contdata=[contdata; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(tempstruct.data)];
        end
        
        % we'll use from 15-45, but split it in half
        params=struct('tapers',[3 5],'Fs',1500,'pad',1,'fpass',[20 40]);
         
         eegData=[];
         for k=1:length(odorTimes)
            eegData(k,:)=contdata(contdata(:,1)>=(odorTimes(k,2)-1) & contdata(:,1)<odorTimes(k,2),2);
         end
        [S,f]=mtspectrumc(eegData',params);
        SuperRat(i).betagamma(LFPtet(j))=mean(mean(S(f<27,:)))/mean(mean(S(f>33,:)));
    end
    SuperRat(i).betagamma(SuperRat(i).betagamma==0)=nan;
    fprintf('Session %d, %s %d done %d seconds \n',i, SuperRat(i).name,SuperRat(i).daynum,toc);
end


allbetagamma=cellfun(@(a) a, {SuperRat.betagamma},'UniformOutput',false);
allbetagamma2=cellfun(@(a) max(a), {SuperRat.betagamma},'UniformOutput',false);

histogram(cell2mat(allbetagamma),0:.5:6);
hold on;
histogram(cell2mat(allbetagamma2),0:.5:6);
% i think what could come next is two things:
% 1. run the cross coherence across the tetrodes, see of betagamma coupling
% pairs with beta cross coherence (may need to consider the amplitude
% confound though
% 2. cross correlate this with the number of odor coding cells on that HPC
% tetrode


%%





%% beta coherence comparing correct and incorrect trials

% currently using tetrode with highest beta to gamma power in HPC

params=struct('tapers',[3 5],'Fs',1500,'pad',1,'fpass',[0 60]);

% want to do this separately for each trial I think
for i=1:length(SuperRat)
    
    try
    eegTime=SuperRat(i).PFCbeta(:,1);
    
    ratname=sprintf('%s_direct',SuperRat(i).name);
    allLFPfiles=dir(fullfile('E:\ClaireData',ratname,'EEG'));
    % get lfp files from today
    sessname=sprintf('%seeg%02d',SuperRat(i).name,SuperRat(i).daynum);
    todayfiles=allLFPfiles(contains({allLFPfiles.name},sessname));
    for j=1:length(todayfiles)
        todayfiles(j).tetrode=str2double(todayfiles(j).name(end-5:end-4));
    end
    [~,HPCtet]=max(SuperRat(i).betagamma);
    loadfiles=todayfiles([todayfiles.tetrode]==HPCtet);
    contdata=[];
    for k=1:length(loadfiles)
        lfpBit=load(fullfile(loadfiles(k).folder,loadfiles(k).name));
        tempstruct=lfpBit.eeg{SuperRat(i).daynum}{k}{HPCtet};
        %tempstruct.tet=LFPtet; tempstruct.filename=loadfiles(k).name;
        %filtered=filtereeg2(tempstruct,respfilter);
        % generate continuous data (ts, amp, instaphase, envelope)
        contdata=[contdata; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(tempstruct.data)];
    end
    contdata=sortrows(contdata,1);
    HPCeeg=zscore(contdata(:,2));
    PFCeeg=zscore(cell2mat({SuperRat(i).PFCEEG.data}'));
    % now get our trial blocks
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10;
    odorTimes=odorTimes(diff(odorTimes')>.5,:);
    trialCorrect=trialCorrect(diff(odorTimes')>.5);
    C=[]; S1=[]; S2=[];
    for k=1:length(odorTimes)
        okinds=eegTime>=odorTimes(k,2)-0.25 & eegTime<=odorTimes(k,2)+.25;
        [C(k,:),phi,S12,S1(k,:),S2(k,:),f]=coherencyc(HPCeeg(okinds),PFCeeg(okinds),params);
    end
    figure
    errorbar(f,nanmean(C(trialCorrect==1,:)),SEM(C(trialCorrect==1,:),1));
    %plot(f,zscore(C(trialCorrect==1,:),1,2),'b');
    hold on;
    errorbar(f,nanmean(C(trialCorrect==0,:)),SEM(C(trialCorrect==0,:),1));
    %plot(f,zscore(C(trialCorrect~=1,:),1,2),'r');
    title(sprintf('PFC HPC coherence Session %d from %s',SuperRat(i).daynum,SuperRat(i).name));
    
    drawnow;
    catch
        fprintf('Session %d with %s didnt work \n',SuperRat(i).daynum,SuperRat(i).name)
    end
end

% Tried this with HPC tetrode with most cells vs PFC tet with most cells
% results did not look good

% then tried this with HPC tetrode with highest beta-gamma ratio with PFC
% tet with most cells, that did not look good either

% I think we need more tets with distal CA1, I have a feeling most of these
% recordings do not have distal tets and therefore are not useful for this
% analysis
%%

params=struct('tapers',[3 5],'Fs',1500,'pad',1,'fpass',[0 60]);

% want to do this separately for each trial I think
for i=1:length(SuperRat)
    
    try
    eegTime=SuperRat(i).PFCbeta(:,1);
    
    ratname=sprintf('%s_direct',SuperRat(i).name);
    allLFPfiles=dir(fullfile('E:\ClaireData',ratname,'EEG'));
    % get lfp files from today
    sessname=sprintf('%seeg%02d',SuperRat(i).name,SuperRat(i).daynum);
    todayfiles=allLFPfiles(contains({allLFPfiles.name},sessname));
    for j=1:length(todayfiles)
        todayfiles(j).tetrode=str2double(todayfiles(j).name(end-5:end-4));
    end
    [~,HPCtet]=max(SuperRat(i).betagamma);
    loadfiles=todayfiles([todayfiles.tetrode]==HPCtet);
    contdata=[];
    for k=1:length(loadfiles)
        lfpBit=load(fullfile(loadfiles(k).folder,loadfiles(k).name));
        tempstruct=lfpBit.eeg{SuperRat(i).daynum}{k}{HPCtet};
        %tempstruct.tet=LFPtet; tempstruct.filename=loadfiles(k).name;
        %filtered=filtereeg2(tempstruct,respfilter);
        % generate continuous data (ts, amp, instaphase, envelope)
        contdata=[contdata; (tempstruct.starttime:(1/tempstruct.samprate):tempstruct.endtime)' double(tempstruct.data)];
    end
    contdata=sortrows(contdata,1);
    HPCeeg=zscore(contdata(:,2));
    PFCeeg=zscore(cell2mat({SuperRat(i).PFCEEG.data}'));
    % now get our trial blocks
    odorTimes=[SuperRat(i).trialdata.sniffstart SuperRat(i).trialdata.sniffend];
    trialCorrect=SuperRat(i).trialdata.CorrIncorr10;
    odorTimes=odorTimes(diff(odorTimes')>.5,:);
    trialCorrect=trialCorrect(diff(odorTimes')>.5);
    C=[]; S1=[]; S2=[];
    for k=1:length(odorTimes)
        okinds=eegTime>=odorTimes(k,2)-0.25 & eegTime<=odorTimes(k,2)+.25;
        [C(k,:),phi,S12,S1(k,:),S2(k,:),f]=coherencyc(HPCeeg(okinds),PFCeeg(okinds),params);
    end
    figure
    errorbar(f,nanmean(C(trialCorrect==1,:)),SEM(C(trialCorrect==1,:),1));
    %plot(f,zscore(C(trialCorrect==1,:),1,2),'b');
    hold on;
    errorbar(f,nanmean(C(trialCorrect==0,:)),SEM(C(trialCorrect==0,:),1));
    %plot(f,zscore(C(trialCorrect~=1,:),1,2),'r');
    title(sprintf('PFC HPC coherence Session %d from %s',SuperRat(i).daynum,SuperRat(i).name));
    
    drawnow;
    catch
        fprintf('Session %d with %s didnt work \n',SuperRat(i).daynum,SuperRat(i).name)
    end
end