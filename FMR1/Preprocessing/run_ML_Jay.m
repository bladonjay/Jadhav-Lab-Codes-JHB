% first generate MDA files using spikegadgets builtin functions

daydir=uigetdir();
olddir=cd;
cd(daydir);
recfiles=dir('**/*.rec');
cd(olddir);

myorder=cellfun(@(a) str2double(a(end-4) ),{recfiles.name});

recstring=[];
for i=1:length(recfiles)
    recpos=find(myorder==i);
    recstring=[recstring ' -rec ' fullfile(recfiles(recpos).folder,recfiles(recpos).name)];
end
% these recfiles are always in epoch folders, so you have to go back two
[~,outname]=fileparts(fileparts(recfiles(1).folder));
outdir=fileparts(recfiles(1).folder);
outstring=[' -output ' outname ' -outputdirectory ' outdir];
%eval(['! /home/jadhav/SpikeGadgets/exportmda' ' -rec ' fullfile(recfiles(myorder==1).folder,recfiles(myorder==1).name) outstring]);

% you must cd into your spikegadgets directory or add to path!!!!!!!!
eval(['! /home/jadhav/SpikeGadgets/exporttime' recstring outstring]);
eval(['! /home/jadhav/SpikeGadgets/exportmda' recstring outstring]);
eval(['! /home/jadhav/SpikeGadgets/exportLFP' recstring outstring]);])


%%


% justins code

% Run mda_util, returns list of tetrode results directories
    
%ml_process_animal(animalprefix,datadir,'tet_list',[12 13 14 17 18 19 20])

tet_list=[3, 5, 10, 17, 18, 20, 24, 28, 31, 32, 34, 36, 39, 40, 4, 9, 64];
%dayDir=uigetdir;

%resDirs = mda_util(dayDir,'tet_list',tet_list,'topDir','/media/jadhav/DATA/Jay','dataDir',fullfile(dayDir,'res'));

animals = {'XFB3'};

for a = 1:length(animals)
    animalprefix = animals{a};
    datadir = sprintf('/media/jadhav/DATA/Jay/%s/',animalprefix);
    ml_process_animal(animalprefix,datadir,'tet_list',tet_list)

end

%% now open my viewer
%{
 open terminal and run the below:
conda activate mlab'
/media/jadhav/DATA/Jay/XFB3_direct/MountainSort/XFB3_04_210827.mountain/XFB3_04_210827.nt39.mountain' );
/media/jadhav/DATA/Jay/_direct/MountainSort/__.mountain/__.nt39.mountain
qt-mountainview --raw=raw.mda --filt=filt.mda.prv --pre=pre.mda.prv --firings=firings_raw.mda --cluster_metrics=metrics_tagged.json --samplerate=30000
%}

%% this is for a single day

daydir=uigetdir();
olddir=cd;
cd(daydir);
recfiles=dir('**/*.rec');
cd(olddir);

myorder=cellfun(@(a) str2double(a(end-4) ),{recfiles.name});

recstring=[];
for i=1:length(recfiles)
    recpos=find(myorder==i);
    recstring=[recstring ' -rec ' fullfile(recfiles(recpos).folder,recfiles(recpos).name)];
end
% these recfiles are always in epoch folders, so you have to go back two
[~,outname]=fileparts(fileparts(recfiles(1).folder));
outdir=fileparts(recfiles(1).folder);
outstring=[' -output ' outname ' -outputdirectory ' outdir];
%eval(['! /home/jadhav/SpikeGadgets/exportmda' ' -rec ' fullfile(recfiles(myorder==1).folder,recfiles(myorder==1).name) outstring]);


eval(['! /home/jadhav/SpikeGadgets/exporttime' recstring outstring]);
eval(['! /home/jadhav/SpikeGadgets/exportmda' recstring outstring]);



% Run mda_util, returns list of tetrode results directories
%resDirs = mda_util(dayDirs,'tet_list',tet_list,'dataDir',dataDir);
dayDirs={daydir};
%dayDirs={uigetdir('/media/jadhav/DATA','Choose a recording day')};
tet_list=[3, 5, 10, 17, 18, 20, 24, 28, 31, 32, 34, 36, 39, 40, 4, 9, 64];

resDirs = mda_util(dayDirs,'tet_list',tet_list);
dayResDirs = cell(numel(dayDirs),1);
dayIdx = 1;
maskErrors = zeros(numel(resDirs),1);

% For each day and tet sort spikes
for k=1:numel(resDirs)
    rD = resDirs{k};
    diary([rD filesep 'ml_sorting.log'])
    fprintf('\n\n------\nBeginning analysis of %s\nDate: %s\n\nBandpass Filtering, Masking out artifacts and Whitening...\n------\n',rD,datestr(datetime('Now')));
    
    % filter mask and whiten
    % TODO: Add check to make sure there is data in the mda files, maybe up in mda_util
    [out,maskErrors(k)] = ml_filter_mask_whiten(rD,'mask_artifacts',mask_artifacts);
    % returns path to pre.mda.prv file
    
    fprintf('\n\n------\nPreprocessing of data done. Written to %s\n------\n',out)
    
    fprintf('\n------\nBeginning Sorting and curation...\n------\n')
    
    % Sort and curate
    out2 = ml_sort_on_segs(rD);
    % returns paths to firings_raw.mda, metrics_raw.json and firings_curated.mda
    %fprintf('\n\nSorting done. outputs saved at:\n    %s\n    %s\n    %s\n',out2{1},out2{2},out2{3})
    fprintf('\n\n------\nSorting done. outputs saved at:\n    %s\n    %s\n------\n',out2{1},out2{2})
    
    % Delete intermediate files
    %if ~keep_intermediates
    %    tmpDir = '/tmp/mountainlab-tmp/';
    %    disp('Removing intermediate processing mda files...')
    %    delete([tmpDir '*.mda'])
    %end
    %
    % Create matclust params file: Not Needed, only for trying to view spikes in matclust 5/8/19
    %fprintf('\n\n------\nCreating Matclust Params and Waves files\n------\n')
    %out3 = generateMatclustFromMountainSort(rD);
    %fprintf('\n\n------\nFile creation done. outputs saved at:\n    %s\n    %s\n------\n',out3{1},out3{2})
    
    if maskErrors(k)
        fprintf('\n######\nMasking error for this day. Masking Artifacts skipped. Spikes may be noisy\n######\n')
    end
    
    diary off
    
    % check if container results folders is already in dayResDirs and add if not
    if rD(end)==filesep
        rD = rD(1:end-1);
    end
    dD = fileparts(rD);
    if ~any(strcmpi(dayResDirs,dD))
        dayResDirs{dayIdx} = dD;
        dayIdx = dayIdx+1;
    end
end
fprintf('Completed automated clustering!\n')



%% here is me renaming my files
 

datadir = '/media/jadhav/DATA/Jay/XFB3/';
    
daydirs=dir(datadir);
daydirs(1:2)=[];
for i=1:length(daydirs)
    subdirs=dir(fullfile(daydirs(i).folder, daydirs(i).name));
    subdirs(1:2)=[];
    for k=1:length(subdirs)
        pat='(?<PI>[A-Z]+)_(?<anim>[A-Z]+[0-9]+)_(?<date>[^_]+)_(?<run>[^_]+)';
        oldname=regexp(subdirs(k).name,pat,'names');
        if ~isempty(oldname)
            datesep=find(oldname.date=='-');
            oldyear=oldname.date(end-1:end); oldmon=sprintf('%02d',str2double(oldname.date(1:datesep(1)-1)));
            oldday=sprintf('%02d',str2double(oldname.date(datesep(1)+1:datesep(2)-1)));
            newname=['XFB3' '_' '05_' oldyear oldmon oldday '_' oldname.run];
            movefile(fullfile(subdirs(k).folder, subdirs(k).name), fullfile(subdirs(k).folder, newname));
        end
    end
end


   
datadir = '/media/jadhav/DATA/Jay/XFB3/';

daydirs=dir(datadir);
daydirs(1:2)=[];
for i=2:length(daydirs)
    subdirs=dir(fullfile(daydirs(i).folder, daydirs(i).name));
    subdirs(1:2)=[];
    for j=2:length(subdirs)
        pat='(?(?<anim>[A-Z]+[0-9]+)_(?<sesn>[^_]+)_(?<date>[^_]+)_(?<run>[^_]+)';
        oldname=regexp(subdirs(j).name,pat,'names');
        allfiles=dir(fullfile(subdirs(j).folder, subdirs(j).name));
        allfiles(1:2)=[];
        for k=1:length(allfiles) 
            postfix=allfiles(k).name(find(allfiles(k).name=='.',1,'first'):end);
            newname=[oldname.anim '_' oldname.sesn '_' oldname.date '_' oldname.run postfix(find(postfix=='.',1,'last'):end)];
            movefile(fullfile(allfiles(k).folder, allfiles(k).name), fullfile(allfiles(k).folder, newname));
        end
    end
end
    
%%

daydirs=dir(datadir);
daydirs(1:2)=[];
for i=1:length(daydirs)
    subdirs=dir(fullfile(daydirs(i).folder, daydirs(i).name));
    subdirs(1:2)=[];
    for k=1:7
        allfiles=dir(fullfile(subdirs(k).folder, subdirs(k).name));
        allfiles(1:2)=[];
        for j=1:length(allfiles)
            oldname=allfiles(j).name;
            newname=[oldname(1:6) '4_' oldname(8:end)];
            movefile(fullfile(allfiles(j).folder, allfiles(j).name), fullfile(allfiles(j).folder, newname));
        end
    end
end

%%
%
%
%
%
%    THis converts the ml data to filter framework...
%     I will be converting filter framework into my own data
%
%
%
%
%
parentdir=uigetdir('/media/jadhav/DATA/Jay','grab a directory of MS result folders');
tetResDirs=dir(parentdir);
dayResDirs={parentdir};

%tetResDir = [resDir parsed.anim '_' parsed.day '_' parsed.date '.nt' parsedF.tet '.mountain' filesep];
%animID=
% dayResDirs are the tetrode results directories for that day
for k=1:numel(dayResDirs)
       dD = dayResDirs{k};
       fprintf('Creating spikes file for %s...\n',dD);
       [remainder,dirName] = fileparts(dD);
       if isempty(dirName)
           [~,dirName] = fileparts(remainder);
       end
       animID=dirName(1:find(dirName=='_',1,'first')-1);
       pat = '\w*_(?<day>[0-9]+)_\w*';
       parsed = regexp(dirName,pat,'names');
       dayNum = str2double(parsed.day);
       fprintf('Identified Day Number as %02i\n',dayNum);
       % db is the animal directory
       spikesFile = convert_ml_to_FF_uncurated(animID,dD,dayNum);
       fprintf('Done! Created %s\n\n',spikesFile)
end

%% another workspace...

convert_ml_to_struct(dayDir);


%% now plotting my spikes!
goodunits=units([units.noise_overlap]<.015 & [units.isolation]>.95);
goodunits=units([units.meanrate]<10);
spikeColors=parula(length(goodunits));
timeframe=[3650 3710];
it=1;
tempspikes=[];
sp=subplot(2,1,1);
for i=1:length(goodunits)
    okspikes=goodunits(i).spikes(goodunits(i).spikes(:,2)>timeframe(1) & ...
        goodunits(i).spikes(:,2)<timeframe(2),2);
    tempspikes=[tempspikes; okspikes];
    if length(okspikes)<10
        continue
    else
        spikeTicsx=repmat(okspikes',3,1);
        spikeTicsy=repmat([it-1; it; nan],1,size(spikeTicsx,2));
        plot(spikeTicsx(:),spikeTicsy(:),'color',spikeColors(i,:));
        hold on;
        it=it+1;
    end
end
box off;
sp(2)=subplot(2,1,2);
step=.2;
[y,x]=histcounts(tempspikes,timeframe(1):step:timeframe(2));
plot(x(2:end)-step,zscore(smoothdata(y,'gaussian',6)));
linkaxes(sp,'x');



