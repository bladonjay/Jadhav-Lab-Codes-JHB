dayfolder=uigetdir('G:\SocialData\','select a day folder');

% now go and find each rec file


% now order them and use the new trodesexport but and designate the outputs
% first we need to rename these files...
% the convention is:
% XFB3
%   02_20210830
%       XFB3_02_20210830_sleep1
%           XFB3_02_20210830_sleep1.rec, etc.
%

% change parent folder:

subdirs=dir(dayfolder);
subdirs=subdirs([subdirs.isdir]==1); subdirs(1:2)=[];

for i=1:length(subdirs)
    namesep=find(subdirs(i).name=='_');
    newname=['XFB3_03_20210827' subdirs(i).name(namesep(end):end)];
    movefile(fullfile(subdirs(i).folder,subdirs(i).name),fullfile(subdirs(i).folder,newname));
    allfiles=dir(fullfile(subdirs(i).folder,newname));
    allfiles=allfiles([allfiles.isdir]==0);
    for j=1:length(allfiles)
        % parse file
        namesep=find(allfiles(j).name=='_');
        fileend=find(allfiles(j).name=='.'); % get all the dots too
        newname=['XFB3_03_20210827' allfiles(j).name([namesep(end):fileend(1) fileend(end)+1:end])];
        movefile(fullfile(allfiles(j).folder,allfiles(j).name),fullfile(allfiles(j).folder,newname));
    end
end



%%
% concatenate in your special order:
olddir=cd;
cd(dayfolder);
% pull all behavior files
recFiles=dir('**/*.rec');
cd(olddir);
runorder=cellfun(@(a) str2double(a(end-4)), {behFiles.name});
[~,fileinds]=sort(runorder);

% cat all the rec files
recstring=[];
for i=1:length(recFiles)
    recstring=[recstring '-rec ' fullfile(recFiles(fileinds(i)).folder,recFiles(fileinds(i)).name)...
        ' '];
end

% and out two because first rec file is in an epoch directory
outstring=['-output XFB3_03_20210827 -outputdirectory ' fileparts(recFiles(1).folder)];


addpath(genpath(('C:\Users\Jadhavlab\Documents\gitRepos\TrodesFullPackage\Trodes_2-2-3_Windows64')))
eval(['! exporttime ' recstring outstring]);
eval(['! exportmda ' recstring outstring]);


rmpath(genpath(('C:\Users\Jadhavlab\Documents\gitRepos\TrodesFullPackage\Trodes_2-2-3_Windows64')))
cd(olddir)

%% if those files were not properly named... lets fix that

parentdir=uigetdir('F:\SocialData\Neural\XFB3\03_20210827');
[~,realname]=fileparts(parentdir);
myfiles=dir(parentdir); myfiles([myfiles.isdir]==1)=[];
for i=1:length(myfiles)
    namesep=find(myfiles(i).name=='.',1,'first');
    newname=[realname myfiles(i).name(namesep:end)];
    movefile(fullfile(myfiles(i).folder,myfiles(i).name),fullfile(myfiles(i).folder,newname));
end
%%

% need to change the names here

mydir=uigetdir();

subdirs=dir(mydir);

for i=3:length(subdirs)
    oldname=subdirs(i).name;
    cutspots=oldname(fnd(oldname=='_'));
    newname='JHB_XFB3_';
end 
    
    