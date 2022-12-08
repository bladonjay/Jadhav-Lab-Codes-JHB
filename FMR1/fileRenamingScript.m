oldDir=uigetdir;
fnames=dir(oldDir);
fnames=fnames(3:end);
for i=1:length(fnames)
    % parse old file info
    % JHB_ANIMAL_DATE-Date_DATE_run.badpostfix.postfix
    %patn='(<person>[A-Z]+)_(<rat>[A-Z]+[1-9]+)_(<date>[^_.]+)_(<epoch>[^_.]+).(<postfix>)';
    patn='log(<date>[^(]+) ';
    parsed=regexp(fnames(i).name,patn,'names');
end


%% fix file extensions:

oldDir=uigetdir;
homedir=cd;
allfiles=getAllFiles(oldDir,0);


%% go throguh all files with too many dots in the name:


mypath=uigetdir('F:\','find a directory with files and folders you want to fix');
allcontent=genpath(mypath); % generate all filenames

% now parse them all
pathsep=find(allcontent==';');

filestruct=struct('name',allcontent(1:pathsep(1)-1));
for i=2:length(pathsep)
    filestruct(i).name=allcontent(pathsep(i-1)+1:pathsep(i)-1);
end
filestruct(i+1).name=allcontent(pathsep(i)+1:end);

% so this generates all the folders where you could have misnamed files

for i=1:length(filestruct)
    allfiles=dir(filestruct(i).name);
    realfiles=allfiles([allfiles.isdir]==0);
    if ~isempty(realfiles)
        for j=1:length(realfiles)
            fixExtensions(fullfile(realfiles(j).folder,realfiles(j).name));
        end
    end
end
     



%% moving files from one folder to another for hosting on figshare


% move each rat folder from hugedata to bigdata and omit lfp data

sourceDir='E:\Brandeis datasets\OdorPlaceAssociation';








