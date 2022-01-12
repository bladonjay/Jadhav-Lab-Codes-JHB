oldDir=uigetdir;
fnames=dir(oldDir);
fnames=fnames(3:end);
for i=1:length(fnames)
    % parse old file info
    % JHB_ANIMAL_DATE-Date_DATE_run.badpostfix.postfix
    %patn='(<person>[A-Z]+)_(<rat>[A-Z]+[1-9]+)_(<date>[^_.]+)_(<epoch>[^_.]+).(<postfix>)';
    patn='log(<date>[^(]+) '
    parsed=regexp(fnames(i).name,patn,'names');
end


%% fix file extensions:

oldDir=uigetdir;
homedir=cd;
allfiles=getAllFiles(oldDir,0);

