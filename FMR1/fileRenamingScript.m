% trying to change 
% data names need to change in :F:\SocialData\Neural\XFB3\XFB3_01-05
% folders
changeDir=uigetdir;
fnames=dir(changeDir);
fnames=fnames(3:end); % remove the parent folders

for i=1:length(fnames)
    % parse old file info
    % JHB_ANIMAL_DATE-Date_DATE_run.badpostfix.postfix
    patn='(?<rat>[A-Z]{2,3}[0-9]+)_(?<runnum>\d{2}+)_(?<date>\d{8}+)_(?<run>[a-zA-Z0-9]+)';
    %patn='log(<date>[^(]+) ';
    parsed=regexp(fnames(i).name,patn,'names');
     if ~isempty(parsed)
         infiles=dir(fullfile(fnames(i).folder,fnames(i).name));
         for j=1:length(infiles)
             patn='(?<person>[A-Z]{3})_(?<rat>[A-Z]{2,3}[0-9]+)_(?<date>[^.]{1,15).(?<extension>\w+)';
             fileparsed=regexp(infiles(j).name,patn,'names');

             if ~isempty(fileparsed) & (isfield(fileparsed,'person'))
                newname=[fnames(i).name, '.', fileparsed.extension];
                movefile(fullfile(fileparsed(i).folder,fnames(i).name),fullfile(fnames(i).folder,newname));
     end
end


%% fix file extensions:

oldDir=uigetdir;
homedir=cd;
allfiles=getAllFiles(oldDir,0);


%% go throguh all files with too many dots in the name:


mypath=uigetdir('G:\','find a directory with files and folders you want to fix');
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


%% running ffmpeg

% first download ffmpeg from a stable source: just google ffmpeg and follow
% the directions.  Drag the folder with the exe file into documents or
% somewhere

% then download the matlab package 'ffmpeg' from the add-on explorer
% then run:
ffmpegsetup

% then to iterate through files:

mypath=uigetdir('G:\','find a directory with files and folders you want to fix');
allcontent=genpath(mypath); % generate all folders in tree

% now parse them all
pathsep=find(allcontent==';');

% turn this into a struct
folderstruct=struct('name',allcontent(1:pathsep(1)-1));
for i=2:length(pathsep)
    folderstruct(i).name=allcontent(pathsep(i-1)+1:pathsep(i)-1);
end
folderstruct(i+1).name=allcontent(pathsep(i)+1:end);

% now go folder to folder converting all the files from h264 to mp4
for i=1:length(folderstruct)
    allfiles=dir(folderstruct(i).name);
    realfiles=allfiles([allfiles.isdir]==0);
    if ~isempty(realfiles)
        for j=1:length(realfiles)
            [~,myname,realext]=fileparts(realfiles(j).name);
            if strcmpi(realext,'.h264') && ~exist(fullfile(realfiles(j).folder,[myname '.mp4']),'file')
                
                % for some reason i cant have it input strings for
                % fileconventions, i think maybe matlab pulls the '' out
                ffmpegexec(['-i ',fullfile(realfiles(j).folder,[myname,realext]),...
                    ' -c copy ' fullfile(realfiles(j).folder,[myname,'.mp4'])])
                %ffmpegtranscode(fullfile(realfiles(j).folder,[myname,realext]),fullfile(realfiles(j).folder,[myname,'.mp4']))
                %eval(['C:\Users\Jadhavlab\Documents\ffmpeg-win64-20221212\bin\ffmpeg -i ''', fullfile(realfiles(j).folder,[myname,realext]),...
                %    ''' -c copy ''', fullfile(realfiles(j).folder,[myname,'.mp4']),'''']);
                fprintf('\n Vidfile %s created \n',[myname,'.mp4']);
            end
            
        end
    end
end

%%


%%



