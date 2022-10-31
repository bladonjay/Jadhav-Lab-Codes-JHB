%cs_copyMSfiles_essentialOnly
animal = 'CS44';
srcDir = ['G:\Data\OdorPlaceAssociation\CS44Expt\CS44_direct\MountainSort'];
%rawDir = ['G:\Data\RAW\CS44'];
destDir = ['D:\OdorPlaceAssociation\',animal,'Expt\',animal,'_direct\MountainSort']; %Destination MountainSort folder

dayFolders = dir(srcDir);
%dayFolders = dir(rawDir);
dayFolders = {dayFolders(3:end).name};

for d = 1:length(dayFolders)
    dayFolder = dayFolders{d};
    tetFolders = dir([srcDir filesep dayFolder filesep '*.mountain']);
    tetFolders = {tetFolders.name};
    disp(['Doing day ',num2str(d), ' of ',num2str(length(dayFolders))])
%     mdaFolder = dir([rawDir filesep dayFolder filesep '*.mda']);
%     mdaFolder = mdaFolder.name;
    destFolder = [destDir filesep animal,'_',dayFolder,'.mountain'];
%     if ~isdir([destDir filesep animal,'_',dayFolder,'.mountain'])
%         continue
%     end
    
%     cd([rawDir filesep dayFolder filesep mdaFolder])
%     timestampsfile = dir('*timestamps.mda');
%     timestampsfile = timestampsfile.name;
    
    
    for t = 1:length(tetFolders)
        tetFolder = tetFolders{t};
        cd([srcDir filesep dayFolder filesep tetFolder])
        
        destFolder = [destDir filesep dayFolder filesep tetFolder];
        if ~isdir(destFolder)
            mkdir(destFolder)
        end
        out1 = copyfile('firings.curated.mda',destFolder);
        %out2 = copyfile('metrics_tagged.json', destFolder);
        if out1 ==0 %|| out2 ==0
            error('Copying failed, check directories');
        end
    end
%     out = copyfile(timestampsfile,destFolder);
%     if out ==0
%         error('Copying failed, check directories');
%     end
end