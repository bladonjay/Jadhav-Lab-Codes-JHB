%change depending on which animal's data is being extraced
animal = 'CS44';
[topDir] = cs_setPaths(); 
%animDir = [topDir,animal,'Expt\',animal,'\'];
%animDir = 'F:\Data\OdorPlaceAssociation\CS41\';
animDir = [topDir,'RAW\',animal,'\'];

cd(animDir);
dayfolders = dir();
dayNames= {dayfolders(3:end).name};

numDays= length(dayNames);
    
    for day = [8]
    
        dayfolder = dayNames{1,day};
        tempdir = strcat([animDir, dayfolder, '\']);
        cd(tempdir);
        sessiondate = dayfolder(3:end);
        filemask = [animal, sessiondate];

         %extractTimeBinaryFile(filemask);
%         extractSpikeBinaryFiles(filemask);
         extractLFPBinaryFiles(filemask);
%        extractMDABinaryFiles(filemask)
         %extractDioBinaryFiles(filemask);
        cd(animDir);
    end