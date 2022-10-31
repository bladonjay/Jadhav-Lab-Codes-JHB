%change depending on which animal's data is being extraced
animal = 'CS41';

%change this to set "topDir" to be your top Data folder, containing the
%ThermocoupleTest folder
[topDir] = cs_setPaths(); 

animDir = [topDir,animal,'Expt\',animal,'\'];

cd(animDir);
dayfolders = dir();
dayNames= {dayfolders(3:end).name};

numDays= length(dayNames);
    
    for day = 1:numDays
    
        dayfolder = dayNames{1,day};
        tempdir = strcat([animDir, dayfolder, '\']);
        cd(tempdir);
        sessiondate = dayfolder(3:end);
        filemask = [animal, sessiondate];

        extractTimeBinaryFile(filemask);
        %extractSpikeBinaryFiles(filemask);
        %createAllMatclustFiles
        extractLFPBinaryFiles(filemask);
        %extractMDABinaryFiles(filemask)
        extractDioBinaryFiles(filemask);

        cd(topDir);
    end