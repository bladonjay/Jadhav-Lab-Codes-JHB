%adds day number to raw day directories and files within them
clear

animal = 'CS42';
animDir = ['Z:\Projects\cs_OdorPlace_Project\RAWDataFolders\',animal,'\']; 
dayDirs = dir(animDir);
dayDirs = {dayDirs(3:end).name};


for d = 1:length(dayDirs)
    dayDir = dayDirs{d};
    daynumstr = dayDir(1:2);
    datestr = dayDir(4:11);
    dayDir = [animDir, dayDir];
    
    dataFolders = dir(dayDir);
    fol = [dataFolders.isdir];
    dataFolders = dataFolders(fol);
    dataFolders = {dataFolders(3:end).name};
    
    [commonStr] = RN_findCommonPrefix(dataFolders);
    
    
    for f = 1:length(dataFolders)
        dataFolder = [dayDir,filesep,dataFolders{f}];
        
        data = dir(dataFolder);
        data = {data(3:end).name};
        
        for t = 1:length(data)
            oldfile = [dataFolder,filesep,data{t}];
            oldstr = data{t};
            sep = strfind(oldstr,'.');
            keepsuffix = oldstr(sep(1):end);
            oldstr = oldstr(1:sep(1)-1);
            %newstr = [commonStr(1:4),'_',daynumstr,commonStr(5:end)];
            %if ~strcmp(oldstr,newstr)
            
                newfile = [dataFolder,filesep,commonStr(1:4),'_',daynumstr,'_',datestr,keepsuffix];
                movefile(oldfile,newfile)
            %end
        end
        
        oldstr = dataFolders{f};
        sep = strfind(oldstr,'.');
            keepsuffix = oldstr(sep(1):end);
        newFolder = [dayDir,filesep,commonStr(1:4),'_',daynumstr,'_',datestr,keepsuffix];
        movefile(dataFolder,newFolder)
    end
end