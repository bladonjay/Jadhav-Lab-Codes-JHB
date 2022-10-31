function cs_fixFilenames(dayDir, oldPrefix, newPrefix)

directoryname = 'E:\Data\PreTraining\CS31\';
directory = dir(directoryname);
folders = directory(3:end);
lengthofdir = length(folders);


 prefix = inputdlg({'Change common prefix from:','To:'},'Change File Prefixes',1,{common_prefix,dayName});
    if ~isempty(prefix)
        if isempty(strfind(common_prefix,prefix{1}))
            h = msgbox('The from field must contain a valid existing common prefix.');
            waitfor(h);
        else
            common_prefix = prefix{1};
            prefix = prefix{2};
            preFlg = 0;
        end
    else
        preFlg = 0;
        prefix = common_prefix;
    end
end

% Loop through all files
allFiles = dir('*.*');
allFiles = {allFiles(~[allFiles.isdir]).name}';
if isempty(allFiles)
    if cdFlg
        cd(currDir);
    end
    return;
end

for i=1:numel(allFiles),
    fn = allFiles{i};
    % Adjust Filename
    fn1 = fn;
    if ~isempty(strfind(fn1,'.1.'))
        fn1 = strrep(fn1,'.1.','.');
    end
    for n=2:9,
        if ~isempty(strfind(fn1,sprintf('.%g.',n)))
            fn1 = strrep(fn1,sprintf('.%g.',n),sprintf('_%g.',n));
        end
    end
    if ~isempty(strfind(fn1,'videoTimeStamps.cameraHWFrameCount'))
        fn1 = strrep(fn1,'videoTimeStamps.cameraHWFrameCount','cameraHWFrameCount');
    end
    if ~isempty(prefix)
        fn1 = strrep(fn1,common_prefix,prefix);
    end

    fn1 = strrep(fn1,' ','_');
    
    % Change filename if needed
    if ~strcmp(fn,fn1)
        disp(['Changing filename ' fn ' to ' fn1])
        movefile(fn,fn1)
    end
end


for i = 1:lengthofdir
    foldername = folders(i).name;
    
    sep = strfind(foldername, '2016');
    
    month1stchar = foldername(sep+4);
    month1stcharnum = str2num(month1stchar);
    if month1stcharnum == 1
        mon = foldername(8:9);
        
        day = foldername(10:end);
    else
        mon = foldername(8);
        day = foldername(9:end);
    end
    newname = [mon,'-',day];
    
    movefile(foldername, newname);
end

    
    