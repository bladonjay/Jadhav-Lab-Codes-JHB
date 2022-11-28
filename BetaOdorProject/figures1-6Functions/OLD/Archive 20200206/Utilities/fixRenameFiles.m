directoryname = 'E:\Data\PreTraining\CS31\';
directory = dir(directoryname);
folders = directory(3:end);
lengthofdir = length(folders);


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

    
    