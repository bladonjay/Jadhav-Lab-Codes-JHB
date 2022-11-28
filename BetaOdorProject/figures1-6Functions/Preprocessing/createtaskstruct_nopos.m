function createtaskstruct_nopos(directoryname,fileprefix,index)

%Creates Task strucutre without needing pos information. Just creates
%structure with "run" type for indicated epochs. 
%Index should be Nx2 matrix where each row is an epoch, first column is day
%and second column is epoch in that day.

task = [];
currday = (index(1,1));

for i = 1:size(index,1)    
    day = index(i,1);
    epoch = index(i,2); 
    newday = 0;
    
    if (day ~= currday) 
        newday = 1;
        daystr = getTwoDigitNumber(currday);
        eval(['save ',directoryname,fileprefix,'task',daystr,' ','task']);
        task = [];
        currday = day;
        
    end
    
   
    task{day}{epoch}.type = 'run';
       
    if (i == size(index,1))
        daystr = getTwoDigitNumber(currday);
        eval(['save ',directoryname,fileprefix,'task',daystr,' ','task']);
    end
    
    
end
