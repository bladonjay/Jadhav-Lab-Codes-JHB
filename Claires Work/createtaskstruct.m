function createtaskstruct(directoryname,fileprefix,index, coordprogram, varargin)
%
%createtaskstruct(directoryname,fileprefix,index, coordprogram,options)
%
% eg.
% createtaskstruct('/data25/sjadhav/RippleInterruption/RCb_direct/','RCb',[1 2; 1 4],'getcoord_wtrack');
%This program creates the 'task' structure for each day's data. This
%structure contains information about the exposure number for each epoch
%and the linearization coordinates for each run session. If a task file
%already exists for one of the days in index, then that file is
%loaded and written to to prevent unwanted deletion of data.
%
%      --INPUTS--
%
%    directoryname - example 'data99/user/animaldatafolder/', a folder 
%                    containing processed matlab data for the animal
%       fileprefix - animal specific prefix for each datafile
%    index - each row describes a run epoch - [day epoch]
%     coordprogram - a string with the name of the function that will
%                    provide the track segment location info. This
%                    program should take three inputs - 1) directoryname
%                    2) pos.data for an epoch, and 3) the epoch index [day epoch]
%                    The output of the program is a cell
%                    array, where each cell describes one trajectory on
%                    the track (do not count foreward and backward
%                    movement as 2 trajectories).  Inside each cell is
%                    a three dimentional M-by-2-by-N matrix.  The
%                    matrix gives the x and y coodinates for the M
%                    segment endpoints, across N timeframes.  N should
%                    equal the length of pos.data.  The
%                    segment end points should start with the position
%                    that the user wants to be 0 on the linear scale and 
%                    progress to the last segment.  All endpoints that
%                    are shared by two segments should only exist
%                    once per trajectory.  The function is called like
%                    this: coords = coordprogram(direcoryname,pos{day}{epoch}.data,[day epoch])
%    ---OPTIONS----
%
%        overwrite - 1 to overwrite old task file (default), 0 to keep the same file
%                    and just change the epochs defined in index
%   lowercasethree - variable prefix. Default is ''.
%    combineepochs - value of 1 combines all the epochs for the day when
%                    calling coordprogram. Default 0;  
%   
%


%set default variables
lowercasethree = '';
overwrite = 0;
combineepochs = 0;

%set variable options
for option = 1:2:length(varargin)-1
    
    switch varargin{option}
        case 'lowercasethree'
            lowercasethree = varargin{option+1};
        case 'overwrite'
            overwrite = varargin{option+1};   
        case 'combineepochs'
            combineepochs = varargin{option+1}; 
        otherwise
            error(['Option ''', varargin{option}, ''' not defined']);
    end
end

index = sortrows(index,1);

eval(['coordProgramHandle = @',coordprogram,';']);
task = [];
currday = 0;
%%%
for i = 1:size(index,1)    
    day = index(i,1);
    epoch = index(i,2); 
    newday = 0;
    if (day ~= currday)
        newday = 1;
        if (currday ~= 0)
            %we have moved on to a new day, so save the previous day's data
            dsz = '';
            if (currday < 10)
                dsz = '0';
            end
            eval([lowercasethree,'task = task;']);
            eval(['save ',directoryname,fileprefix,'task',dsz,num2str(currday),' ',lowercasethree,'task']);
        end
        currday = day;
        dsz = '';
        if (day < 10)
            dsz = '0';
        end
        %load the new day's data
        eval(['load ',directoryname,fileprefix,'pos', dsz, num2str(day), '.mat']);
        eval(['pos = ',lowercasethree,'pos;'])
        
        try
            if (overwrite == 0)
               %load the day's task data if it exists
               eval(['load ',directoryname,fileprefix,'task', dsz, num2str(day),'.mat']);
               eval(['task = ',lowercasethree,'task;'])
            else
               task = [];
            end
        catch
            task = [];
        end
    end


    if (combineepochs & newday) %we are using the same track coordinates for all epochs of the day
        posdata = [];
        for dayepoch = index((index(:,1) == day),2)'
            posdata = [posdata; pos{day}{dayepoch}.data];         
        end
        tmplinearcoord = feval(coordProgramHandle,directoryname,posdata);
        for dayepoch = index((index(:,1) == day),2)'
            S = size(pos{day}{dayepoch}.data,1);
            task{day}{dayepoch}.type = 'run';
            %all time points are the same, so just grab the number of
            %points for the epoch
            for reshape = 1:length(tmplinearcoord)
                epochlinearcoord{reshape} = tmplinearcoord{reshape}(:,:,1:S);
            end
            task{day}{dayepoch}.linearcoord = epochlinearcoord;
        end
    end
            
    if ~(combineepochs)
        %fill the structure with data
        task{day}{epoch}.type = 'run';
        if ~isempty(coordprogram)
            %run the user-defined coordinate program
            if isstruct(pos); 
                temppos = pos; clear pos;
                pos = {}; pos{day} = {};
                pos{day}{epoch} = temppos;
                clear temppos;
            end;
            task{day}{epoch}.linearcoord = feval(coordProgramHandle,directoryname,pos{day}{epoch}.data,[day epoch]);
        end
    end
end
%save the final day's task file
dsz = '';
if (currday < 10)
    dsz = '0';
end
eval([lowercasethree,'task = task;']);
eval(['save ',directoryname,fileprefix,'task',dsz,num2str(currday),' ',lowercasethree,'task']);
