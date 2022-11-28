function coords = getcoord_Tmaze(directoryname, pos, index)
%this program is called by CREATETASKSTRUCT to produce the trajectory
%coodinates for a T maze
%Click locations in the following order:
%        
%        1
%        |
%        |
%        |
%        |
%   2----3----4
% 
%

fid = figure;
plot(pos(:,2),pos(:,3));
[x,y] = ginput(4);
lincoord{1} = [x([1 3 2]) y([1 3 2])];
lincoord{2} = [x([1 3 4]) y([1 3 4])];
    
numtimes = size(pos,1);
for i = 1:length(lincoord)
    coords{i} = repmat(lincoord{i},[1 1 numtimes]);
end
close(fid);