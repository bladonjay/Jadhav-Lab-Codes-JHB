function out = slidingAvg(x,win,step,dim)
    % takes a vector or 2D matrix and averages across it with window size win and steps of
    % step. Win should be odd, if it is even, then the window will be set to one
    % less than win. Default dimension (dim) to slide across is 1.
    % Default step is win, in order to not overlap windows. 
    
    % one helpful addition woudl be to take care of the two edges
    
    % from Roshan Nanu
   
    if ~exist('dim','var')
        dim = 1;
    end
    if ~exist('step','var')
        step = win;
    end
    n = size(x,dim);
    if dim == 2
        x = x';
    end
    if mod(win,2)~=0
        win = win-1;
    end
    
    nWin = floor((n-win)/step);
    out = zeros(nWin,size(x,2));

    % the step stops one window from the end, and will progressively grow
    % from start until max window
    for i=1:nWin
        out(i,:) = mean(x((i-1)*step+1:(i-1)*step+win+1,:),1);
    end
    if dim == 2
        out = out';
    end
end

