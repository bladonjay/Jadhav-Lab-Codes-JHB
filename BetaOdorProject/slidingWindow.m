function [Y,winIdx] = slidingWindow(X,winSize,winStep,varargin)
    % [Y,winIdx] = slidingWindow(X,winSize,winStep,NAME,VALUE) takes vector X
    % and returns a matrix Y where each column is a window of winSize of X and
    % each window moves step winStep along X. winIdx is a matrix when each
    % column contains the start and end indices of each window. 
    % NAME-VALUE Pairs:
    %       fullIdxOutput - if set to 0 (default) winIdx is just the start and
    %                       end indices of each window, if set to 1 then winIdx will contain
    %                       all indices in each window
    fullIdxOutput = 0;
    assignVars(varargin{:});

    N = numel(X);
    winIdx = 1:winStep:N-winSize+1;
    winIdx = [winIdx;winIdx+winSize-1];
    idxMat = arrayfun(@(x,y) (x:y)',winIdx(1,:),winIdx(2,:),'UniformOutput',false);
    idxMat = cell2mat(idxMat);
    Y = X(idxMat);
    if fullIdxOutput
        winIdx = idxMat;
    end

    if isrow(Y)
        Y = Y';
    end
