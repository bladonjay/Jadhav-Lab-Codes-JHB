%% quick tetinfo addition to the data

% adds tetinfo field to all the sessions.
% the parent directory is just whatever directory has all of
parentDir='E:\Brandeis datasets\Claire Data';
subfolders=dir(parentDir); subfolders(1:2)=[];

% for each rat, pull metadata (name, day, import all epochs)
for i=1:length(SuperRat)
    ratName=SuperRat(i).name;
    daynum=SuperRat(i).daynum;
    ratDir=subfolders(contains({subfolders.name},ratName)).name;
    tetinfoFF=load(fullfile(parentDir,ratDir,sprintf('%stetinfo.mat',ratName)));
    tetdata=tetinfoFF.tetinfo{daynum};
    
    % initialize the struct with the first day
    firstrun=SuperRat(i).RunEpochs(1);
    tetstruct=cellfun(@(a) a, tetdata{firstrun},'UniformOutput',false);
    for k=2:length(tetstruct)
        tetstruct{k}.tetnum=k;
    end
    tetstruct(cellfun(@(a) isempty(a), tetstruct))=[]; % kill empty data (sometimes happens)
    tetinfo=tetstruct{1};
    for j=1:length(tetstruct)-1

        tetinfo=MergeStructs(tetinfo,tetstruct{j+1},'addall');
    end
    % i think i can get away with just using the working epochs
    
            
    SuperRat(i).tetinfo=tetinfo;
end

    