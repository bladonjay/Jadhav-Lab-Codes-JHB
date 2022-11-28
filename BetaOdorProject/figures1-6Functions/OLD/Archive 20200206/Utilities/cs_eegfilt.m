topDir = cs_setPaths();
filtPath= 'D:\DataAnalysis\usrlocal\filtering\';


animID= 'CS31'; 
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory
days = 3;

for sessionNum = 1:3
    disp(['doing CS31 day ',num2str(sessionNum), ' resp']);
mcz_respdayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'respfilter.mat']);
end

animID= 'CS33'; 
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory

for sessionNum = 1:4
    disp(['doing CS33 day ',num2str(sessionNum), ' resp']);
mcz_respdayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'respfilter.mat']);
end


animID= 'CS34'; 
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory

for sessionNum = 4
    disp(['doing CS34 day ',num2str(sessionNum), ' resp']);
mcz_respdayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'respfilter.mat']);
end


animID= 'CS35'; 
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory

for sessionNum = 1:3
    disp(['doing CS35 day ',num2str(sessionNum), ' resp']);
mcz_respdayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'respfilter.mat']);
end

animID= 'CS39'; 
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory

for sessionNum = 1:7
    disp(['doing CS39 day ',num2str(sessionNum), ' resp']);
mcz_respdayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'respfilter.mat']);
end


animID= 'CS41'; 
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory

for sessionNum = [1 2 3 5 6 7 8 9 10]
    disp(['doing CS41 day ',num2str(sessionNum), ' resp']);
mcz_respdayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'respfilter.mat']);
end


animID= 'CS42'; 
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory

for sessionNum = 1:9
    disp(['doing CS42 day ',num2str(sessionNum), ' resp']);
mcz_respdayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'respfilter.mat']);
end


animID= 'CS44'; 
dataDir= [topDir, animID, 'Expt\',animID, '_direct\'];    % data directory

for sessionNum = 1:7
    disp(['doing CS44 day ',num2str(sessionNum), ' resp']);
mcz_respdayprocess(dataDir, animID, sessionNum, 'f', [filtPath 'respfilter.mat']);
end