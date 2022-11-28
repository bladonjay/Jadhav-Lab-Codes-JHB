function periods = cs_parseDIOtimes(animal, day, epoch, dionum)

topDir = cs_setPaths();
animDir = [topDir,animal,'Expt\',animal,'_direct\'];

dio = loaddatastruct(animDir, animal, 'dio', day);

dio = dio{day}{epoch}{dionum};

if dio.state(1) == 1
    on = dio.time(1:2:end);
    off = dio.time(2:2:end);
elseif dio.state(1) == 0
    on = dio.time(2:2:end);
    off = dio.time(3:2:end);
end

if length(on) > length(off)
   on = on(1:end-1);
end
periods = [on,off];