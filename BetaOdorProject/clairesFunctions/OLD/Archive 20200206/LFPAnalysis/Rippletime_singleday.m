clc;
close all;
clear all;
%%
addpath('C:\Users\wbtang\Desktop\HMM\SleepCode')
addpath(genpath('C:\Users\wbtang\Desktop\HMM\Src_Matlab\'));
%%
animaldir = ('D:\SingledayExp\JS21_direct\');
animalprefix = ('JS21');
day = 1;
% eps = 2:2:16; %Run epoch
eps = 1:17; %all epoch

minstd = 3;
minnum = 1;
velfilter = 4; %velthresh
nn = 0;
epochs = zeros(length(eps),2);
for ep = eps
   nn = nn+1;
   epochs(nn,1) = day;
   epochs(nn,2) = ep;
end

%find riptets
tetinfo = loaddatastruct(animaldir, animalprefix, 'tetinfo');
%get ripple times
rippletime = DFTFsj_getriptimes(animaldir,animalprefix, epochs,'tetfilter', '(isequal($descrip, ''riptet''))','minthresh',minstd);
%% find periods of immobility
if ~isempty(velfilter)
    pos = loaddatastruct(animaldir, animalprefix, 'pos', day); % get animal [time, positions, velocity]
end
for ep = eps
    riplist = vec2list(rippletime{day}{ep}.nripples >= minnum,rippletime{day}{ep}.time);
    if ~isempty(velfilter)
        if size(pos{day}{ep}.data,2) > 5  % get velocity
            velocity = pos{day}{ep}.data(:,9);% smoothed velocity
        else
            velocity = pos{day}{ep}.data(:,5);
        end
        postimes = pos{day}{ep}.data(:,1);  % get time stamps
        immobile = vec2list(velocity < velfilter,postimes); % generate [start end] list of immobile epochs
    end
    times_filteeg = rippletime{day}{ep}.time;
    duration = (times_filteeg(end) - times_filteeg(1)).*1000; % ms
    [imtime,im_vec] = wb_list2vec(immobile,times_filteeg);
    [riptime,ripvec] = wb_list2vec(riplist,times_filteeg);
    nonripvec = ~ripvec;
    if ~isempty(velfilter)
        ripvec = ripvec & im_vec;   
        riplist = vec2list(ripvec,times_filteeg);
        nonripvec = ~ripvec & ~im_vec;
    end    
    nriplist = vec2list(nonripvec,times_filteeg);
    if ~isempty(riplist)
        rip.timevec = times_filteeg;
        rip.starttime = riplist(:,1);
        rip.endtime = riplist(:,2);
        rip.ripvec = ripvec;
        rip.nonripvec = nonripvec;
        rip.nstarttime =  nriplist(:,1);
        rip.nendtime =  nriplist(:,2);
        ripduration = round(sum(riplist(:,2)-riplist(:,1)));
        rip.total_duration = ripduration;
        
         im.timevec = times_filteeg;
         im.starttime = immobile(:,1);
         im.endtime = immobile(:,2);
         im.total_duration =  round(sum(immobile(:,2)-immobile(:,1)));
    else
        rip.timevec = [];
        rip.starttime =[];
        rip.endtime = [];
        rip.ripvec = [];
        rip.total_duration = 0 ;
    end
    rip.minthresh = minstd;
    rip.minnumber = minnum;
    rip.velthresh = velfilter;
    ripple{day}{ep} = rip;
    
    im.velthresh = velfilter;
    immobility{day}{ep} = im;
    disp(sprintf('d %d e %d: %d s of ripples',day,ep,ripduration))
    disp(sprintf('d %d e %d: %d s of immobility',day,ep,round(sum(immobile(:,2)-immobile(:,1)))))
    clear rip;clear im;
end
save(sprintf('%s/%srippletime%02d.mat', animaldir, animalprefix, day), 'ripple');
%save(sprintf('%s/%simmobiletime%02d.mat', animaldir, animalprefix, day), 'immobility');

clear ripple;clear immobility


