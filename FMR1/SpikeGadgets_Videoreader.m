%SpikeGadgets Video Editor
% first lets see if we can load h264 videos


[myvid,viddir]=uigetfile;
vr=VideoReader(fullfile(viddir,myvid));

myframes=read(vr,[100,300]);



image(uint8(nanmean(myframes,4)));

otherframe=read(vr,40000);

figure; image(uint8(otherframe-uint8(nanmean(myframes,4))));

%% sketchpad

set(gca,'ylim',[-.05 .2]);

hold on;
% wt-fx pval
plot([1 1 3 3],[.15 .16 .16 .15],'k');
text(2,.17,'**')
plot([2 2 3 3],[.17 .18 .18 .17],'k');
text(2.5,.145,'**')