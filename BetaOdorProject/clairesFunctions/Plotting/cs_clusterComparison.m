epoch = 2;
tet =31;

load('D:\OdorPlaceAssociation\CS33Expt\CS33spikes01_MS.mat')
spikes_ms = spikes;

load('D:\OdorPlaceAssociation\CS33Expt\CS33_direct\CS33spikes01.mat');

spikes_old = spikes{1}{epoch}{tet};
spikes_old = spikes_old(~cellfun(@isempty,spikes_old));

spikes_new = spikes_ms{1}{epoch}{tet};
spikes_new = spikes_new(~cellfun(@isempty,spikes_new));


meanrates_o = cell2mat(cellfun(@(x) x.meanrate, spikes_old,'UniformOutput',false));
meanrates_n = cell2mat(cellfun(@(x) x.meanrate, spikes_new,'UniformOutput',false));

% [~,i_o] = sort(meanrates_o);
% [~, i_n] = sort(meanrates_n);
i_o = [7 2 8 1 5 3 4];
i_n = [5 1 4 6 2 3];

figure,
subplot(2,1,1)
hold on
for i = 1:length(i_o)
    data = spikes_old{i_o(i)}.data(:,1);
    plot(data, repmat(i,length(data),1), '.','MarkerSize',10)
    axis([2195, 2235, 0, 8])
    title('Matclust')
end
vline([2200,2210,2220,2230],{'k','k','k','k'})
xticks([2200,2210,2220,2230])
xticklabels({5 15 25 35})
xlabel('Time (seconds)')
yticks(1:7)
ylabel('Cell Number')
subplot(2,1,2)
hold on
for i = 1:length(i_n)
    data = spikes_new{i_n(i)}.data(:,1);
    plot(data, repmat(i,length(data),1), '.','MarkerSize',10)
    axis([2195, 2235, 0, max(i_n)+1])
    title('MountainSort')
end
vline([2200,2210,2220,2230],{'k','k','k','k'})
xticks([2200,2210,2220,2230])
xticklabels({5 15 25 35})
yticks(1:6)
xlabel('Time (seconds)')
ylabel('Cell Number')
