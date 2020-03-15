%MakePlacePlots




%% run a quick place plot script this script makes rather crappy place fields
for rt=1:length(SuperRat)
    tracking=SuperRat(rt).tracking;
    Runs=SuperRat(rt).RunEpochs;
    if ~isempty(tracking) && ~isempty(Runs) && ~isempty(SuperRat(rt).units)
        for i=1:length(SuperRat(rt).units)
            % for each unit, grab the spikes and make a place plot
            figure;
            spikedata=SuperRat(rt).units(i).ts;
            % for each day, do a top stick and dot plot, and a bottom place
            % plot
            for k=1:length(Runs)
                thistrack=tracking.data(tracking.data(:,6)==Runs(k),:);
                thesespikes=spikedata(spikedata>thistrack(1,1) & spikedata<thistrack(end,1));
                subplot(2,length(Runs),k);
                plot(thistrack(:,2),thistrack(:,3),'k');
                hold on;
                spikex=interp1(thistrack(:,1),thistrack(:,2),thesespikes);
                spikey=interp1(thistrack(:,1),thistrack(:,3),thesespikes);
                plot(spikex,spikey,'r*');
                % not make place plot here
                subplot(2,length(Runs),k+length(Runs));
                session=struct('edit_coords',tracking.data(tracking.data(:,6)==Runs(k),[1:3]));
                
                % this is the really basic open field place plotting script
                [ratemap,~,finalcolormap]=cell_SmoothPlacePlot(session,SuperRat(rt).units(i),'Factor',4,'suppress',1);
                image(finalcolormap); set(gca,'YDir','normal');
                title(sprintf('Max: %.2f',max(linearize(ratemap))));
                
                % you can use another one as follows;
                
            end
            title([SuperRat(rt).units(i).type ' ' SuperRat(1).units(i).area ' ' num2str(i)]);
        end
        
    end
end


%% this does it much better




