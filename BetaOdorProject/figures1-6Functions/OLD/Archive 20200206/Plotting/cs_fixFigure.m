function cs_fixFigure(figfile, figtype, varargin)

region = '';

for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'region'
            region = varargin{option+1};
        case 'win'
            win = varargin{option+1};
    end
end

figdir = 'D:\Figures\';
%figdir = 'F:\Figures\';


switch figtype
    case 'specgram' 
        
        openfig([figdir, 'Specgrams\',figfile]);
        figtitle = erase(figfile,'.fig'); 
        
%         if strcmp(region, 'OB')
%             coloraxis = [-0.3 2.3];
%         else
%             coloraxis = [-0.28 0.37];
%         end

        coloraxis = caxis;
        y = ylim;
        %caxis(coloraxis)
        colorbar('YTick', [coloraxis(1) 0 coloraxis(2)]); 
        
        %xlabel('Time (seconds)')
        %ylabel('Frequency (Hz)')
        xticks([]);
        yticks([]);
        
        title('');
        set(gca,'fontsize',30);
        
        line = findobj(gca,'Type','Line');
        delete(line)
        plot([0 0],[y(1) y(2)],'k--','LineWidth',6);
        
        set(gcf, 'Position', [2000 50 900 600]);
        %set(gcf, 'Position', [50 50 1200 900]);
        %set(gcf, 'Position', [1500 50 1200 250]);
        
        %newfigfile = [figdir,'NicePPTFigures\NRSA_OBspec'];
        newfigfile = [figdir,'NicePPTFigures\',figtitle];
        saveas(gcf,newfigfile,'fig');
        print('-dpdf', newfigfile);
        print('-djpeg', newfigfile);
    
    case 'cohgram'
        openfig([figdir, 'Cohgrams\',figfile]);
        figtitle = erase(figfile,'.fig'); 
        
        coloraxis = caxis;
        y = ylim;
        
        %coloraxis = [0.51 0.85];
        %y = coloraxis;
        caxis(coloraxis)
        colorbar('YTick', [coloraxis(1) coloraxis(2)]); 
        colormap(hot)
        xlabel('Time (seconds)')
        ylabel('Frequency (Hz)')
        
        title('');
        set(gca,'fontsize',30);
        
        line = findobj(gca,'Type','Line');
        delete(line)
        plot([0 0],[y(1) y(2)],'k--','LineWidth',6);
        

        set(gcf, 'Position', [50 50 1200 900]);
        set(gca, 'XLim', [-win(1) win(2)]) %change depending on window
        
        newfigfile = [figdir,'NicePPTFigures\',figtitle];
        saveas(gcf,newfigfile,'fig');
        print('-dpdf', newfigfile);
        print('-djpeg',newfigfile);
        
    case 'phaselock'
        openfig([figdir, 'PhaseLocking\',figfile]);
        figtitle = erase(figfile,'.fig'); 
        
         h = gcf; %current figure handle
         axesObjs = get(h, 'Children');  %axes handles
        dataObjs = get(axesObjs, 'Children'); %handles to low-level graphics objects in axes
        objTypes = get(dataObjs, 'Type');  %type of low-level graphics object
        xdata = get(dataObjs, 'XData');  %data from low-level grahics objects
        ydata = get(dataObjs, 'YData');
        
        disp('Done')
        
end

