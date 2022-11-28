function cs_errorbar(x,y,err,varargin)

%x = xvlaue to plot
%y = mean of the data
%err = error bars


color = 'red';
lineWidth = 2;


for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'color'
            color = varargin{option+1};
        case 'lineWidth'
            lineWidth = varargin{option+1};
    end
end


gcf;
plot(x,y,'o','MarkerEdgeColor',color,'MarkerFaceColor',color)
hold on
plot([x x],[y+err y-err],'-','LineWidth',lineWidth,'Color',color)

