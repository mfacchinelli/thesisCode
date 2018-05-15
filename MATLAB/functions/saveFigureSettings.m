function [figSizeLarge,figSizeMedium,figSizeSmall] = saveFigureSettings(saveFigure)

%...Select figure sizes based on number of outputs
if nargout == 1
    if saveFigure % set plot visibility
        set(0,'DefaultFigureVisible','off'), figSizeLarge = [1,1,1440,827];
    else
        set(0,'DefaultFigureVisible','on'), figSizeLarge = [440,378,840,630];
    end
elseif nargout == 2
    if saveFigure % set plot visibility
        set(0,'DefaultFigureVisible','off');
        figSizeLarge = [1,1,1440,827];
        figSizeMedium = [440,378,840,630];
    else
        set(0,'DefaultFigureVisible','on');
        figSizeLarge = [440,378,840,630];
        figSizeMedium = figSizeLarge;
    end
elseif nargout == 3
    if saveFigure % set plot visibility
        set(0,'DefaultFigureVisible','off');
        figSizeLarge = [1,1,1440,827];
        figSizeMedium = [440,378,840,630];
        figSizeSmall = [440,378,560,420];
    else
        set(0,'DefaultFigureVisible','on');
        figSizeLarge = [440,378,840,630];
        figSizeMedium = figSizeLarge;
        figSizeSmall = figSizeLarge;
    end
end