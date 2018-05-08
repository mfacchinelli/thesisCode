function [figSizeLarge,figSizeSmall] = saveFigureSettings(saveFigure)

%...Select figure sizes based on number of outputs
if nargout == 1
    if saveFigure % set plot visibility
        set(0,'DefaultFigureVisible','off'), figSizeLarge = [1,1,1440,827];
    else
        set(0,'DefaultFigureVisible','on'), figSizeLarge = [440,378,840,630];
    end
elseif nargout == 2
    if saveFigure % set plot visibility
        set(0,'DefaultFigureVisible','off'), figSizeLarge = [1,1,1440,827]; figSizeSmall = [440,378,560,420];
    else
        set(0,'DefaultFigureVisible','on'), figSizeLarge = [440,378,840,630]; figSizeSmall = figSizeLarge;
    end
end