%  MATLAB Function < subplotTitle >
%
%  Author:      Michele Facchinelli
%  Date:        10th April 2017
%  Purpose:     give overall title to subplots; reference: 
%  Input:
%   - title:    string with title name
%  Output:
%   - N/A

function subplotTitle(title)

%...Add title to plot
axes('Position',[0,0.95,1,0.05])
set(gca,'Color','None','XColor','None','YColor','None')
text(0.5,0.25,title,'FontSize',14','FontWeight','Bold','HorizontalAlignment','Center','VerticalAlignment','Bottom')