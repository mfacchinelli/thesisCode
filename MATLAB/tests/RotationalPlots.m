if loadRotational
    F = figure('rend','painters','pos',figSizeSmall);
    for i = 1:4
        subplot(2,2,i)
        plot(interpolatedTime,rotationalPropagatedResults(:,i)-commandedQuaternions(1,i+1),'LineWidth',1.25)
%         plot(interpolatedTime,rotationalEstimatedResults(:,i)-rotationalPropagatedResults(:,i),'LineWidth',1.25)
        xlabel(timeLabel)
        ylabel(rotationLabelsDifference{i})
        set(gca,'FontSize',15)
        grid on
    end
    if saveFigures, saveas(F,'../../Report/figures/aero_rot_diff','epsc'), end
    
    F = figure('rend','painters','pos',figSizeSmall);
    for i = 1:3
        subplot(3,1,i)
        plot(interpolatedTime,rotationalPropagatedResults(:,i+5),'LineWidth',1.25)
        xlabel(timeLabel)
        ylabel(rotationLabels{i+5})
        set(gca,'FontSize',15)
        grid on
    end
    if saveFigures, saveas(F,'../../Report/figures/aero_rot_vel','epsc'), end
end

saveRot = false;
if loadFilter
    plotFreq = 1;
    if saveRot, plotFreq = 10; end
    
    %...Plot Cartesian translational motion
    F = figure('rend','painters','pos',figSizeLarge);
    for i = 1:6
        subplot(2,3,i)
        hold on
        if saveIMAN
            plot(interpolatedTime(1:plotFreq:end),abs(filterStateEstimatedResults(1:plotFreq:end,i)-...
                CartesianPropagatedResults(1:plotFreq:end,i)),'LineWidth',1.25)
            plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--')
        else
            plot(interpolatedTime(1:plotFreq:end),filterStateEstimatedResults(1:plotFreq:end,i)-...
                CartesianPropagatedResults(1:plotFreq:end,i),'LineWidth',1.25)
            plot(interpolatedTime(3:plotFreq:end),sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--')
            plot(interpolatedTime(3:plotFreq:end),-sqrt(filterCovarianceEstimatedResults(3:plotFreq:end,i)),...
                'LineWidth',1.25,'LineStyle','--','Color',[0.85,0.325,0.098])
        end
        hold off
        xlabel(timeLabel)
        if saveIMAN
            ylabel(CartesianLabelsAbsoluteDifference{i})
        else
            ylabel(CartesianLabelsDifference{i})
        end
        if saveIMAN
            set(gca,'FontSize',15,'YScale','log')
        else
            set(gca,'FontSize',15)
        end
        grid on
    end
    subplotLegend({'Difference','STD'})
    if ( saveFigures && saveIMAN ), saveas(F,'../../Report/figures/aero_filt_cart','epsc'),
    elseif ( saveFigures && saveRot ), saveas(F,'../../Report/figures/aero_withrot_cart','epsc'), end
end

if loadDependent && loadRotational
    styles = {'-','--','-.'};
    aeroAngles = dependentVariables(:,7:9);
    
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    for i = 1:3
        plot(interpolatedTime,aeroAngles(:,i),'LineWidth',1.25,'LineStyle',styles{i})
    end
    hold off
    xlabel(timeLabel)
    ylabel('Aerodynamic Angle [deg]')
    set(gca,'FontSize',15)
    grid on
    legend('Attack','Side-slip','Bank','Location','Best')
    if saveFigures, saveas(F,'../../Report/figures/aero_angles','epsc'), end
end

if loadRotational
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    plot(onboardTime,controlTorques(:,3),'LineWidth',1.25)
    plot(onboardTime,controlTorques(:,1),'LineWidth',1.25,'LineStyle','--')
    hold off
    xlabel(timeLabel)
    ylabel('Torque [N m]')
    set(gca,'FontSize',15)
    grid on
    legend('z_B','x_B','Location','Best')
    if saveFigures, saveas(F,'../../Report/figures/aero_contr_xz','epsc'), end
    
    F = figure('rend','painters','pos',figSizeSmall);
    subplot(2,1,1)
    plot(onboardTime,controlTorques(:,2),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel('Torque [N m]')
    set(gca,'FontSize',15)
    grid on
    subplot(2,1,2)
    plot(onboardTime,controlTorques(:,2),'LineWidth',1.25)
    xlabel(timeLabel)
    ylabel('Torque [N m]')
    xlim([0.7,0.77])
    set(gca,'FontSize',15)
    grid on
    if saveFigures, saveas(F,'../../Report/figures/aero_contr_y','epsc'), end
end