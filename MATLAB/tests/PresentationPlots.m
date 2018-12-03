%% Plot Orbit

if loadEstimated
    %...Load results
    filename = ['/Users/Michele/GitHub/tudat/tudatBundle/tudatApplications/Thesis/',outputFolder,'/cartesianPropagated.dat'];
    fileID = fopen(filename,'r');
    presentationCartesianResults = textscan(fileID,repmat('%f',[1,7]),'Delimiter',',','CollectOutput',true);
    presentationTime = ( presentationCartesianResults{1}(:,1) - initialTime ) / timeConversion;
    presentationCartesianResults = presentationCartesianResults{1}(:,2:end);
    presentationCartesianResults(:,1:3) = presentationCartesianResults(:,1:3) / 1e3;
    fclose(fileID);
    
    %...Compute apo- and periapses
    presentationAltitude = rssq( presentationCartesianResults(:,1:3), 2 ) - marsRadius/1e3;
    [presentationApoapses,locPresentationApoapses] = findpeaks(presentationAltitude,'MinPeakProminence',100);
    [presentationPeriapses,locPresentationPeriapses] = findpeaks(-actualAltitude,'MinPeakProminence',25);
    presentationPeriapses = -presentationPeriapses;
    
    %...Plot trajectory
    rotatedTrajectory = (roty(-45) * presentationCartesianResults(1:end,1:3)')';
    reducedRotatedTrajectory = rotatedTrajectory(1:100:end,:);
    
    F = figure('rend','painters','pos',figSizeWideHAR);
    hold on
    plot3(reducedRotatedTrajectory(:,1),reducedRotatedTrajectory(:,2),reducedRotatedTrajectory(:,3),'LineWidth',1.5)
    [x,y,z] = sphere; surf(marsRadius/1e3*x,marsRadius/1e3*y,marsRadius/1e3*z)
    hold off
    view([0,0])
    axis equal tight off
    if saveFigures, saveas(F,'../../Presentation/figures/trajectory','epsc'), end
  
    %...Plot first and last orbits
    firstOrbit = rotatedTrajectory(1:locPresentationApoapses(1),:);
    lastOrbit = rotatedTrajectory(locPresentationApoapses(end)+83025:end,:);

    F = figure('rend','painters','pos',figSizeWideHAR);
    hold on
    plot3(firstOrbit(:,1),firstOrbit(:,2),firstOrbit(:,3),'LineWidth',1.5,'Color',[0.85,0.325,0.098])
    plot3(lastOrbit(:,1),lastOrbit(:,2),lastOrbit(:,3),'LineWidth',1.5,'Color',[0.494,0.184,0.556])
    surf(marsRadius/1e3*x,marsRadius/1e3*y,marsRadius/1e3*z)
    hold off
    view([0,0])
    axis equal tight off
    if saveFigures, saveas(F,'../../Presentation/figures/configuration','epsc'), end
    
    %...Plot last orbit
    F = figure('rend','painters','pos',figSizeSmall);
    hold on
    plot3(lastOrbit(:,1),lastOrbit(:,2),lastOrbit(:,3),'LineWidth',1.5,'Color',[0.494,0.184,0.556])
    surf(marsRadius/1e3*x,marsRadius/1e3*y,marsRadius/1e3*z)
    hold off
    view([0,0])
    axis equal tight off
    if saveFigures, saveas(F,'../../Presentation/figures/final','epsc'), end
end

%% Plot Control Results

if loadDependent && loadRotational
    styles = {'-','--','-.'};
    aeroAngles = dependentVariables(:,7:9);
    
    F = figure('rend','painters','pos',figSizeWideLAR);
    hold on
    for i = 1:3
        plot(interpolatedTime,aeroAngles(:,i),'LineWidth',1.75,'LineStyle',styles{i})
    end
    hold off
    xlabel(timeLabel)
    ylabel('Aerodynamic Angle [deg]')
    set(gca,'FontSize',19)
    grid on
    legend('Attack','Side-slip','Bank','Location','Best')
    if saveFigures, saveas(F,'../../Presentation/figures/aerodynamic','epsc'), end
end

if loadRotational
    F = figure('rend','painters','pos',figSizeWideLAR);
    for i = 1:3
        subplot(1,3,i)
        plot(interpolatedTime,rotationalPropagatedResults(:,i+5),'LineWidth',1.25)
        xlabel(timeLabel)
        ylabel(rotationLabels{i+5})
        set(gca,'FontSize',15)
        grid on
    end
    if saveFigures, saveas(F,'../../Presentation/figures/rotation','epsc'), end
end

%% Periapsis Corridor

if loadEstimated
    if ~isempty(actualPeriapses)
        %...Plot periapses
        F = figure('rend','painters','pos',figSizeWideLAR);
        hold on
        scatter(periapsesLocs,actualPeriapses,50)
        scatter(upManeuverLoc,upManeuverValue,50,'^','filled')
        scatter(downManeuverLoc,downManeuverValue,50,'v','filled')
        plot(periapsisCorridors(2:end,1),periapsisCorridors(1:end-1,2)/1e3,'LineWidth',1.25,'LineStyle','--')
        plot(periapsisCorridors(2:end,1),periapsisCorridors(1:end-1,3)/1e3,'LineWidth',1.25,'LineStyle','-.')
        scatter(periapsesLocsFail,actualPeriapsesFail,50,[0,0.447,0.741],'o','filled')
        hold off
        xlabel('Orbit Number [-]')
        ylabel('Altitude [km]')
        xlim([0,210])
        set(gca,'FontSize',15)
        [L,icons] = legend('Periapsis','Up Maneuver','Down Maneuver','Lower Bound','Upper Bound','Location','Best');
        for i = 1:12
            if i >= 7 && i < 9
                icons(i).Children.MarkerSize = 10;
            end
        end
        grid on
        if saveFigures
            saveas(F,'../../Presentation/figures/corridor','epsc')
            F.Position = figSizeSmall;
            xlim([0,8])
            saveas(F,'../../Presentation/figures/corridor_walkin','epsc')
            xlim([145,210])
            L.Location = 'northwest';
            saveas(F,'../../Presentation/figures/corridor_walkout1','epsc')
            xlim([185,195])
            legend off
            saveas(F,'../../Presentation/figures/corridor_walkout2','epsc')
        end
    end
end