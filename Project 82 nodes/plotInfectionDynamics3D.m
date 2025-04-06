function plotInfectionDynamics3D(Coord, CoordTable, c_sol, t_sol, plotTime)
%PLOTINFECTIONDYNAMICS3D Visualizes the spatial distribution of the initial and
% infection at a specified time.
%
%   plotInfectionDynamics3D(Coord, CoordTable, c_sol, t_sol, plotTime)
%
%   INPUTS:
%       Coord     - Matrix of node coordinates (Nx3).
%       CoordTable- Table containing node data, including region information in 
%                   column 4.
%       c_sol     - Matrix of infection concentrations over time (T x N).
%       t_sol     - Time vector corresponding to c_sol (length T).
%       plotTime  - Desired time at which to plot the infection distribution.
%
%   This function internally computes the initial infection distribution (c0)
%   by seeding nodes in the 'Entorhinal' region with 0.1 and all others with 0.
%   It then creates a figure with two subplots:
%       1. The initial infection distribution.
%       2. The infection distribution at the time closest to plotTime.
%
%   A custom colormap from light grey to bordeaux red is applied, and the nodes 
%   initially infected are highlighted with a thin black circle.
%
%   Example usage:
%       plotInfectionDynamics3D(Coord, CoordTable, c_sol, t_sol, 20);
%

    %% Compute initial infection distribution (c0)
    N = size(Coord, 1);
    c0 = zeros(N, 1);
    infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
    c0(infected_mask) = 0.1;
    
    %% Create a new figure and define common color limits
    figure;
    common_clims = [0 1];
    
    %% Define custom colormap: from light grey to bordeaux red
    nColors = 256;
    startColor = [0.9, 0.9, 0.9];   % Light grey for 0
    endColor   = [0.5, 0, 0.13];      % Bordeaux red for 1
    customCmap = [linspace(startColor(1), endColor(1), nColors)', ...
                  linspace(startColor(2), endColor(2), nColors)', ...
                  linspace(startColor(3), endColor(3), nColors)'];
              
    %% Plot the initial infection distribution (subplot 1)
    subplot(1,2,1);
    scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, c0, 'filled');
    title('Initial Infection Distribution');
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on;
    clim(common_clims);
    colormap(customCmap);
    
    % Highlight nodes that are initially infected (c0 > 0) with a thin black circle
    infectedIdx = find(c0 > 0);
    circleRadius = 2;  % Adjust based on your coordinate scale
    theta = linspace(0, 2*pi, 50);
    hold on;
    for i = 1:length(infectedIdx)
        idx = infectedIdx(i);
        x_center = Coord(idx, 1);
        y_center = Coord(idx, 2);
        z_center = Coord(idx, 3);
        x_circle = x_center + circleRadius * cos(theta);
        y_circle = y_center + circleRadius * sin(theta);
        z_circle = repmat(z_center, size(theta));
        plot3(x_circle, y_circle, z_circle, 'k-', 'LineWidth', 1);
    end
    hold off;
    
    %% Find the index corresponding to the desired plotTime
    [~, idx_plot] = min(abs(t_sol - plotTime));
    
    %% Plot the infection distribution at the desired time (subplot 2)
    subplot(1,2,2);
    scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, c_sol(idx_plot, :), 'filled');
    title(['Infection Distribution at t = ' num2str(t_sol(idx_plot))]);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on;
    clim(common_clims);
    colormap(customCmap);
    
    % Add a single colorbar for both subplots
    hcb = colorbar('Position', [0.92 0.15 0.02 0.7]);
    ylabel(hcb, 'Infection Concentration');
end
