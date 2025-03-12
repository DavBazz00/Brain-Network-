%% Simulation of Infection Propagation in Human Brain Network
% This script models the spread of an infection within a human brain network
% using a reaction-diffusion model. It leverages the Brain Connectivity Toolbox
% for network analysis and visualization.

clear; close all; clc;

%% 1.Load Data
L = load('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/Brain-Network--main/Laplacian.csv'); % Laplacian
CoordTable = readtable('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/Brain-Network--main/Coordinates.csv');   % Coordinates

% Extract numeric coordinates
Coord = table2array(CoordTable(:, 4:6));  % X, Y, Z are in columns 4, 5, 6

% Find the Adjacency matrix from the Laplacian
A = -L;
A(A < 0) = 0;

% Remove small values from the Adjacency matrix (only for graph clarity)
A_tresh = threshold_absolute(A, 0.1);

%% 2. 3D Brain Map
figure;

% ====================== Parameters to Adjust ======================
edgeLineWidth = 0.1;     % Reduce for thinner edges (original was 0.5)
nodeSize = 100;           % Increase for greater node size (original was 50)
% ===================================================================

% Plot edges first
[X, Y, Z] = adjacency_plot_und(A_tresh, Coord);
plot3(X, Y, Z, 'y-', 'LineWidth', edgeLineWidth); 
hold on;

% Prepare node colors based on brain regions
brainRegions = unique(CoordTable{:, 7});  % Get unique region names
colors = lines(length(brainRegions));     % Color palette

% Initialize array to store scatter objects for legend
hScatter = gobjects(length(brainRegions), 1); 

% Plot nodes with different colors
for i = 1:length(brainRegions)
    region = brainRegions{i};
    isRegion = strcmp(CoordTable{:, 7}, region);
    
    % Store scatter handle for legend
    hScatter(i) = scatter3(Coord(isRegion,1), Coord(isRegion,2), Coord(isRegion,3),...
        nodeSize, colors(i,:), 'filled', 'DisplayName', region);
end

% Finalize plot
title('3D Brain Network');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend(hScatter, 'Location', 'eastoutside'); % Use SCATTER handles for legend
grid on; 
axis equal;
hold off;

%% 3. Infection Model

%%%%%%%%%%% Initial Conditions %%%%%%%%%%
N = size(A, 1);
c0 = zeros(N, 1);

% Identify and seed ONLY the selected brain regions
infected_mask = strcmp(CoordTable{:, 3}, 'entorhinal') | ...
                strcmp(CoordTable{:, 3}, 'temporalpole'); % as said in the paper

% Assign initial concentration (0.1) to all selected nodes
c0(infected_mask) = 0.1;
%%%%%%%%%%%%%%%%% END initial conditions %%%%%%%%%%%%%%%%%%%

% Model parameters
a = 0.5;        % Logistic growth rate parameter (adjust as needed)
dt = 0.4;       % Time step size in years
num_steps = 100; % Number of time steps
tmax = dt * num_steps; % Total simulation time

% Define the ODE Function (Reaction-Diffusion Model)
% dc/dt = - L*c + a*c.*(1 - c)
f = @(t, c) -L*c*0.05 + a*c.*(1-c);

%% 4. Simulate Infection Propagation

% Time vector for implicit time integration
t_sol = linspace(0, tmax, num_steps + 1);

% Preallocate solution matrix
c_sol = zeros(length(t_sol), N);
c_sol(1, :) = c0;

disp('Starting implicit time integration...');
try
    for i = 2:length(t_sol)
        t = t_sol(i);
        temp = c_sol(i-1, :)';              % Convert to column vector
        temp = temp + dt * f(t, temp);        % Compute update
        c_sol(i, :) = temp';                  % Convert back to row vector if needed
    end

    disp('Implicit time integration completed successfully.');
catch ME
    rethrow(ME);
end

%% 5. Plot Infection Dynamics

% Spatial Visualization of Initial and Final Infection Distribution

figure;

% Define the common color limits for both subplots (0 to 1) and a colormap
common_clims = [0 1];

% Plot the initial infection distribution
subplot(1,2,1);
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, c0, 'filled');
title('Initial Infection Distribution');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on;
% Apply the common color limits and colormap
clim(common_clims);

% Plot the final infection distribution
subplot(1,2,2);
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, c_sol(end, :), 'filled');
title(['Final Infection Distribution at t = ' num2str(tmax)]);
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on;
% Apply the common color limits and colormap to this subplot as well
clim(common_clims);

% Single colorbar for both subplots
hcb = colorbar('Position', [0.92 0.15 0.02 0.7]);
ylabel(hcb, 'Infection Concentration');

%% 6. Plot Average Infection Concentration Over Time

% Compute the average infection concentration at each time point
avg_conc = mean(c_sol, 2);

% Create a new figure for the plot
figure;
plot(t_sol, avg_conc, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); % MATLAB default blue
xlabel('Time');
ylabel('Average Infection Concentration');
title('Average Infection Concentration Over Time');
grid on;
hold off;

%% 7. Compute Average Infection Concentration for Selected Brain Regions Over Time

% Define custom groups for each region based on labels in CoordTable{:,3}
group_temporal = {'inferiortemporal', 'middletemporal', 'Right-Hippocampus', 'Left-Hippocampus'};
group_frontal  = {'precentral', 'caudalmiddlefrontal'};
group_parietal = {'superiorparietal', 'supramarginal', 'precuneus', 'insula'};
group_occipital = {'lateraloccipital', 'cuneus'};

% Combine the groups into a cell array and define region names
region_groups = {group_temporal, group_frontal, group_parietal, group_occipital};
region_names = {'Temporal', 'Frontal', 'Parietal', 'Occipital'};
num_groups = length(region_names);

% Define colors for plotting (same order as region_names)
region_colors = [0 1 0;    % Green for Temporal
                 1 0 0;    % Red for Frontal
                 1 0.5 0;  % Orange for Parietal
                 0 0 1];   % Blue for Occipital

% Initialize a matrix to store average concentrations
% Rows correspond to time points (length(t_sol)), columns to region groups
avg_conc_desired = zeros(length(t_sol), num_groups);

% Compute the average concentration for each region group at each time point
for r = 1:num_groups
    % Get the list of labels for this group
    group = region_groups{r};
    
    % Create a mask: true for nodes whose label in CoordTable{:,3} matches any label in the group.
    % Using lower-case conversion for robustness.
    nodes_in_region = ismember(lower(CoordTable{:,3}), lower(group));
    
    % Remove restriction to only left nodes for the Frontal group
    
    % If any nodes are found in the group, compute the mean over those nodes at each time point
    if any(nodes_in_region)
        avg_conc_desired(:, r) = mean(c_sol(:, nodes_in_region), 2);
    else
        avg_conc_desired(:, r) = 0;  % If no nodes found, set average to 0
    end
end

% Plot the results
figure('Name', 'Average Infection Concentration for Selected Brain Regions Over Time', 'NumberTitle', 'off');
hold on;
for r = 1:num_groups
    plot(t_sol, avg_conc_desired(:, r), 'LineWidth', 2, 'Color', region_colors(r, :));
end
hold off;
xlabel('Time');
ylabel('Average Infection Concentration');
title('Average Infection Concentration per Selected Brain Region Over Time');
legend(region_names, 'Location', 'bestoutside');
grid on;



%% 8. Create Video of Infection Propagation on 3D Brain Map

% Setup VideoWriter
videoFilename = 'InfectionPropagation.avi';
v = VideoWriter(videoFilename);
v.FrameRate = 5;  % Adjust frame rate as desired
open(v);

figure('Color','w'); % Create a new figure with white background

% Loop over each time step in the simulation
for i = 1:length(t_sol)
    clf;  % Clear current figure
    
    hold on;
    % Plot the edges (network structure) in thin yellow lines
    plot3(X, Y, Z, 'y-', 'LineWidth', edgeLineWidth);
    
    % Get infection concentration for each node at current time step
    infectionLevel = c_sol(i, :);
    
    % Plot nodes with marker color determined by infection level
    h = scatter3(Coord(:,1), Coord(:,2), Coord(:,3), ...
                 nodeSize, infectionLevel, 'filled');
    
    % Set colormap and color limits (assume infection is between 0 and 1)
    colormap(jet);
    clim([0 1]);
    colorbar; % Optional: show colorbar to indicate infection intensity
    
    % Label the plot
    title(sprintf('Infection Propagation at t = %.2f', t_sol(i)));
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on;
    view(3);  % Ensure 3D view
    
    % Capture frame and write to video
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    hold off;
end

% Close the video file
close(v);
disp(['Video saved as ', videoFilename]);

%%9. selezione nodi in base alla centralità nella matrice L
%vogliamo per ogni sottomatrice definita in base alla tipologia di nodo, come ad esempio frontale, occipatle, etc i primi 10 nodi più centrali
%per fare ciò dobbiamo calcolare la centralità di ogni nodo e poi selezionare i primi 10 nodi più centrali
% Calculate node centrality using eigenvector centrality
centrality_scores = eigenvector_centrality_und(A);

% Initialize a cell array to store the top 10 central nodes for each region
top_central_nodes = cell(num_groups, 1);

% Loop through each region group to find the top 10 central nodes
for r = 1:num_groups
    % Get the list of labels for this group
    group = region_groups{r};
    
    % Create a mask: true for nodes whose label in CoordTable{:,3} matches any label in the group.
    nodes_in_region = ismember(lower(CoordTable{:,3}), lower(group));
    
    % If any nodes are found in the group, compute the centrality for those nodes
    if any(nodes_in_region)
        % Extract centrality scores for nodes in the region
        region_centrality_scores = centrality_scores(nodes_in_region);
        
        % Get the indices of the top 10 central nodes in the region
        [~, sorted_indices] = sort(region_centrality_scores, 'descend');
        top_indices = sorted_indices(1:min(10, length(sorted_indices)));
        
        % Store the indices of the top central nodes
        top_central_nodes{r} = find(nodes_in_region);
        top_central_nodes{r} = top_central_nodes{r}(top_indices);
    else
        top_central_nodes{r} = [];  % If no nodes found, store an empty array
    end
end

% Display the top 10 central nodes for each region
for r = 1:num_groups
    fprintf('Top central nodes for %s region:\n', region_names{r});
    disp(top_central_nodes{r}');
end

