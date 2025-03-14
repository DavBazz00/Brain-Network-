%% Simulation of Infection Propagation in Human Brain Network
% This script models the spread of an infection within a human brain network
% using a reaction-diffusion model. It leverages the Brain Connectivity Toolbox
% for network analysis and visualization.

clear; close all; clc;
addpath('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/matlab/Custom_Scripts');

%% 1. Load Data
L = load('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/matlab/Laplacian.csv'); % Laplacian
CoordTable = readtable('//Users/Speranza/Desktop/VSCODE/matlab/biomedicina/matlab/Coordinates.csv');   % Coordinates

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
edgeLineWidth = 0.5;     % Reduce for thinner edges (original was 0.5)
nodeSize = 100;           % Increase for greater node size (original was 50)
% ===================================================================

% Plot edges first
% Define the adjacency_plot_und function if not available
if ~exist('adjacency_plot_und1', 'file')
    % Function not available
    warning('Function adjacency_plot_und not found.');
end


%function [X, Y, Z] = adjacency_plot_und(A, Coord)
%    [i, j] = find(A);
%    X = [Coord(i, 1) Coord(j, 1)]';
%    Y = [Coord(i, 2) Coord(j, 2)]';
%    Z = [Coord(i, 3) Coord(j, 3)]';
%end


[X, Y, Z] = adjacency_plot_und1(A_tresh, Coord);
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
% Estrarre le etichette dalla colonna 3 della tabella
% Eliminate "insula" group because it can be recasted in more then 1 group
group_temporal = {'inferiortemporal', 'middletemporal', 'Right-Hippocampus', 'Left-Hippocampus','fusiform','temporalpole','superiortemporal','transversetemporal','bankssts'};
group_frontal  = {'precentral', 'caudalmiddlefrontal','lateralorbitofrontal', 'rostralmiddlefrontal', 'superiorfrontal', 'parsopercularis', 'parstriangularis', 'parsorbitalis', 'rostralanteriorcingulate', 'caudalanteriorcingulate', 'frontalpole','medialorbitofrontal'};
group_parietal = {'superiorparietal', 'supramarginal', 'precuneus','inferiorparietal','postcentral','paracentral','isthmuscingulate'};
group_occipital = {'lateraloccipital', 'cuneus', 'pericalcarine', 'lingual'};

%It could be better to choose directly from column 7 of CoordTable (idk)

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

figure('Color','w', 'Position', [100, 100, 1120, 840]); % Create a new figure with white background and specified size

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

%% 9. selezione nodi in base alla centralità nella matrice L
% vogliamo per ogni sottomatrice definita in base alla tipologia di nodo, come ad esempio frontale, occipitale, etc i primi 10 nodi più centrali
% per fare ciò dobbiamo calcolare la centralità di ogni nodo e poi selezionare i primi 10 nodi più centrali

% Calcolare la centralità di grado per tutti i nodi
degree_centrality = sum(A, 2);

% Inizializzare una struttura per memorizzare i nodi più centrali per ciascuna regione
top_central_nodes = struct();

% Iterare su ciascun gruppo di regioni
for r = 1:num_groups
    % Ottenere l'elenco delle etichette per questo gruppo
    group = region_groups{r};
    
    % Creare una maschera: true per i nodi la cui etichetta in CoordTable{:,3} corrisponde a qualsiasi etichetta nel gruppo
    nodes_in_region = ismember(lower(CoordTable{:,3}), lower(group));
    
    % Estrarre le centralità dei nodi per questa regione
    centrality_in_region = degree_centrality(nodes_in_region);
    
    % Estrarre gli indici dei nodi in questa regione
    region_indices = find(nodes_in_region);
    
    % Ordinare i nodi per centralità in ordine decrescente e selezionare i primi 10
    [~, sorted_indices] = sort(centrality_in_region, 'descend');
    top_10_indices = region_indices(sorted_indices(1:min(10, length(sorted_indices))));
    
    % Memorizzare i nodi più centrali per questa regione
    top_central_nodes.(region_names{r}) = top_10_indices;
end

% Visualizzare i nodi più centrali per ciascuna regione
for r = 1:num_groups
    fprintf('Top 10 central nodes in %s region:\n', region_names{r});
    disp(CoordTable(top_central_nodes.(region_names{r}), :));
end
