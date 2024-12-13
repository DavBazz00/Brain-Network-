%% Simulation of Infection Propagation in Human Brain Network
% This script models the spread of an infection within a human brain network
% using a reaction-diffusion model. It leverages the Brain Connectivity Toolbox
% for network analysis and visualization.

%% 1. Clear Workspace and Initialize
clear; close all; clc;

%% 2. Load Data
% Load the coactivation matrix and coordinates
data = load('/home/davbaz/MATLAB/Toolboxes/BCT/2019_03_03_BCT/data_and_demos/Coactivation_matrix.mat');
Cij = data.Coactivation_matrix; % Adjacency matrix (NxN)
Coord = data.Coord;             % Coordinates (Nx3)

%% 3. Preprocess Adjacency Matrix
% Adjacency matrix 
Atemp = Cij;

% Thresholding: Choose between absolute or proportional thresholding
% Option 1: Absolute Thresholding
absolute_threshold = 0.1; % Adjust based on data
A_thr_abs = threshold_absolute(Atemp, absolute_threshold);

% Option 2: Proportional Thresholding
proportional_threshold = 0.1; % Preserve top 10% connections
A_thr_prop = threshold_proportional(Atemp, proportional_threshold);

% Choose the thresholded adjacency matrix to use
A = A_thr_prop; % Change to A_thr_abs if absolute thresholding is preferred


%% 4. Visualize the Network in 3D
% Use the pre-built adjacency_plot_und function to get edge data for 3D visualization
[X, Y, Z] = adjacency_plot_und(A, Coord);

% Plot the edges in 3D with green color
figure;
plot3(X, Y, Z, 'g-', 'LineWidth', 0.5); % Edges as green lines
hold on;

% Plot the nodes
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, 'b', 'filled'); % Nodes as blue dots

title('3D Brain Network');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
grid on;
axis equal;
hold off;

%% 5. Define the Infection Model

% 5.1. Compute Graph Laplacian
%   L = D - A
% where:
%   - A is the adjacency matrix of the graph.
%   - D is the degree matrix, a diagonal matrix where D_ii = sum over j of A_ij.

D = diag(sum(A, 2));
L = D - A;
L = sparse(L); % Convert to sparse for efficiency

% 5.2. Model Parameters
a = 0.005;        % Logistic growth rate parameter (adjust as needed)
tmax = 1e+10;      % Max simulation time

% The initial infected nodes are fixed in space (i.e. cluster of infected nodes), 
% and they are not random instead

% 5.3. Initial Conditions
N = size(A, 1);
c0 = zeros(N, 1);

% Initialize infection 
desired_coord = [-43, -26, -23];
% Compute distances from each node to the desired coordinate
dists = pdist2(Coord, desired_coord); 

% Sort all nodes based on their distance to desired_coord
[sorted_dists, sorted_idx] = sort(dists, 'ascend');

% Determine the number of nodes to infect
cluster_size = 5;
num_infected = min(5, cluster_size);

% Select the indices of the five closest nodes
infected_nodes = sorted_idx(1:num_infected);
% Set initial concentration of infection in this cluster to 0.1
c0(infected_nodes) = 0.1;


% 6.4. Define the ODE Function (Reaction-Diffusion Model)
% dc/dt = L*c + a*c.*(1 - c)
%f = @(t, c) L * c + a * c .* (1 - c); (MODIFICATION)
% Introduce a time scaling factor to slow down the dynamics
time_scaling_factor = 1.535e-10; 
% Redefine the ODE function to include the time scaling factor
f = @(t, c) time_scaling_factor * (L * c + a * c .* (1 - c));

%% 6. Simulate Infection Propagation

options = odeset(...
    'RelTol', 1e-4, ...
    'AbsTol', 1e-6, ...
    'Vectorized', 'on', ...
    'Stats', 'on' ...
);


% 6.1. Solve the ODE using ode15s
disp('Starting ODE solver...');
try
    [t_sol, c_sol] = ode15s(f, [0 tmax], c0, options);
    disp('ODE solver completed successfully.');
catch ME
    disp('Error during ODE solving:');
    disp(ME.message);
end

% 6.2. Verify Solution
% Ensure concentrations are within [0,1]
c_sol(c_sol < 0) = 0;
c_sol(c_sol > 1) = 1;

%% 7. Plot Infection Dynamics

% 7.1. Spatial Visualization of Initial and Final Infection Distribution

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

% To have a single colorbar for both subplots, create it at the figure level
hcb = colorbar('Position', [0.92 0.15 0.02 0.7]);
ylabel(hcb, 'Infection Concentration');

%% 8. Plot Average Infection Concentration Over Time

% 8.1. Compute the average infection concentration at each time point
avg_conc = mean(c_sol, 2);

% Create a new figure for the plot
figure;
plot(t_sol, avg_conc, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); % MATLAB default blue
xlabel('Time');
ylabel('Average Infection Concentration');
title('Average Infection Concentration Over Time');
grid on;
hold off;

%% 9. Compute Average Infection Concentration per Community Over Time
% Determine the Number of Communities
gamma = 1; % Classic modularity. 
% >1 detects smaller communities, <1 detects larger communities

%Run modularity_und to Detect Communities
[Ci, Q] = modularity_und(A, gamma);

num_communities = max(Ci);
fprintf('Number of Communities Detected: %d\n', num_communities);
fprintf('Modularity Score: %.4f\n', Q);

% Initialize a matrix to store average concentrations
% Rows correspond to time points, columns to communities
avg_conc_per_comm = zeros(length(t_sol), num_communities);

% Compute the average concentration for each community at each time point
for comm = 1:num_communities
    % Find nodes belonging to the current community
    nodes_in_comm = find(Ci == comm);
    
    % Compute the mean concentration for these nodes across all time points
    if ~isempty(nodes_in_comm)
        avg_conc_per_comm(:, comm) = mean(c_sol(:, nodes_in_comm), 2);
    else
        avg_conc_per_comm(:, comm) = 0; % Handle empty communities
    end
end


% Generate a distinct colormap for the communities
cmap = lines(num_communities); % 'lines' colormap provides distinguishable colors

% Create a new figure for the plot
figure('Name', 'Average Infection Concentration per Community Over Time', 'NumberTitle', 'off');
hold on;
for comm = 1:num_communities
    plot(t_sol, avg_conc_per_comm(:, comm), 'LineWidth', 2, 'Color', cmap(comm, :));
end
hold off;

% Customize the plot
xlabel('Time');
ylabel('Average Infection Concentration');
title('Average Infection Concentration per Community Over Time');
legend(arrayfun(@(x) sprintf('Community %d', x), 1:num_communities, 'UniformOutput', false), ...
       'Location', 'best');
grid on;

%% 9. Visualization of Communities in 3D Space

% Assign colors to nodes based on community assignments
node_colors = cmap(Ci, :);

% Plot the brain network with nodes colored by community
figure('Name', '3D Brain Network - Community Visualization', 'NumberTitle', 'off');
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 60, node_colors, 'filled'); % Nodes as filled circles
title('3D Brain Network - Community Visualization using modularity\_und');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
grid on;
axis equal;
hold off;

%% End of Script
disp('Simulation complete.');
