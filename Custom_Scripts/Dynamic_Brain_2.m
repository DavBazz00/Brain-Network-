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
a = 0.05;        % Logistic growth rate parameter (adjust as needed)
tmax = 100;      % Max simulation time

% The initial infected nodes are fixed in space (i.e. cluster of infected nodes), 
% and they are not random instead

% 5.3. Initial Conditions
N = size(A, 1);
c0 = zeros(N, 1);

% Initialize infection 
desired_coord = [0, 0, 0]; % Replace with desired coordinates
% Compute distances from each node to the desired coordinate
dists = pdist2(Coord, desired_coord); 

% Find the closest node to the desired coordinate
[min_val, min_ind] = min(dists);
seed_node = min_ind;

% Find neighbors of the seed node by checking which nodes have a nonzero edge from node 1.
neighbors = find(A(seed_node, :) > 0);
% Cluster of infected nodes
cluster_size = min(5, length(neighbors)); % Infect up to 5 neighbors, or fewer if less are available
infected_nodes = [seed_node, neighbors(1:cluster_size)]; 

% Set initial concentration of infection in this cluster to 0.1
c0(infected_nodes) = 0.1;

% 6.4. Define the ODE Function (Reaction-Diffusion Model)
% dc/dt = L*c + a*c.*(1 - c)
f = @(t, c) L * c + a * c .* (1 - c);

%% 6. Simulate Infection Propagation


% 6.1. Solve the ODE using ode15s (stiff solver)
disp('Starting ODE solver...');
try
    [t_sol, c_sol] = ode15s(f, [0 tmax], c0);
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

%% End of Script
disp('Simulation complete.');
