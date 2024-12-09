%% Community Detection in Brain Network Using modularity_und
% This script detects communities within a brain network using the modularity_und
% function and visualizes the results in both 3D space and as a sorted adjacency matrix.

%% 1. Clear Workspace and Initialize
clear; close all; clc;

%% 2. Load Data
% Specify the path to your data file
data_path = '/home/davbaz/MATLAB/Toolboxes/BCT/2019_03_03_BCT/data_and_demos/Coactivation_matrix.mat';

% Load the coactivation matrix and coordinates
data = load(data_path);
Cij = data.Coactivation_matrix; % Adjacency matrix (NxN)
Coord = data.Coord;             % Coordinates (Nx3)

%% 3. Preprocess Adjacency Matrix (Optional)

% For this script, we'll use the original adjacency matrix
A = Cij;

%% 4. Community Detection Using modularity_und

% 4.1. Define the Resolution Parameter Gamma
gamma = 1; % Classic modularity. 
% >1 detects smaller communities, <1 detects larger communities

% 4.2. Run modularity_und to Detect Communities
[Ci, Q] = modularity_und(A, gamma);

% 4.3. Determine the Number of Communities
num_communities = max(Ci);
fprintf('Number of Communities Detected: %d\n', num_communities);
fprintf('Modularity Score: %.4f\n', Q);

%% 5. Visualization

% 5.1. 3D Visualization of Communities

% Generate a distinct colormap for communities
cmap = lines(num_communities); % 'lines' colormap provides distinguishable colors

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
disp('Community detection and visualization complete.');
