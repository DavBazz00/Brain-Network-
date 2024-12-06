% Load the coactivation matrix and coordinates
data = load('/home/davbaz/MATLAB/Toolboxes/BCT/2019_03_03_BCT/data_and_demos/Coactivation_matrix.mat');
Cij = data.Coactivation_matrix; % Adjacency matrix
Coord = data.Coord;             % Coordinates (x, y, z)

% Threshold the matrix to include only significant connections (optional)
threshold = 0.1; % Adjust as needed
Cij_thresh = Cij > threshold;

% Use adjacency_plot_und to get edge data for 3D visualization
[X, Y, Z] = adjacency_plot_und(Cij_thresh, Coord);

% Plot the edges in 3D
figure;
plot3(X, Y, Z, 'k-', 'LineWidth', 0.5); % Edges as black lines
hold on;

% Plot the nodes
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, 'b', 'filled'); % Nodes as blue dots
title('3D Brain Network');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
grid on;
hold off;
