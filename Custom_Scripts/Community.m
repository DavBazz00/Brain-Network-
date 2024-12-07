% Load the coactivation matrix and coordinates
data = load('/home/davbaz/MATLAB/Toolboxes/BCT/2019_03_03_BCT/data_and_demos/Coactivation_matrix.mat');
Cij = data.Coactivation_matrix; % Adjacency matrix
Coord = data.Coord;             % Coordinates (x, y, z)

% Number of times to repeat community detection algorithm
num_iter = 100;

% Number of nodes
n_nodes = length(Cij);

% Empty array for storing the community labels
Ci = zeros(n_nodes, num_iter);

% Run the community detection algorithm num_iter times
for iter = 1:num_iter
    Ci(:, iter) = community_louvain(Cij);
end

% Calculate the module coassignment matrix
Coassignment = agreement(Ci) / num_iter;

% Perform consensus clustering
thr = 0.5;
cicon = consensus_und(Coassignment, thr, num_iter); % Final community assignments

% Visualize the coassignment matrix
figure;
imagesc(Coassignment);
colorbar;
title('Coassignment Matrix');
xlabel('Node Index');
ylabel('Node Index');

% Use adjacency_plot_und for visualization
[X, Y, Z] = adjacency_plot_und(Cij, Coord);

% 3D Visualization of Communities
figure;
% Plot edges
plot3(X, Y, Z, 'k-', 'LineWidth', 0.5); % Edges as black lines
hold on;
% Plot nodes colored by community
num_communities = max(cicon); % Determine the number of communities
cmap = lines(num_communities); % Generate a distinct colormap
node_colors = cmap(cicon, :); % Assign colors to nodes based on communities
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, node_colors, 'filled'); % Nodes as filled circles
title('3D Brain Network - Community Visualization');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
grid on;
hold off;
