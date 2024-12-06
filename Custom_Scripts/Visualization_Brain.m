% Load the coactivation matrix and coordinates
data = load('/home/davbaz/MATLAB/Toolboxes/BCT/2019_03_03_BCT/data_and_demos/Coactivation_matrix.mat');
Cij = data.Coactivation_matrix; % Adjacency matrix
Coord = data.Coord;             % Coordinates (x, y, z)

% Convert weights into costs for weighted distance calculations
Cost = 1 ./ Cij; % Transform weights into costs (inverse of connectivity)

% create a figure and use imagesc to visualize the matrix
f = figure; 
subplot(1,2,1)
title('Unweighted Graph');
imagesc(Cij);
grid on;
subplot(1,2,2)
title('Weighted Graph');
imagesc(Cost);
