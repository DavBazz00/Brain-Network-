% Load the coactivation matrix and coordinates
data = load('/home/davbaz/MATLAB/Toolboxes/BCT/2019_03_03_BCT/data_and_demos/Coactivation_matrix.mat');
Cij = data.Coactivation_matrix; % Adjacency matrix
Coord = data.Coord;             % Coordinates (x, y, z)

% Threshold the matrix to include only significant connections
threshold = 0.07; % Adjust as needed
Cij_thresh = threshold_absolute(Cij, threshold); % BCT function to apply threshold

% Normalize the weights in the thresholded adjacency matrix
Cij_norm = Cij_thresh / max(Cij_thresh(:)); % Normalize all weights to [0, 1]

% Convert the normalized weights to distances
Cost = weight_conversion(Cij_norm, 'lengths'); % BCT function to convert weights to distances

% Calculate the weighted distance matrix
D_wei = distance_wei(Cost); % BCT function to calculate the shortest path matrix


% Use adjacency_plot_und to get edge data for visualization
[X, Y, Z] = adjacency_plot_und(Cij_norm, Coord);

% Initialize the dynamic state
n_nodes = size(Cij_norm, 1);
state = zeros(n_nodes, 1); % 0: Susceptible (Blue), 1: Infected (Red)
initial_node = randi(n_nodes); % Randomly pick an initial infected node
state(initial_node) = 1; % Set the initial infected node

% Number of transitions and infection probability
n_steps = 10;
p_infect = 0.085; % Probability of infection per connection

% Visualize the initial state
figure;
% Plot edges
plot3(X, Y, Z, 'k-', 'LineWidth', 0.5); % Edges as black lines
hold on;
% Plot nodes
node_colors = [0 0 1] .* (1 - state) + [1 0 0] .* state; % Blue for susceptible, red for infected
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, node_colors, 'filled');
title('3D Brain Network - Dynamic States');
xlabel('X Coordinate');
ylabel('Y Coordinate');
zlabel('Z Coordinate');
grid on;

% Update the states dynamically using the normalized weighted graph
for t = 1:n_steps
    % Identify newly infected nodes
    infected_nodes = find(state == 1); % Indices of infected nodes
    susceptible_nodes = find(state == 0); % Indices of susceptible nodes
    susceptible_nodes = susceptible_nodes(:); % Force column vector
   
    for s = susceptible_nodes'
        % Check if the node is connected to any infected node using normalized weighted distances
        connections = D_wei(s, infected_nodes); % Weighted distances to infected nodes
        probabilities = exp(-connections); % Higher connection weight -> higher infection probability
        if any(rand(size(probabilities)) < probabilities)
            state(s) = 1; % Transition to infected
        end
    end

    % Update the visualization
    cla; % Clear previous plot
    % Re-plot edges
    plot3(X, Y, Z, 'k-', 'LineWidth', 0.5); % Edges as black lines
    hold on;
    % Re-plot nodes
    node_colors = [0 0 1] .* (1 - state) + [1 0 0] .* state; % Blue for susceptible, red for infected
    scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, node_colors, 'filled');
    title(['3D Brain Network - Dynamic States (Step ', num2str(t), ')']);
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    zlabel('Z Coordinate');
    grid on;
    pause(1); % Pause for visualization
end

