clear; close all; clc;

%% 1.Load Data
A = load('/home/davbaz/MATLAB/Custom_Scripts/Project2/Main/Edge.csv'); % Adjacency
A(1, :) = []; % Remove the first row

CoordTable = readtable('/home/davbaz/MATLAB/Custom_Scripts/Project2/Main/Node.csv'); % Coordinates
% Extract numeric coordinates
Coord = table2array(CoordTable(:, 1:3));

% Degree matrix
D = diag(sum(A, 2));
% Laplacian matrix
L = D - A;

%% 3. Infection Model

%%%%%%%%%%% Initial Conditions %%%%%%%%%%
N = size(A, 1);
c0 = zeros(N, 1);

% Identify and seed ONLY the selected brain regions
infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal'); % as said in the paper

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
f = @(t, c) -L*c*5e-4 + a*c.*(1-c);

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


%% 7. 
% Define regions of interest and colors
region_names = {'Temporal', 'Frontal', 'Parietal', 'Occipital'};
region_colors = [0 1 0;    % Green for Temporal
                 1 0 0;    % Red for Frontal
                 1 0.5 0;  % Orange for Parietal
                 0 0 1];   % Blue for Occipital

num_groups = length(region_names);

% Extract region labels from CoordTable
all_regions = CoordTable{:, 4};

% Initialize figure
figure;
hold on;

% Loop through each region and compute average concentration over time
for i = 1:num_groups
    region = region_names{i};
    
    % Logical index for nodes in this region
    idx = strcmp(all_regions, region);
    
    % Extract the relevant columns from c_sol and compute the mean
    avg_conc_region = mean(c_sol(:, idx), 2);
    
    % Plot
    plot(t_sol, avg_conc_region, 'LineWidth', 2, 'Color', region_colors(i, :));
end

% Formatting
xlabel('Time');
ylabel('Average Misfolded Protein Concentration');
title('Average Concentration by Brain Region Over Time');
legend(region_names, 'Location', 'northeast');
grid on;
hold off;
