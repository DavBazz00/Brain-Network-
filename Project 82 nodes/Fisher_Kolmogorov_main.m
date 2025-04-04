clear; close all; clc;
%% 1. Load file

edgeFile = '/home/davbaz/MATLAB/Custom_Scripts/Project2/Main/Edge.csv';
nodeFile = '/home/davbaz/MATLAB/Custom_Scripts/Project2/Main/Node.csv';

[A, CoordTable, Coord] = AdjacencyMatrix(edgeFile, nodeFile);

%% 2. Infection propagation
% Model parameters
diffusion = 5e-4; 
a = 0.5;         % Logistic growth rate parameter
dt = 0.4;        % Time step size in years
num_steps = 100;

% Call the infection simulation function
[t_sol, c_sol] = FK_propagation(A, CoordTable,diffusion, a, dt, num_steps);

%% 3. Plot Average Infection Concentration Over Time
plotAvgInfectionConcentration(t_sol, c_sol, 'Average concentration');

%% 4. Plot average concentratio by brain region overtime
plotRegionalConcentration(t_sol, c_sol, CoordTable, 'Average concentration per region');

%% 5. Infection propagation dynamic
edge_reduction = 0.20;
[t_sol_dinamic, c_sol_dinamic] = FK_propagation_dynamic(A, CoordTable, diffusion, a, dt, num_steps, edge_reduction);

%% 6. Plot Average Infection Concentration Over Time
plotAvgInfectionConcentration(t_sol_dinamic, c_sol_dinamic, 'Average concentration with dynamic edge variation');

%% 7. Comparison between first and after dynamic variation
plotComparisonAvgConcentration(t_sol, c_sol, t_sol_dinamic, c_sol_dinamic, 'Average Infection Concentration Comparison');

%% 8. Plot average concentration with different values of alpha
plotMultiAlpha_FK(A, CoordTable, diffusion, dt, num_steps);