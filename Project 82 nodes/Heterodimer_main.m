clear; close all; clc;
%% 1. Load file

edgeFile = '/home/davbaz/MATLAB/Custom_Scripts/Project2/Main/Edge.csv';
nodeFile = '/home/davbaz/MATLAB/Custom_Scripts/Project2/Main/Node.csv';

[A, CoordTable, Coord] = AdjacencyMatrix(edgeFile, nodeFile);

%% 2. Infection Propagation
% Set model parameters for the heterodimer infection model
k0 = 1.0;       % Production rate of healthy protein
k1 = 0.5;       % Clearance rate of healthy protein
ktilde1 = 0.5;  % Clearance rate of misfolded protein
k12 = 0.5;      % Conversion rate from healthy to misfolded protein
diffusion_coeff = 5e-4;  % Diffusion coefficient
dt = 0.4;       % Time step size in years
num_steps = 100; % Number of time steps

% Run the heterodimer infection
[t_sol, p_sol, pt_sol] = HeterodimerInfection(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps);

%% 3. Plot Average Infection Concentration Over Time
plotAvgInfectionConcentration(t_sol, pt_sol, 'Average concentration');

%% 4. Plot average concentratio by brain region overtime
plotRegionalConcentration(t_sol, pt_sol, CoordTable, 'Average concentration per region');

%% 5. Run the simulation with dynamic edge variation
edge_reduction = 0.20;
[t_sol_dynamic, p_sol_dynamic, pt_sol_dynamic] = HeterodimerInfection_dynamic(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt, num_steps);

%% 6. Plot Average Infection Concentration Over Time
plotAvgInfectionConcentration(t_sol_dynamic, pt_sol_dynamic, 'Average concentration with dynamic edge variation');

%% 7. Comparison between first and after dynamic variation
plotComparisonAvgConcentration(t_sol, pt_sol, t_sol_dynamic, pt_sol_dynamic, 'Average Infection Concentration Comparison');

%% 8. Plot average concentration with different values of k12
plotMultiK12_Heterodimer(A, CoordTable, k0, k1, ktilde1, diffusion_coeff, dt, num_steps);