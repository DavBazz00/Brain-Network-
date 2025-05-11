clear; close all; clc;

%% 1. Load file
edgeFile = 'Edge.csv';
nodeFile = 'Node.csv';

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


%% 3.5 3D plot of the infection distribution both at time t=0 and at a desired time
time_desired = 20;
plotInfectionDynamics3D(Coord, CoordTable, pt_sol, t_sol, time_desired);


%% 4. Plot Average Concentration Over Time for Specific Node Pairs
plotSpecificNodePairs(t_sol, pt_sol, CoordTable, 'Average Concentration Over Time for Specific Node Pairs');


%% 5. Run the simulation with dynamic edge variation
edge_reduction = 0.20;
[t_sol_dynamic, p_sol_dynamic, pt_sol_dynamic] = HeterodimerInfection_dynamic(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt, num_steps);


%% 6. Plot Average Infection Concentration Over Time
plotAvgInfectionConcentration(t_sol_dynamic, pt_sol_dynamic, 'Average concentration with dynamic edge variation');


%% 7. Comparison between first and after dynamic variation
plotComparisonAvgConcentration(t_sol, pt_sol, t_sol_dynamic, pt_sol_dynamic, 'Average Infection Concentration Comparison');


%% 8. Plot average concentration with different values of k12
plotMultiK12_Heterodimer(A, CoordTable, k0, k1, ktilde1, diffusion_coeff, dt, num_steps);


%% 9. Baseline vs. Cure
% parameters for cure scenario
dt_aging    = 0.4;
dt_treat    = 1;
num_steps_t = 40;
edge_reduction = 0.20;
t_switch    = 15;
t_end       = 40;

% 9.1 run baseline (static network) up to t_end
num_steps_base = round(t_end / dt_aging);
[t_base, ~, pt_base] = HeterodimerInfection(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff,dt_aging, num_steps_base);

% 9.2 run cure scenario (treatment)
[t_treat, ~, pt_treat] = HeterodimerInfection_treatment(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, ...
    edge_reduction, dt_aging, dt_treat, t_switch, t_end);

% 9.3 plot
plotPropagationWithTreatment(t_base, pt_base,t_treat, pt_treat,t_switch, 'Baseline vs. Cure (Heterodimer)' );


%% 10. Comparison of Four Scenarios (re-computed for clarity)
plotFourScenariosHeterodimer( ...
    A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, ...
    edge_reduction, dt_aging, dt_treat, t_switch, t_end );



