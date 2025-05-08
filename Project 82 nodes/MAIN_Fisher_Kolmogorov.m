clear; close all; clc;
addpath('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/project_82_withfunction/functions')

%% 1. Load file

edgeFile = 'Edge.csv';
nodeFile = 'Node.csv';

[A, CoordTable, Coord] = AdjacencyMatrix(edgeFile, nodeFile);


%% 2. Infection propagation
% Model parameters
diffusion = 5e-4; 
a = 0.5;         % Logistic growth rate parameter, determine how fast the infection can grow
dt = 0.4;        % Time step size in years
num_steps = 100;

% Call the infection simulation function
[t_sol, c_sol] = FK_propagation(A, CoordTable,diffusion, a, dt, num_steps);

%% 3. Plot Average Infection Concentration Over Time
plotAvgInfectionConcentration(t_sol, c_sol, 'Average concentration');

%% 3.5 3D plot of the infection distribution both at time t=0 and at a desired time
time_desired = 20;
plotInfectionDynamics3D(Coord, CoordTable, c_sol, t_sol, time_desired);


%% 4. Plot average concentration by brain region overtime
plotRegionalConcentration(t_sol, c_sol, CoordTable, 'Average concentration per region');

%% 4.5 Plot Average Concentration Over Time for Specific Node Pairs
plotSpecificNodePairs(t_sol, c_sol, CoordTable, 'Average Concentration Over Time for Specific Node Pairs');


%% 5. Infection propagation dynamic
edge_reduction = 0.20;
[t_sol_dinamic, c_sol_dinamic] = FK_propagation_dynamic(A, CoordTable, diffusion, a, dt, num_steps, edge_reduction);

%% 6. Plot Average Infection Concentration Over Time
plotAvgInfectionConcentration(t_sol_dinamic, c_sol_dinamic, 'Average concentration with dynamic edge variation');

%% 7. Comparison between first and after dynamic variation
plotComparisonAvgConcentration(t_sol, c_sol, t_sol_dinamic, c_sol_dinamic, 'Average Infection Concentration Comparison');

%% 8. Plot average concentration with different values of alpha
plotMultiAlpha_FK(A, CoordTable, diffusion, dt, num_steps);


%% 9. Simulate treatment scenario
t_switch = 15.2;      % Time at which treatment starts (years), realistically could be 10y or more
dt1 = 0.4; dt2 = 1; num_steps = 40;
edge_reduction = 0.20;
[t_baseline, c_baseline, t_treatment, c_treatment] = simulateFKPropagationWithTreatment(A, CoordTable, diffusion, a, dt1, dt2, num_steps, t_switch, edge_reduction);

plotPropagationWithTreatment(t_baseline, c_baseline, t_treatment, c_treatment, t_switch, 'Baseline vs. Treatment (post-switch) FK Propagation');

%% 10. Simulazione treatment con nuova cura (senza usare più FK_propagation_combined, addio combined, da togliere)

% Parametri
dt1        = 0.4;      % passo prima dello switch (anni)
dt2        = 1.0;      % passo dopo lo switch (anni)
num_steps  = 40;       % numero totale di step a dt2
t_switch   = 10;       % anno di inizio trattamento
edge_reduction = 0.20; % fattore di riduzione per la cura

% Lancio la simulazione: 
% - baseline non trattata viene calcolata internamente via FK_propagation
% - post-switch usa FK_propagation_treatment_dynamic con nuova cura
[t_baseline, c_baseline, t_treatment, c_treatment] = ...
    simulateFKPropagationWithTreatment( ...
        A, CoordTable, diffusion, a, ...
        dt1, dt2, num_steps, ...
        t_switch, edge_reduction );

% (Opzionale) Plot di controllo della sola curva Aging+Cure
plotAvgInfectionConcentration(t_comb, c_comb, ...
    'Average conc. con nuova cura');

%% 11. Comparison of Four Scenarios (FK)

plotFourScenariosFK( ...
  A, CoordTable, diffusion, a, edge_reduction, ...
  dt_aging, dt_treat, t_switch, t_end );

