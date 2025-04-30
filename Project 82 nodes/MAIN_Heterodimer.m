clear; close all; clc;
addpath('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/project_82_withfunction/functions')

%% 1. Load file

edgeFile = '/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Edge.csv';
nodeFile = '/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Node.csv';

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

%% 4. Plot average concentratio by brain region overtime
plotRegionalConcentration(t_sol, pt_sol, CoordTable, 'Average concentration per region');

%% 4.5 Plot Average Concentration Over Time for Specific Node Pairs
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

%% 9. Plot average misfolded infection with starting point different nodes (Heterodimer)
plotSelectedNodes_Heterodimer(A, CoordTable, 1.0, 0.5, 0.5, 0.5, 5e-4, 0.4, 100);

%% 10. Simulate treatment scenario (Heterodimer)
t_switch    = 10;      % anno inizio trattamento
dt_aging    = 0.4; 
dt_treat    = 1;   % passo per simulateHeterodimerPropagationWithTreatment
num_steps_t = 40;    % per coprire ad esempio 40 anni
edge_reduction = 0.20;

[t_base, p_base, pt_base, t_treat, p_treat, pt_treat] = ...
    simulateHeterodimerPropagationWithTreatment( ...
        A, CoordTable, ...
        k0, k1, ktilde1, k12, diffusion_coeff, ...
        dt_aging, dt_treat, num_steps_t, t_switch, edge_reduction ); 

plotPropagationWithTreatment( ...
    t_base,   pt_base, ...
    t_treat,  pt_treat, ...
    t_switch, ...
    'Baseline vs. Treatment (post-switch) Heterodimer Propagation' );

%% 11. Simulazione treatment con nuova cura (Heterodimer, senza combined)

% Parametri
dt1           = 0.4;   % passo aging (anni)
dt2           = 1.0;   % passo annuale per cura
num_steps     = 40;    % anni totali
t_switch      = 10;    % anno di inizio trattamento
edge_reduction= 0.20;  % fattore di riduzione per la cura

% Lancio simulazione:
% - baseline non trattata è calcolata internamente via HeterodimerInfection
% - da t_switch in poi uso HeterodimerInfection_treatment_dynamic
[ t_baseline_hd, p_baseline_hd, pt_baseline_hd, ...
  t_treatment_hd, p_treatment_hd, pt_treatment_hd ] = ...
    simulateHeterodimerPropagationWithTreatment( ...
        A, CoordTable, ...
        k0, k1, ktilde1, k12, diffusion_coeff, ...
        dt1, dt2, num_steps, t_switch, edge_reduction );

% Plot della sola componente "infected" media (pt) post-switch
plotAvgInfectionConcentration( ...
    t_treatment_hd, pt_treatment_hd, ...
    'Average concentration with Aging + Annual Treatment (Heterodimer)' );

%% 12. Panoramica dei grafici generati (Heterodimer)

% 12.1 – Average concentration: Baseline
plotAvgInfectionConcentration( ...
    t_sol, pt_sol, ...
    '1) Average concentration – Baseline');

% 12.2 – Average concentration: Aging only (dynamic edge variation)
plotAvgInfectionConcentration( ...
    t_sol_dynamic, pt_sol_dynamic, ...
    '2) Average concentration – Aging only');

% 12.3 – Baseline vs Treatment (post-switch)
plotPropagationWithTreatment( ...
    t_base,   pt_base, ...
    t_treat,  pt_treat, ...
    t_switch, ...
    '3) Baseline vs. Treatment (post-switch)');

% 12.4 – Average concentration: Aging + Annual Treatment
plotAvgInfectionConcentration( ...
    t_treatment_hd, pt_treatment_hd, ...
    '4) Average concentration – Aging + Annual Treatment');

% 12.5 – Plot comparativo finale con tutte e quattro le curve
figure; hold on;
h1 = plot( t_sol,           mean(pt_sol,        2), 'b-', 'LineWidth',2 );
h2 = plot( t_sol_dynamic,   mean(pt_sol_dynamic,2), 'g-', 'LineWidth',2 );
h3 = plot( t_treat,         mean(pt_treat,      2), 'r-', 'LineWidth',2 );
h4 = plot( t_treatment_hd,  mean(pt_treatment_hd,2), 'm-', 'LineWidth',2 );
xline(t_switch, '--k', 'Switch time', 'LineWidth',1.5);

xlabel('Time (years)');
ylabel('Average Misfolded-Protein Concentration');
legend( ...
  [h1 h2 h3 h4], ...
  {'Baseline','Aging only','Treatment only','Aging + Treatment'}, ...
  'Location','best' );
grid on; hold off;

