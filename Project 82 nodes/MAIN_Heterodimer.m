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
dt_treat    = 0.4;    % passo per simulateHeterodimerPropagationWithTreatment
num_steps_t = 100;    % per coprire ad esempio 40 anni
edge_reduction = 0.20;

[t_base, p_base, pt_base, t_treat, p_treat, pt_treat] = ...
    simulateHeterodimerPropagationWithTreatment( ...
        A, CoordTable, ...
        k0, k1, ktilde1, k12, diffusion_coeff, ...
        dt_treat, num_steps_t, t_switch, edge_reduction ); 

plotPropagationWithTreatment( ...
    t_base,   pt_base, ...
    t_treat,  pt_treat, ...
    t_switch, ...
    'Baseline vs. Treatment (post-switch) Heterodimer Propagation' );

%% 11. Combined Aging + Annual Treatment (Heterodimer)
dt_aging    = 0.4;   % stesso passo aging del dinamico
dt_treat    = 1;     % annuale per cura
T_end       = 40;    % anni totali
t_switch    = 10;     % inizio cura  
aging_red   = 0.20;  % fattore moltiplicativo aging
cure_red    = 0.20;  % fattore moltiplicativo cura

[t_comb_hd, pt_comb_hd] = HeterodimerInfection_combined( ...
    A, CoordTable, ...
    k0, k1, ktilde1, k12, diffusion_coeff, ...
    dt_aging, dt_treat, T_end, t_switch, ...
    aging_red, cure_red );

plotAvgInfectionConcentration( ...
    t_comb_hd, pt_comb_hd, ...
    'Average concentration with Aging + Annual Treatment (Heterodimer)' );

%% 12. Comparative plot (Heterodimer: Baseline / Aging / Cure / Aging+Cure)

% — riutilizzo dei risultati già calcolati —
% section 2:    [t_sol,           p_sol,           pt_sol]
% section 5:    [t_sol_dynamic,   p_sol_dynamic,   pt_sol_dynamic]
% section 10:   [t_baseline,      p_baseline,      pt_baseline, ...
%                t_treat,         p_treat,         pt_treat]
% section 11:   [t_comb_hd,       pt_comb_hd]

% concentrazioni medie
avg_baseline = mean(pt_sol,         2);   % no aging, no cure
avg_aging    = mean(pt_sol_dynamic, 2);   % aging only
avg_cure     = mean(pt_treat,       2);   % cure only
avg_combined = mean(pt_comb_hd,     2);   % aging + cure

% vettori tempo
t_base   = t_sol;
t_aging  = t_sol_dynamic;
t_cure   = t_treat;
t_comb   = t_comb_hd;

% plot
figure; hold on;
plot(t_base,  avg_baseline, 'b-', 'LineWidth',2);
plot(t_aging, avg_aging,    'g-', 'LineWidth',2);
plot(t_cure,  avg_cure,     'r-', 'LineWidth',2);
plot(t_comb,  avg_combined, 'm-', 'LineWidth',2);
xline(t_switch,'--k','Switch time','LineWidth',1.5);

xlabel('Time (years)');
ylabel('Average Misfolded-Protein Concentration');
legend( ...
  '1) Baseline', ...
  '2) Aging only', ...
  '3) Cure only', ...
  '4) Aging + Cure', ...
  'Switch time', ...
  'Location','best' );
grid on; hold off;

