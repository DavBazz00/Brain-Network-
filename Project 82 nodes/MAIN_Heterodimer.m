clear; close all; clc;
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

% --- 2) Definisci le etichette selezionate per ciascuna regione ---
temp_labels    = {'22L', '22R'};   % Region Temporale
frontal_labels = {'10L', '10R'};     % Region Frontale
parietal_labels= {'5L',  '5R'};       % Region Parietale
occipital_labels = {'19L', '19R'};    % Region Occipitale

% Combina le etichette (si tratta di 8 nodi in totale)
selected_labels = [temp_labels, frontal_labels, parietal_labels, occipital_labels];

% --- 3) Estrai le etichette dei nodi dalla CoordTable
% Nel modello eterodimero si assume che le etichette (regioni o nodi) siano nella colonna 5.
all_labels = cellstr(CoordTable{:,5});
selectedIdx = find(ismember(all_labels, selected_labels));

if length(selectedIdx) ~= 8
    error('Expected to find 8 nodes with the specified labels, but found %d.', length(selectedIdx));
end

% --- 4) Definisci la mappatura fissa dei colori per le 4 regioni ---
% I colori fissi (es. Temporale: Green, Frontale: Red, Parietale: Orange, Occipitale: Blue)
region_names  = {'Temporal', 'Frontal', 'Parietal', 'Occipital'};
region_colors = [0   1   0;    % Green per Temporal
                 1   0   0;    % Red per Frontal
                 1   0.5 0;    % Orange per Parietal
                 0   0   1];   % Blue per Occipital
             
regionMap = containers.Map(region_names, num2cell(region_colors, 2));

% --- 5) Parametri del modello eterodimero ---
% Questi parametri devono essere gli stessi usati nelle sezioni iniziali del MAIN_Heterodimer
k0 = 1.0;
k1 = 0.5;
ktilde1 = 0.5;
k12 = 0.5;
diffusion_coeff = 5e-4;
dt = 0.4;
num_steps = 100;

% Prealloca memoria per eventuale salvataggio dei risultati
all_times = cell(length(selectedIdx), 1);
all_solutions = cell(length(selectedIdx), 1);

% Mappa per controllare se per una regione è già stata assegnata la legenda
legendSet = containers.Map(region_names, {false, false, false, false});

figure; 
hold on;
title('Heterodimer: Single Node 10% Seeding for Selected 8 Nodes');
xlabel('Time (years)');
ylabel('Average Misfolded Protein Concentration (pt)');

for i = 1:length(selectedIdx)
    idx = selectedIdx(i);
    
    thisLabel = strtrim(char(all_labels{idx}));

    % Determina la regione (Temporal, Frontal, Parietal, Occipital) in base all'etichetta
    if ismember(thisLabel, temp_labels)
        regionCategory = 'Temporal';
    elseif ismember(thisLabel, frontal_labels)
        regionCategory = 'Frontal';
    elseif ismember(thisLabel, parietal_labels)
        regionCategory = 'Parietal';
    elseif ismember(thisLabel, occipital_labels)
        regionCategory = 'Occipital';
    else
        regionCategory = 'Unknown';
    end

    % Recupera il colore della regione
    if isKey(regionMap, regionCategory)
        thisColor = regionMap(regionCategory);
    else
        thisColor = [0 0 0];
    end

    % Se la regione non è ancora stata visualizzata in legenda, assegna DisplayName
    if ~legendSet(regionCategory)
        dispName = regionCategory;       % comparirà in legenda
        handleVis = 'on';
        legendSet(regionCategory) = true; % ora la regione è "già assegnata"
    else
        dispName = '';                   % nessun nome in legenda
        handleVis = 'off';               % e la curva non appare in legenda
    end

    % Esegui la simulazione
    [t_sol, p_sol, pt_sol] = HeterodimerInfection_singleSeed( ...
        A, CoordTable, k0, k1, ktilde1, k12, ...
        diffusion_coeff, dt, num_steps, idx);

    % Calcola l'infezione media (misfolded) e plottala
    avgInfection = mean(pt_sol, 2);
    plot(t_sol, avgInfection, ...
        'Color', thisColor, ...
        'LineWidth', 2, ...         % ispessisce la linea
        'DisplayName', dispName, ...
        'HandleVisibility', handleVis);
end

legend('show','Location','best');
hold off;
