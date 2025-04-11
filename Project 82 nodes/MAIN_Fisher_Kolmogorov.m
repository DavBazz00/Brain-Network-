clear; close all; clc;
addpath('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/project_82_withfunction/functions')

%% 1. Load file

edgeFile = '/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Edge.csv';
nodeFile = '/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Node.csv';

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


%% 9. Plot average infection with starting point different nodes (FK)

% --- 2) Definisci le etichette selezionate per ciascuna regione ---
temp_labels    = {'22L', '22R'};   % Temporal
frontal_labels = {'10L', '10R'};   % Frontal
parietal_labels= {'5L',  '5R'};     % Parietal
occipital_labels = {'19L','19R'};   % Occipital

% Combina tutte le etichette (8 nodi totali)
selected_labels = [temp_labels, frontal_labels, parietal_labels, occipital_labels];

% --- 3) Estrai le etichette dei nodi dalla CoordTable ---
% (Si assume che la quinta colonna di CoordTable contenga le etichette dei nodi)
all_labels = cellstr(CoordTable{:,5});
selectedIdx = find(ismember(all_labels, selected_labels));

if length(selectedIdx) ~= 8
    error('Expected to find 8 nodes with the specified labels, but found %d.', length(selectedIdx));
end

% --- 4) Definisci la mappatura fissa dei colori per le 4 regioni ---
region_names  = {'Temporal', 'Frontal', 'Parietal', 'Occipital'};
region_colors = [0   1   0;    % Green per Temporal
                 1   0   0;    % Red per Frontal
                 1   0.5 0;    % Orange per Parietal
                 0   0   1];   % Blue per Occipital
regionMap = containers.Map(region_names, num2cell(region_colors, 2));

% --- Crea una mappa per controllare la visualizzazione unica in legenda ---
legendSet = containers.Map(region_names, {false, false, false, false});

% --- 5) Parametri del modello FK ---
diffusion = 5e-4; 
a         = 0.5;
dt        = 0.4;     
num_steps = 100;

% Prealloca se necessario
all_times = cell(length(selectedIdx), 1);
all_solutions = cell(length(selectedIdx), 1);

% --- 6) Crea una figura per plottare i risultati ---
figure; 
hold on;
title('FK: Single Node 10% Seeding for Selected 8 Nodes');
xlabel('Time (years)');
ylabel('Average Infection');

% --- 7) Loop sui nodi selezionati ---
for i = 1:length(selectedIdx)
    idx = selectedIdx(i);
    
    % Recupera l'etichetta del nodo (pulita)
    thisLabel = strtrim(char(all_labels{idx}));
    
    % Determina la regione in base alle etichette fornite
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
    
    % Recupera il colore fissato per la regione
    if isKey(regionMap, regionCategory)
        thisColor = regionMap(regionCategory);
    else
        thisColor = [0 0 0]; % default nero
    end
    
    % Se per questa regione non è ancora stata assegnata la voce in legenda, impostala;
    % altrimenti, disabilita la visualizzazione del handle.
    if ~legendSet(regionCategory)
        dispName = regionCategory;
        handleVis = 'on';
        legendSet(regionCategory) = true;
    else
        dispName = '';
        handleVis = 'off';
    end
    
    % Esegui la simulazione per il nodo corrente usando la funzione FK_propagation_singleSeed
    [t_sol, c_sol] = FK_propagation_singleSeed(A, CoordTable, diffusion, a, dt, num_steps, idx);
    
    % Calcola l'infezione media (media delle concentrazioni)
    avgInfection = mean(c_sol, 2);
    
    % Plotta la curva con linea ispessita (LineWidth=2) e con DisplayName e HandleVisibility impostati
    plot(t_sol, avgInfection, 'Color', thisColor, 'LineWidth', 2, ...
         'DisplayName', dispName, 'HandleVisibility', handleVis);
     
    % Salva, se necessario, i risultati
    all_times{i} = t_sol;
    all_solutions{i} = c_sol;
end

legend('show', 'Location', 'best');
hold off;

