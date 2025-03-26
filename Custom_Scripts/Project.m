%% Simulation of Infection Propagation in Human Brain Network
% This script models the spread of an infection within a human brain network
% using a reaction-diffusion model. It leverages the Brain Connectivity Toolbox
% for network analysis and visualization.

clear; close all; clc;
addpath('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/matlab/Custom_Scripts');

%% 1. Load Data
L = load('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/matlab/Laplacian.csv'); % Laplacian
CoordTable = readtable('//Users/Speranza/Desktop/VSCODE/matlab/biomedicina/matlab/Coordinates.csv');   % Coordinates

% Extract numeric coordinates
Coord = table2array(CoordTable(:, 4:6));  % X, Y, Z are in columns 4, 5, 6

% Find the Adjacency matrix from the Laplacian
A = -L;
A(A < 0) = 0;

% Remove small values from the Adjacency matrix (only for graph clarity)
A_tresh = threshold_absolute(A, 0.1);

%% 2. 3D Brain Map
figure;

% ====================== Parameters to Adjust ======================
edgeLineWidth = 0.5;     % Reduce for thinner edges (original was 0.5)
nodeSize = 100;           % Increase for greater node size (original was 50)
% ===================================================================

% Plot edges first
% Define the adjacency_plot_und function if not available
if ~exist('adjacency_plot_und1', 'file')
    % Function not available
    warning('Function adjacency_plot_und not found.');
end


%function [X, Y, Z] = adjacency_plot_und(A, Coord)
%    [i, j] = find(A);
%    X = [Coord(i, 1) Coord(j, 1)]';
%    Y = [Coord(i, 2) Coord(j, 2)]';
%    Z = [Coord(i, 3) Coord(j, 3)]';
%end


[X, Y, Z] = adjacency_plot_und1(A_tresh, Coord);
plot3(X, Y, Z, 'y-', 'LineWidth', edgeLineWidth); 
hold on;

% Prepare node colors based on brain regions
brainRegions = unique(CoordTable{:, 7});  % Get unique region names
colors = lines(length(brainRegions));     % Color palette

% Initialize array to store scatter objects for legend
hScatter = gobjects(length(brainRegions), 1); 

% Plot nodes with different colors
for i = 1:length(brainRegions)
    region = brainRegions{i};
    isRegion = strcmp(CoordTable{:, 7}, region);
    
    % Store scatter handle for legend
    hScatter(i) = scatter3(Coord(isRegion,1), Coord(isRegion,2), Coord(isRegion,3),...
        nodeSize, colors(i,:), 'filled', 'DisplayName', region);
end

% Finalize plot
title('3D Brain Network');
xlabel('X'); ylabel('Y'); zlabel('Z');
legend(hScatter, 'Location', 'eastoutside'); % Use SCATTER handles for legend
grid on; 
axis equal;
hold off;

%% 3. Infection Model

%%%%%%%%%%% Initial Conditions %%%%%%%%%%
N = size(A, 1);
c0 = zeros(N, 1);

% Identify and seed ONLY the selected brain regions
infected_mask = strcmp(CoordTable{:, 3}, 'entorhinal')
                % strcmp(CoordTable{:, 3}, 'temporalpole'); % as said in the paper, nel paper non vedo dove parli anche di temporalphole

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
f = @(t, c) -L*c*0.05 + a*c.*(1-c);

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

%% 5. Plot Infection Dynamics
% Visualizzazione spaziale della distribuzione iniziale e finale dell'infezione
% con una colormap personalizzata (da grigio chiaro a rosso bordeaux) e
% evidenziazione dei nodi inizialmente infetti tramite un cerchio sottile attorno.

figure;

% Definisci i limiti comuni dei colori (0 a 1)
common_clims = [0 1];

% Definizione della colormap personalizzata: da grigio chiaro (0) a rosso bordeaux (1)
nColors = 256;
startColor = [0.9, 0.9, 0.9];   % Grigio chiaro per valore 0
endColor   = [0.5, 0, 0.13];      % Rosso bordeaux per valore 1
customCmap = [linspace(startColor(1), endColor(1), nColors)', ...
              linspace(startColor(2), endColor(2), nColors)', ...
              linspace(startColor(3), endColor(3), nColors)'];

% Plot del primo subplot: distribuzione iniziale
subplot(1,2,1);
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, c0, 'filled');
title('Initial Infection Distribution');
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on;
clim(common_clims);
colormap(customCmap);   % Applica la colormap personalizzata

% Aggiungi un cerchio sottile attorno ai nodi inizialmente infetti (c0 > 0)
infectedIdx = find(c0 > 0);  % indice dei nodi infetti
circleRadius = 2;  % Puoi regolare questo valore in base alla scala delle coordinate
theta = linspace(0, 2*pi, 50);  % Risoluzione del cerchio
hold on;
for i = 1:length(infectedIdx)
    idx = infectedIdx(i);
    x_center = Coord(idx, 1);
    y_center = Coord(idx, 2);
    z_center = Coord(idx, 3);
    % Calcola le coordinate del cerchio nel piano XY, mantenendo lo stesso z
    x_circle = x_center + circleRadius * cos(theta);
    y_circle = y_center + circleRadius * sin(theta);
    z_circle = repmat(z_center, size(theta));
    % Disegna il cerchio in nero con linea sottile
    plot3(x_circle, y_circle, z_circle, 'k-', 'LineWidth', 1);
end
hold off;

% Plot del secondo subplot: distribuzione finale
subplot(1,2,2);
scatter3(Coord(:,1), Coord(:,2), Coord(:,3), 50, c_sol(end, :), 'filled');
title(['Final Infection Distribution at t = ' num2str(tmax)]);
xlabel('X'); ylabel('Y'); zlabel('Z');
axis equal; grid on;
clim(common_clims);
colormap(customCmap);   % Applica la stessa colormap

% Aggiungi una colorbar unica per entrambi i subplot
hcb = colorbar('Position', [0.92 0.15 0.02 0.7]);
ylabel(hcb, 'Infection Concentration');

%% 5.5 Change nodes used (use only 10 central nodes per area)
%selezioniamo i nodi in base alla centralità data dalla matrice L
% vogliamo per ogni sottomatrice definita in base alla tipologia di nodo, come ad esempio frontale, occipitale, etc i primi 10 nodi più centrali
% per fare ciò dobbiamo calcolare la centralità di ogni nodo e poi selezionare i primi 10 nodi più centrali

group_temporal = {'inferiortemporal', 'middletemporal', 'Right-Hippocampus', 'Left-Hippocampus','fusiform','temporalpole','superiortemporal','transversetemporal','bankssts'};
group_frontal  = {'precentral', 'caudalmiddlefrontal','lateralorbitofrontal', 'rostralmiddlefrontal', 'superiorfrontal', 'parsopercularis', 'parstriangularis', 'parsorbitalis', 'rostralanteriorcingulate', 'caudalanteriorcingulate', 'frontalpole','medialorbitofrontal'};
group_parietal = {'superiorparietal', 'supramarginal', 'precuneus','inferiorparietal','postcentral','paracentral','isthmuscingulate'};
group_occipital = {'lateraloccipital', 'cuneus', 'pericalcarine', 'lingual'};
region_groups = {group_temporal, group_frontal, group_parietal, group_occipital};
region_names = {'Temporal', 'Frontal', 'Parietal', 'Occipital'};
num_groups = length(region_names);

% Calcolare la centralità di grado per tutti i nodi
degree_centrality = sum(A, 2);

% Inizializzare una struttura per memorizzare i nodi centrali per ogni regione
top_central_nodes = struct();

% Per ciascun gruppo, selezionare i 10 nodi più centrali
for r = 1:num_groups
    group = region_groups{r};
    nodes_in_region = ismember(lower(CoordTable{:,3}), lower(group));
    centrality_in_region = degree_centrality(nodes_in_region);
    region_indices = find(nodes_in_region);
    [~, sorted_indices] = sort(centrality_in_region, 'descend');
    top_10_indices = region_indices(sorted_indices(1:min(10, length(sorted_indices))));
    top_central_nodes.(region_names{r}) = top_10_indices;
end

% Visualizzare i nodi centrali per ogni area
for r = 1:num_groups
    fprintf('Top 10 central nodes in %s region:\n', region_names{r});
    disp(CoordTable(top_central_nodes.(region_names{r}), [1 3 7]));
end

% Creare un vettore unificato di tutti i nodi centrali
central_nodes = unique([top_central_nodes.Temporal; top_central_nodes.Frontal; top_central_nodes.Parietal; top_central_nodes.Occipital]);

%% 6. Plot Average Infection Concentration Over Time (Central Nodes Only)
% Calcola la concentrazione media dell'infezione considerando solo i nodi centrali
avg_conc_central = mean(c_sol(:, central_nodes), 2);

figure;
plot(t_sol, avg_conc_central, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); % Blu default MATLAB
xlabel('Time');
ylabel('Average Infection Concentration (Central Nodes)');
title('Average Infection Concentration Over Time (Central Nodes Only)');
grid on;
hold off;

%% 7. Compute Average Infection Concentration for Selected Brain Regions Over Time (Central Nodes Only)
% Per ciascuna regione, calcolare la concentrazione media usando solo i nodi centrali selezionati

avg_conc_desired_central = zeros(length(t_sol), num_groups);
for r = 1:num_groups
    nodes_in_region = top_central_nodes.(region_names{r});
    if isempty(nodes_in_region)
        warning('No central nodes found for region %s. Setting average to 0.', region_names{r});
        avg_conc_desired_central(:, r) = 0;
    else
        avg_conc_desired_central(:, r) = mean(c_sol(:, nodes_in_region), 2);
    end
end

figure('Name', 'Average Infection Concentration for Central Nodes per Brain Region Over Time', 'NumberTitle', 'off');
hold on;
for r = 1:num_groups
    % Nota: region_colors deve essere definito in anticipo; se non lo è, è possibile definirlo qui:
    region_colors = [0 1 0; 1 0 0; 1 0.5 0; 0 0 1];
    plot(t_sol, avg_conc_desired_central(:, r), 'LineWidth', 2, 'Color', region_colors(r, :));
end
hold off;
xlabel('Time');
ylabel('Average Infection Concentration');
title('Average Infection Concentration per Selected Brain Region (Central Nodes) Over Time');
legend(region_names, 'Location', 'bestoutside');
grid on;

%% 8. Create Video of Infection Propagation on 3D Brain Map (Central Nodes Only)
% Visualizzazione della propagazione dell'infezione limitata ai soli nodi centrali

% Estrai la sottomatrice della matrice di adiacenza per i nodi centrali
A_central = A_tresh(central_nodes, central_nodes);
% Ottieni le coordinate per il plot degli spigoli utilizzando i nodi centrali
[X_c, Y_c, Z_c] = adjacency_plot_und1(A_central, Coord(central_nodes,:));

figure('Color','w', 'Position', [100, 100, 1120, 840]); 

% Imposta il VideoWriter per salvare la propagazione
videoFilename = 'InfectionPropagation_CentralNodes.avi';
v = VideoWriter(videoFilename);
v.FrameRate = 5;  % Regola il frame rate secondo necessità
open(v);

for i = 1:length(t_sol)
    clf;
    hold on;
    % Plot degli spigoli per i nodi centrali
    plot3(X_c, Y_c, Z_c, 'y-', 'LineWidth', edgeLineWidth);
    
    % Seleziona il livello di infezione corrente per i nodi centrali
    infectionLevel_central = c_sol(i, central_nodes);
    
    % Plot dei nodi centrali
    scatter3(Coord(central_nodes,1), Coord(central_nodes,2), Coord(central_nodes,3), ...
             nodeSize, infectionLevel_central, 'filled');
    
    colormap(jet);
    clim([0 1]);
    colorbar;
    
    title(sprintf('Infection Propagation (Central Nodes) at t = %.2f', t_sol(i)));
    xlabel('X'); ylabel('Y'); zlabel('Z');
    axis equal; grid on;
    view(3);
    
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
    
    hold off;
end

close(v);
disp(['Video saved as ', videoFilename]);

%% Modello eterodimero

%% 9 da qui stoc.
% Simulazione con variazione dinamica degli archi su A (40 anni)
% Impostazioni iniziali
L = load('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/matlab/Laplacian.csv'); % Laplacian

dt = 0.4;             % passo temporale in anni
num_steps = 100;      % numero di passi (T = 0.4*100 = 40 anni)
t_sol = linspace(0, dt*num_steps, num_steps+1);
N = size(A, 1);       % numero di nodi
c_sol = zeros(num_steps+1, N);

% Condizioni iniziali: ad esempio, se "entorhinal" è la regione di interesse
c0 = zeros(N, 1);
infected_mask = strcmp(CoordTable{:, 3}, 'entorhinal');
c0(infected_mask) = 0.1;
c_sol(1, :) = c0';

% Parametri del modello
a = 0.5;                % parametro di crescita logistica
edge_reduction = 0.8;   % fattore di riduzione degli archi (es. 0.8 riduce del 20%)
diffusion_coeff = 0.05; % coefficiente di diffusione

for i = 2:length(t_sol)
    % --- Aggiornamento casuale degli archi su A ---
    % Estrai gli indici degli archi esistenti (considera solo la parte superiore per simmetria)
    [edgeI, edgeJ] = find(triu(A, 1));
    numEdgesTotal = length(edgeI);
    
    % Definisci il numero massimo di archi da modificare (ad es. il 10% degli archi esistenti)
    maxEdgesToModify = ceil(0.15 * numEdgesTotal);
    % Scegli casualmente un numero (minimo 1) di archi da modificare
    numEdgesToModify = randi([1, maxEdgesToModify]);
    numEdgesToModify
    
    % Seleziona casualmente gli archi da modificare
    selectedIndices = randperm(numEdgesTotal, numEdgesToModify);
    for k = 1:numEdgesToModify
        idx1 = edgeI(selectedIndices(k));
        idx2 = edgeJ(selectedIndices(k));
        % Riduci il peso dell'arco selezionato (ricorda la simmetria)
        A(idx1, idx2) = A(idx1, idx2) * edge_reduction;
        A(idx2, idx1) = A(idx2, idx1) * edge_reduction;
    end
    
    % --- Aggiornamento della Laplaciana ---
    % Calcola la Laplaciana corrente basata su A aggiornato
    L_current = diag(sum(A, 2)) - A;
    
    % --- Aggiornamento del modello Reaction-Diffusion ---
    % dc/dt = - diffusion_coeff * L_current * c + a * c .* (1 - c)
    temp = c_sol(i-1, :)';         % vettore colonna della concentrazione al passo precedente
    diffusion = -diffusion_coeff * L_current * temp;
    growth = a * temp .* (1 - temp);
    temp = temp + dt * (diffusion + growth);  % aggiornamento con passo implicito
    c_sol(i, :) = temp';           % salva la concentrazione al passo corrente
end

% A questo punto c_sol contiene l'evoluzione temporale della concentrazione 
% d'infezione per ogni nodo, mentre A è stata modificata ad ogni iterazione.

%% 10. plot della propagazione dell'infezione con variazione dinamica degli archi
% Plot Average Infection Concentration Over Time (Central Nodes Only)
% Calcola la concentrazione media dell'infezione considerando solo i nodi centrali
avg_conc_central = mean(c_sol(:, central_nodes), 2);

figure;
plot(t_sol, avg_conc_central, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); % Blu default MATLAB
xlabel('Time');
ylabel('Average Infection Concentration (Central Nodes)');
title('Average Infection Concentration Over Time (Central Nodes Only)');
grid on;
hold off;
