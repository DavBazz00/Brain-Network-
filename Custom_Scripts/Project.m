clear; close all; clc;

%% 1.Load Data
A = load('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Edge.csv'); % Adjacency
A(1, :) = []; % Remove the first row

CoordTable = readtable('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Node.csv'); % Coordinates
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
f = @(t, c) -L*c*8e-4 + a*c.*(1-c); 

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


%% 7. infection by brain region
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

%% Heterodimer Infection Model (Adattato)

% Parametri cinetici del modello eterodimero
k0 = 1.0;       % Tasso di produzione della proteina sana
k1 = 0.5;       % Tasso di clearance della proteina sana
ktilde1 = 0.5;  % Tasso di clearance della proteina misfolddata
k12 = 0.5;      % Tasso di conversione da proteina sana a misfolddata
diffusion_coeff = 0.05;  % Coefficiente di diffusione

%%%%%%%%%%% Condizioni Iniziali %%%%%%%%%%
N = size(A, 1);
% Concentrazione iniziale della proteina sana: p0 = k0/k1 (stato sano)
p0 = ones(N, 1) * (k0/k1);

% Concentrazione iniziale della proteina misfolddata: inizialmente zero
pt0 = zeros(N, 1);
% Seeding: imposta un seme (0.1) nella regione "Entorhinal" (stessa colonna usata nel modello FK)
infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
pt0(infected_mask) = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prealloca le matrici soluzione
p_sol = zeros(num_steps + 1, N);    % Concentrazione di proteina sana nel tempo
pt_sol = zeros(num_steps + 1, N);   % Concentrazione di proteina misfolddata nel tempo
p_sol(1, :) = p0';
pt_sol(1, :) = pt0';

disp('Inizio simulazione del modello eterodimero...');
for i = 2:length(t_sol)
    % Estrai le concentrazioni dal passo precedente
    p_old = p_sol(i-1, :)';      % Proteina sana
    pt_old = pt_sol(i-1, :)';    % Proteina misfolddata

    % Calcola la derivata combinando diffusione e reazione
    dpdt = - diffusion_coeff * L * p_old + k0 - k1 * p_old - k12 * (p_old .* pt_old);
    dptdt = - diffusion_coeff * L * pt_old - ktilde1 * pt_old + k12 * (p_old .* pt_old);
    
    % Aggiornamento esplicito (Forward Euler)
    p_new = p_old + dt * dpdt;
    pt_new = pt_old + dt * dptdt;
    
    % Salva le nuove concentrazioni
    p_sol(i, :) = p_new';
    pt_sol(i, :) = pt_new';
end

disp('Simulazione eterodimero completata con successo.');

% Calcola il biomarcatore totale (somma su tutti i nodi)
biomarker_pt = sum(pt_sol(end, :));
fprintf('Biomarker (totale proteina misfolddata) al tempo finale: %.3f\n', biomarker_pt);

% Calcola la concentrazione media su TUTTI i nodi
avg_pt_all = mean(pt_sol, 2);
figure;
plot(t_sol, avg_pt_all, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
xlabel('Time [years]');
ylabel('Average Misfolded Protein Concentration (All Nodes)');
title('Heterodimer Model: p_{tilde} Evolution (All Nodes)');
grid on;

%% 9. Simulazione con variazione dinamica degli archi (FK)
%ricarica L
A = load('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Edge.csv'); % Adjacency
A(1, :) = []; % Remove the first row

CoordTable = readtable('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Node.csv'); % Coordinates
% Extract numeric coordinates
Coord = table2array(CoordTable(:, 1:3));

% Degree matrix
D = diag(sum(A, 2));
% Laplacian matrix
L = D - A;

dt = 0.4;             % passo temporale (anni)
num_steps = 100;      % numero di passi (T = 40 anni)
t_sol = linspace(0, dt*num_steps, num_steps+1);
N = size(A, 1);
c_sol = zeros(num_steps+1, N);

% Condizioni iniziali: seeding nella regione "Entorhinal" (stessa come nel modello FK)
c0 = zeros(N, 1);
infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
c0(infected_mask) = 0.1;
c_sol(1, :) = c0';

% Parametri del modello
a = 0.5;                % parametro di crescita logistica
edge_reduction = 0.75;  % fattore di riduzione degli archi
diffusion_coeff = 0.05; % coefficiente di diffusione

for i = 2:length(t_sol)
    % --- Aggiornamento casuale degli archi su A ---
    [edgeI, edgeJ] = find(triu(A, 1));
    numEdgesTotal = length(edgeI);
    
    maxEdgesToModify = ceil(0.2 * numEdgesTotal);
    numEdgesToModify = randi([1, maxEdgesToModify]);
    
    selectedIndices = randperm(numEdgesTotal, numEdgesToModify);
    for k = 1:numEdgesToModify
        idx1 = edgeI(selectedIndices(k));
        idx2 = edgeJ(selectedIndices(k));
        A(idx1, idx2) = A(idx1, idx2) * edge_reduction;
        A(idx2, idx1) = A(idx2, idx1) * edge_reduction;
    end
    
    % --- Aggiornamento della Laplaciana ---
    L_current = diag(sum(A, 2)) - A;
    
    % --- Aggiornamento del modello Reaction-Diffusion ---
    temp = c_sol(i-1, :)';
    diffusion = -diffusion_coeff * L_current * temp * 5e-4;
    growth = a * temp .* (1 - temp);
    temp = temp + dt * (diffusion + growth);
    c_sol(i, :) = temp';
end

% Calcola la concentrazione media su TUTTI i nodi (anziché solo i nodi centrali)
avg_conc_all = mean(c_sol, 2);
figure;
plot(t_sol, avg_conc_all, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
xlabel('Time');
ylabel('Average Infection Concentration (All Nodes)');
title('Average Infection Concentration Over Time (All Nodes)');
grid on;

% d'infezione per ogni nodo, mentre A è stata modificata ad ogni iterazione.

%% 10. Modello Eterodimero con Variazione Dinamica degli Archi

% Ricarica la matrice di adiacenza e ricalcola la Laplaciana
A_dynamic = load('/Users/Speranza/Desktop/VSCODE/matlab/biomedicina/82_nodes/Edge.csv');
A_dynamic(1, :) = [];  % Rimuovi la prima riga
D_dynamic = diag(sum(A_dynamic, 2));
L_dynamic = D_dynamic - A_dynamic;

% Parametri cinetici del modello eterodimero (già definiti in precedenza)
% k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps, t_sol

% Reinizializza le condizioni iniziali per il modello eterodimero
N = size(A_dynamic, 1);
p0 = ones(N, 1) * (k0/k1);  % Concentrazione di proteina sana (stato sano)
pt0 = zeros(N, 1);          % Concentrazione di proteina misfolddata
% Seeding nella regione "Entorhinal" (stessa colonna usata nel modello FK)
infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
pt0(infected_mask) = 0.1;

% Prealloca le matrici soluzione per il modello eterodimero dinamico
p_sol_dyn = zeros(num_steps + 1, N);
pt_sol_dyn = zeros(num_steps + 1, N);
p_sol_dyn(1, :) = p0';
pt_sol_dyn(1, :) = pt0';

% Parametri per l'aggiornamento dinamico degli archi
edge_reduction_dyn = 0.5;  % Fattore di riduzione (es. 0.80 per ridurre dell'80%)

disp('Avvio simulazione del modello eterodimero con variazione dinamica degli archi...');
for i = 2:length(t_sol)
    % --- Aggiornamento casuale degli archi su A_dynamic ---
    [edgeI, edgeJ] = find(triu(A_dynamic, 1));
    numEdgesTotal = length(edgeI);
    
    maxEdgesToModify = ceil(1 * numEdgesTotal);
    numEdgesToModify = randi([1, maxEdgesToModify]);
    
    selectedIndices = randperm(numEdgesTotal, numEdgesToModify);
    for k = 1:numEdgesToModify
        idx1 = edgeI(selectedIndices(k));
        idx2 = edgeJ(selectedIndices(k));
        A_dynamic(idx1, idx2) = A_dynamic(idx1, idx2) * edge_reduction_dyn;
        A_dynamic(idx2, idx1) = A_dynamic(idx2, idx1) * edge_reduction_dyn;
    end
    
    % --- Aggiornamento della Laplaciana corrente ---
    L_current_dyn = diag(sum(A_dynamic, 2)) - A_dynamic;
    
    % --- Aggiornamento del modello eterodimero con L dinamica ---
    % Estrai le concentrazioni dal passo precedente
    p_old = p_sol_dyn(i-1, :)';
    pt_old = pt_sol_dyn(i-1, :)';
    
    % Calcola le derivate combinando diffusione (con L_current_dyn) e reazione
    dpdt_dyn = - diffusion_coeff * L_current_dyn * p_old + k0 - k1 * p_old - k12 * (p_old .* pt_old);
    dptdt_dyn = - diffusion_coeff * L_current_dyn * pt_old - ktilde1 * pt_old + k12 * (p_old .* pt_old);
    
    % Aggiornamento esplicito (Forward Euler)
    p_new_dyn = p_old + dt * dpdt_dyn;
    pt_new_dyn = pt_old + dt * dptdt_dyn;
    
    % Salva le nuove concentrazioni
    p_sol_dyn(i, :) = p_new_dyn';
    pt_sol_dyn(i, :) = pt_new_dyn';
end

disp('Simulazione eterodimero dinamica completata.');

% Calcola e visualizza la concentrazione media di proteina misfolddata su tutti i nodi
avg_pt_all_dyn = mean(pt_sol_dyn, 2);
figure;
plot(t_sol, avg_pt_all_dyn, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
xlabel('Time [years]');
ylabel('Average Misfolded Protein Concentration (All Nodes)');
title('Heterodimer Model with Dynamic Edge Variation: p_{tilde} Evolution');
grid on;

%% 11 grafico per singolo nodo del modello FK
% Grafico dell'evoluzione temporale dell'infezione per ciascun nodo (Modello FK)
figure;
hold on;
% Prealloca un vettore di handle per ciascuna regione
h = gobjects(num_groups,1);
for i = 1:num_groups
    region = region_names{i};
    idx = strcmp(all_regions, region);
    region_nodes = find(idx);
    for j = 1:length(region_nodes)
        h_temp = plot(t_sol, c_sol(:, region_nodes(j)), 'Color', region_colors(i, :), 'LineWidth', 1);
        % Salva il handle della prima linea per la regione i
        if j == 1
            h(i) = h_temp;
        end
    end
end
xlabel('Time');
ylabel('Infection Concentration (FK Model)');
title('Evoluzione dell''Infezione per Nodo (Modello FK)');
legend(h, region_names, 'Location', 'northeast');
grid on;
hold off;

%% 12. Grafico per singolo nodo del modello eterodimero
% Grafico dell'evoluzione temporale della proteina misfolddata per ciascun nodo (Modello Eterodimero)
figure;
hold on;
% Prealloca un vettore di handle per ciascuna regione
h = gobjects(num_groups,1);
for i = 1:num_groups
    region = region_names{i};
    idx = strcmp(all_regions, region);
    region_nodes = find(idx);
    for j = 1:length(region_nodes)
        h_temp = plot(t_sol, pt_sol(:, region_nodes(j)), 'Color', region_colors(i, :), 'LineWidth', 1);
        % Salva il handle della prima linea per la regione i
        if j == 1
            h(i) = h_temp;
        end
    end
end
xlabel('Time');
ylabel('Misfolded Protein Concentration (Heterodimer Model)');
title('Evoluzione dell''Infezione per Nodo (Modello Eterodimero)');
legend(h, region_names, 'Location', 'northeast');
grid on;
hold off;

