function plotSelectedNodes_Heterodimer(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps)
    % plotSelectedNodes_Heterodimer Plot average misfolded protein concentration
    % over time (modelo eterodimero) for selected nodes based on their labels.
    %
    %   plotSelectedNodes_Heterodimer(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps)
    %
    % Input arguments:
    %   - A             : Matrice di adiacenza.
    %   - CoordTable    : Tabella contenente le informazioni sui nodi (la 5a colonna
    %                     deve contenere le etichette dei nodi).
    %   - k0            : Tasso di produzione della proteina sana.
    %   - k1            : Tasso di clearance della proteina sana.
    %   - ktilde1       : Tasso di clearance della proteina misfolded.
    %   - k12           : Tasso di conversione da proteina sana a misfolded.
    %   - diffusion_coeff: Coefficiente di diffusione (e.g., 5e-4).
    %   - dt            : Intervallo temporale (in anni).
    %   - num_steps     : Numero di passi temporali.
    %
    % La funzione esegue le seguenti operazioni:
    %   1. Definisce le etichette dei nodi selezionati per le regioni:
    %         - Temporale: '22L', '22R'
    %         - Frontale : '10L', '10R'
    %         - Parietale: '5L',  '5R'
    %         - Occipitale: '19L', '19R'
    %   2. Estrae gli indici dei nodi corrispondenti dalla colonna 5 di CoordTable.
    %   3. Definisce una mappatura fissa dei colori per ciascuna regione:
    %         - Temporale: verde, Frontale: rosso, Parietale: arancione, Occipitale: blu.
    %   4. Per ciascun nodo selezionato, esegue la simulazione del modello eterodimero
    %      a singolo seme (richiamando la funzione HeterodimerInfection_singleSeed) e plotta
    %      la concentrazione media della proteina misfolded (pt) nel tempo.
    %   5. Visualizza una legenda che associa ad ogni regione il relativo colore.
    %
    % Assicurati che la funzione HeterodimerInfection_singleSeed sia nel path.
    
    % 1. Definisci le etichette per le regioni
    temp_labels    = {'22L', '22R'};      % Regione Temporale
    frontal_labels = {'10L', '10R'};      % Regione Frontale
    parietal_labels= {'5L',  '5R'};        % Regione Parietale
    occipital_labels = {'19L','19R'};      % Regione Occipitale
    selected_labels = [temp_labels, frontal_labels, parietal_labels, occipital_labels];
    
    % 2. Estrai gli indici dei nodi dalla quinta colonna di CoordTable
    all_labels = cellstr(CoordTable{:,5});
    selectedIdx = find(ismember(all_labels, selected_labels));
    
    if numel(selectedIdx) ~= 8
        error('Si attendevano 8 nodi con le etichette specificate, ma sono stati trovati %d.', numel(selectedIdx));
    end
    
    % 3. Definisci i colori per le regioni
    region_names  = {'Temporal', 'Frontal', 'Parietal', 'Occipital'};
    region_colors = [0   1   0;    % Verde per Temporale
                     1   0   0;    % Rosso per Frontale
                     1   0.5 0;    % Arancione per Parietale
                     0   0   1];   % Blu per Occipitale
    regionMap = containers.Map(region_names, num2cell(region_colors, 2));
    
    % Crea una mappa per gestire la legenda (visualizza un'unica voce per regione)
    legendSet = containers.Map(region_names, {false, false, false, false});
    
    % 4. Prealloca, se necessario, le celle per salvare i risultati
    all_times = cell(numel(selectedIdx), 1);
    all_solutions = cell(numel(selectedIdx), 1);
    
    % 5. Crea la figura per il plotting
    figure; 
    hold on;
    title('Heterodimer: Single Node Seeding per Regioni Selezionate');
    xlabel('Time (years)');
    ylabel('Average Misfolded Protein Concentration (pt)');
    
    % 6. Loop sui nodi selezionati
    for i = 1:numel(selectedIdx)
        idx = selectedIdx(i);
        thisLabel = strtrim(all_labels{idx});
        
        % Determina la regione basata sull'etichetta
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
        
        % Recupera il colore corrispondente alla regione
        if isKey(regionMap, regionCategory)
            thisColor = regionMap(regionCategory);
        else
            thisColor = [0 0 0]; % default nero
        end
        
        % Gestisci la legenda: visualizza il nome della regione solo la prima volta che appare
        if ~legendSet(regionCategory)
            dispName = regionCategory;
            handleVis = 'on';
            legendSet(regionCategory) = true;
        else
            dispName = '';
            handleVis = 'off';
        end
        
        % Esegui la simulazione per il nodo corrente tramite HeterodimerInfection_singleSeed
        [t_sol, p_sol, pt_sol] = HeterodimerInfection_singleSeed(...
            A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps, idx);
        
        % Calcola la concentrazione media della proteina misfolded (pt) lungo il tempo
        avgInfection = mean(pt_sol, 2);
        
        % Plotta la curva con specifiche impostazioni di colore e legenda
        plot(t_sol, avgInfection, 'Color', thisColor, 'LineWidth', 2, ...
             'DisplayName', dispName, 'HandleVisibility', handleVis);
         
        % Salva i risultati se necessario
        all_times{i} = t_sol;
        all_solutions{i} = pt_sol;
    end
    
    legend('show','Location','best');
    hold off;
    end
    