function plotSelectedNodes_FK(A, CoordTable, diffusion, a, dt, num_steps)
    % plotSelectedNodes_FK Plot average infection concentration over time
    % for selected nodes based on their labels in the FK model.
    %
    %   plotSelectedNodes_FK(A, CoordTable, diffusion, a, dt, num_steps)
    %
    % Input arguments:
    %   - A         : Matrice di adiacenza.
    %   - CoordTable: Tabella contenente informazioni sui nodi (la 5a colonna deve
    %                 contenere le etichette dei nodi).
    %   - diffusion : Coefficiente di diffusione (es. 5e-4).
    %   - a         : Parametro di crescita logistica (es. 0.5).
    %   - dt        : Intervallo temporale (es. 0.4 anni).
    %   - num_steps : Numero di step di simulazione (es. 100).
    %
    % La funzione esegue le seguenti operazioni:
    %   1. Definisce le etichette dei nodi per le regioni:
    %         - Temporal: '22L', '22R'
    %         - Frontal : '10L', '10R'
    %         - Parietal: '5L', '5R'
    %         - Occipital: '19L', '19R'
    %
    %   2. Estrae gli indici dei nodi corrispondenti.
    %   3. Definisce una mappatura fissa dei colori per ciascuna regione:
    %         - Temporal: verde, Frontal: rosso, Parietal: arancione, Occipital: blu.
    %   4. Per ciascun nodo selezionato, esegue la simulazione FK con seme singolo
    %      (richiamando la funzione FK_propagation_singleSeed) e plotta la concentrazione
    %      media dell'infezione (media delle concentrazioni lungo i nodi simulati) nel tempo.
    %   5. Visualizza una legenda che associa ad ogni regione il relativo colore.
    
    % 1. Definisci le etichette per ciascuna regione
    temp_labels    = {'22L', '22R'};      % Regione temporale
    frontal_labels = {'10L', '10R'};      % Regione frontale
    parietal_labels= {'5L', '5R'};        % Regione parietale
    occipital_labels = {'19L','19R'};      % Regione occipitale
    selected_labels = [temp_labels, frontal_labels, parietal_labels, occipital_labels];
    
    % 2. Estrai le etichette dei nodi dalla CoordTable (colonna 5)
    all_labels = cellstr(CoordTable{:,5});
    selectedIdx = find(ismember(all_labels, selected_labels));
    
    if numel(selectedIdx) ~= 8
        warning('Si attendevano 8 nodi, ma ne sono stati trovati %d. Continuo comunque.', numel(selectedIdx));
    end
    
    % 3. Definisci i colori per ciascuna regione
    region_names  = {'Temporal', 'Frontal', 'Parietal', 'Occipital'};
    region_colors = [0   1   0;    % Verde per Temporal
                     1   0   0;    % Rosso per Frontal
                     1   0.5 0;    % Arancione per Parietal
                     0   0   1];   % Blu per Occipital
    regionMap = containers.Map(region_names, num2cell(region_colors, 2));
    
    % Inizializza una mappa per controllare la visualizzazione unica in legenda
    legendSet = containers.Map(region_names, {false, false, false, false});
    
    % 4. Prealloca celle per salvare risultati (se necessario)
    all_times = cell(numel(selectedIdx), 1);
    all_solutions = cell(numel(selectedIdx), 1);
    
    % 5. Crea una figura per plottare i risultati
    figure; 
    hold on;
    title('FK: Single Node Seeding per regioni selezionate');
    xlabel('Time (years)');
    ylabel('Average Infection Concentration');
    
    % 6. Loop sui nodi selezionati
    for i = 1:numel(selectedIdx)
        idx = selectedIdx(i);
        thisLabel = strtrim(all_labels{idx});
        
        % Determina la regione in base all'etichetta
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
        
        % Recupera il colore assegnato alla regione
        if isKey(regionMap, regionCategory)
            thisColor = regionMap(regionCategory);
        else
            thisColor = [0 0 0]; % default: nero
        end
        
        % Gestisci la visualizzazione della legenda
        if ~legendSet(regionCategory)
            dispName = regionCategory;
            handleVis = 'on';
            legendSet(regionCategory) = true;
        else
            dispName = '';
            handleVis = 'off';
        end
        
        % Esegui la simulazione per il nodo corrente usando FK_propagation_singleSeed.
        % Si assume che la funzione FK_propagation_singleSeed sia presente nel path.
        [t_sol, c_sol] = FK_propagation_singleSeed(A, CoordTable, diffusion, a, dt, num_steps, idx);
        
        % Calcola la concentrazione media per ogni istante
        avgInfection = mean(c_sol, 2);
        
        % Plotta la curva per il nodo corrente
        plot(t_sol, avgInfection, 'Color', thisColor, 'LineWidth', 2, ...
             'DisplayName', dispName, 'HandleVisibility', handleVis);
         
        % Salva i risultati se necessario
        all_times{i} = t_sol;
        all_solutions{i} = c_sol;
    end
    
    legend('show', 'Location', 'best');
    hold off;
    end
    