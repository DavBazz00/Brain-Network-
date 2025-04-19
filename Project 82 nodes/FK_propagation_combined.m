function [t, c] = FK_propagation_combined(A, CoordTable, diffusion, a, dt_aging, dt_treat, T_end, t_switch, aging_reduction, cure_reduction)
    % FK_PROPAGATION_COMBINED Simula invecchiamento dinamico e trattamento annuale
    %
    % INPUTS:
    %   A               - matrice di adiacenza iniziale
    %   CoordTable      - tabella con informazioni sui nodi
    %   diffusion       - coefficiente di diffusione
    %   a               - parametro di crescita logistica
    %   dt_aging        - intervallo di invecchiamento (anni)
    %   dt_treat        - intervallo di applicazione del trattamento (anni)
    %   T_end           - durata totale simulazione (anni)
    %   t_switch        - anno di inizio trattamento
    %   aging_reduction - fattore moltiplicativo per aging (0-1)
    %   cure_reduction  - frazione di riduzione globale per cura (0-1)
    %
    % OUTPUTS:
    %   t - vettore tempi da 0 a T_end con passo dt_aging
    %   c - matrice concentrazioni (length(t) x N)
    
    % Costruzione vettore tempo
    t = 0:dt_aging:T_end;
    num_steps = length(t);
    N = size(A,1);
    
    % Preallocazione conservazione concentrazioni
    c = zeros(num_steps, N);
    
    % Condizione iniziale: semina entorinale
    c0 = zeros(N,1);
    infected_mask = strcmp(CoordTable{:,4}, 'Entorhinal');
    c0(infected_mask) = 0.1;
    c(1,:) = c0';
    
    % Copia matrice di adiacenza per aggiornamenti
    A_curr = A;
    
    % Inizializza tempo del primo trattamento
    next_cure = t_switch;
    
    % Funzione di reazione-diffusione
    dyn = @(L, cvec) -L * cvec * diffusion + a * cvec .* (1 - cvec);
    
    % Loop principale
    for i = 2:num_steps
        t_now = t(i);
    
        % --- Aging dinamico (ogni dt_aging) ---
        [ei, ej] = find(triu(A_curr,1));
        M = numel(ei);
        Kmax = ceil(0.8 * M);
        if Kmax >= 2
            if mod(Kmax,2) == 1
                Kmax = Kmax - 1;
            end
            K = 2 * randi([1, Kmax/2]);
            sel = randperm(M, K);
            for k = 1:2:K
                u = ei(sel(k));
                v = ej(sel(k+1));
                A_curr(u,v) = A_curr(u,v) * aging_reduction;
                A_curr(v,u) = A_curr(v,u) * aging_reduction;
            end
        end
    
        % --- Trattamento annuale (ogni dt_treat anni), dopo t_switch ---
        if t_now >= next_cure
            A_curr = A_curr * cure_reduction;
            next_cure = next_cure + dt_treat;
        end
    
        % --- Integrazione forward Euler ---
        L = diag(sum(A_curr,2)) - A_curr;
        prev = c(i-1,:)';
        cnew = prev + dt_aging * dyn(L, prev);
        c(i,:) = cnew';
    end
    end
    