function [t_sol, p_sol, pt_sol] = HeterodimerInfection_singleSeed(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps, seedNode)
    % HETERODIMERINFECTION_SINGLESEED simula il modello eterodimero con semina
    % iniziale in uno (o in più) nodi specificati.
    %
    % INPUT:
    %   A              - Matrice di adiacenza (NxN).
    %   CoordTable     - Tabella dei nodi (si assume che la colonna 4 contenga le etichette).
    %   k0             - Tasso di produzione delle proteine sane.
    %   k1             - Tasso di clearance delle proteine sane.
    %   ktilde1       - Tasso di clearance delle proteine misfolded.
    %   k12            - Tasso di conversione da proteina sana a misfolded.
    %   diffusion_coeff- Coefficiente di diffusione.
    %   dt             - Time step (in anni).
    %   num_steps      - Numero di passi temporali.
    %   seedNode       - Indice (o vettore di indici) del nodo da seminare (pt=0.1).
    %
    % OUTPUT:
    %   t_sol  - Vettore dei tempi.
    %   p_sol  - Matrice (num_steps+1 x N) delle concentrazioni di proteine sane.
    %   pt_sol - Matrice (num_steps+1 x N) delle concentrazioni di proteine misfolded.
    %
    % Se seedNode è fornito, il nodo (o i nodi) corrispondente verrà
    % seminato con 0.1, altrimenti (default) si utilizza la semina per "Entorhinal".
    
        N = size(A, 1);
    
        % Calcola il Laplaciano
        D = diag(sum(A,2));
        L = D - A;
    
        % Stato iniziale: proteine sane in p0 e proteine misfolded a zero
        p0 = ones(N,1) * (k0 / k1);
        pt0 = zeros(N,1);
        
        if nargin >= 10 && ~isempty(seedNode)
            pt0(seedNode) = 0.1;
        else
            % Default: semina nella regione "Entorhinal"
            infected_mask = strcmp(CoordTable{:,4}, 'Entorhinal');
            pt0(infected_mask) = 0.1;
        end
    
        % Vettore temporale
        t_sol = linspace(0, dt*num_steps, num_steps+1);
    
        % Prealloca le matrici soluzioni
        p_sol = zeros(num_steps+1, N);
        pt_sol = zeros(num_steps+1, N);
        p_sol(1,:) = p0';
        pt_sol(1,:) = pt0';
    
        % Definisci le funzioni ODE:
        % dp/dt = -diffusion_coeff * L * p + k0 - k1 * p - k12*(p .* pt)
        f_p = @(p, pt) -diffusion_coeff * L * p + k0 - k1 * p - k12 * (p .* pt);
        % dpt/dt = -diffusion_coeff * L * pt - ktilde1 * pt + k12*(p .* pt)
        f_pt = @(p, pt) -diffusion_coeff * L * pt - ktilde1 * pt + k12 * (p .* pt);
    
        % Integrazione temporale (metodo di Eulero)
        for i = 2:length(t_sol)
            p_old = p_sol(i-1,:)';
            pt_old = pt_sol(i-1,:)';
    
            dpdt = f_p(p_old, pt_old);
            dptdt = f_pt(p_old, pt_old);
    
            p_new = p_old + dt * dpdt;
            pt_new = pt_old + dt * dptdt;
    
            p_sol(i,:) = p_new';
            pt_sol(i,:) = pt_new';
        end
    end
    
