function [t_sol, c_sol] = FK_propagation_treatment(...
    A, CoordTable, diffusion, a, edge_reduction, dt1, dt2, t_switch, t_end)
% FK_PROPAGATION_TREATMENT
%   Simula il modello F-K in due fasi:
%     Fase 1: propagazione (FK) con passo dt1 fino a t = t_switch
%     Fase 2: trattamento con passo dt2 da t = t_switch a t = t_end
%
% INPUT:
%   A               - matrice di adiacenza
%   CoordTable      - tabella con coordinate e regione (colonna 4)
%   diffusion       - coefficiente di diffusione
%   a               - parametro di crescita logistica
%   edge_reduction  - fattore di riduzione dei pesi degli archi in trattamento
%   dt1             - passo temporale fase 1
%   dt2             - passo temporale fase 2
%   t_switch        - tempo (anni) inizio trattamento
%   t_end           - tempo (anni) fine simulazione
%
% OUTPUT:
%   t_sol           - vettore tempi [0, t_end] (colonna)
%   c_sol           - matrice (numel(t_sol)×N) concentrazioni

    N = size(A,1);

    % numero di step
    num_steps1 = round(t_switch    / dt1);
    num_steps2 = round((t_end - t_switch) / dt2);

    % --- Fase 1: propagazione F-K ---
    [t1, c1] = FK_propagation(...
        A, CoordTable, diffusion, a, dt1, num_steps1);
    t1 = t1(:);  % assicura colonna

    % condizioni iniziali per fase 2
    c_init = c1(end, :)';

    % preallocazione fase 2
    c2 = zeros(num_steps2+1, N);
    t2 = zeros(num_steps2+1, 1);
    c2(1, :) = c_init';
    t2(1)    = t1(end);

    % parametro di crescita effettivo
    a_eff = a;

    % --- Fase 2: trattamento F-K ---
    for step = 2:(num_steps2+1)
        % 1) dynamic edge variation (trattamento)
        [iList,jList] = find(triu(A,1));
        E = numel(iList);
        if E>0
            minX = floor(0.2*E);
            maxX = ceil(0.5*E);
            X = minX + randi(max(1, maxX-minX+1)) - 1;
            if X>0
                linIdx = sub2ind([N,N], iList, jList);
                w      = A(linIdx);
                [~, ord] = sort(w,'ascend');
                topE     = ord(end-X+1:end);
                for m = 1:numel(topE)
                    u = iList(topE(m));
                    v = jList(topE(m));
                    A(u,v) = A(u,v) * edge_reduction;
                    A(v,u) = A(v,u) * edge_reduction;
                end
            end
        end

        % 2) stochastic reduction of growth rate
        a_eff = a_eff * (rand()*0.1 + 0.9);

        % 3) Forward Euler integration of F-K
        L      = diag(sum(A,2)) - A;
        c_old  = c2(step-1, :)';
        dcdt   = - diffusion * (L * c_old) + a_eff * (c_old .* (1 - c_old));
        c_new  = c_old + dt2 * dcdt;

        % salva
        c2(step, :) = c_new';
        t2(step)    = t2(step-1) + dt2;
    end

    % --- concatena fasi (evita duplicato di t1(end)) ---
    t_sol = [t1;      t2(2:end)];
    c_sol = [c1;      c2(2:end, :)];
end
