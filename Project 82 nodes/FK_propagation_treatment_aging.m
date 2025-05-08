function [t_sol, c_sol] = FK_propagation_treatment_aging(...
    A, CoordTable, diffusion, a, edge_reduction, dt1, dt2, t_switch, t_end)
% FK_PROPAGATION_TREATMENT_AGING
%   Simula il modello F-K su rete con due fasi:
%     Fase 1: “aging” dinamico (modifica archi casuali) + crescita logistica
%             passo dt1 fino a t = t_switch
%     Fase 2: trattamento + aging
%             (riduzione pesi top‐X, modifica dinamica, stocastica su a)
%             passo dt2 da t = t_switch a t = t_end
%
% INPUT:
%   A               - matrice di adiacenza (N×N)
%   CoordTable      - tabella con info sui nodi (regione in colonna 4)
%   diffusion       - coefficiente diffusione
%   a               - tasso di crescita logistica
%   edge_reduction  - fattore riduzione peso archi in trattamento
%   dt1             - passo temporale fase 1
%   dt2             - passo temporale fase 2
%   t_switch        - tempo (anni) inizio trattamento
%   t_end           - tempo (anni) fine simulazione
%
% OUTPUT:
%   t_sol           - vettore colonna tempi [0 … t_end]
%   c_sol           - matrice (numel(t_sol)×N) di concentrazioni

    N = size(A,1);

    % numero di step
    num_steps1 = round(t_switch    / dt1);
    num_steps2 = round((t_end - t_switch) / dt2);
    % step dopo lo switch
    num2 = round((t_end - t_switch) / dt2);
    % --- Fase 1: propagazione F-K ---
    [t1, c1] = FK_propagation_dynamic(...
        A, CoordTable, diffusion, a, dt1, num_steps1, edge_reduction);
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

    % --- Fase 2: treatment + aging ---
    for i = 2:num2+1
        % 1.1) ridge treatment: riduci pesi top‐X
        [iList,jList] = find(triu(A,1));
        E = numel(iList);
        if E>0
            minX = floor(0.2*E);
            maxX = ceil(0.5*E);
            X = minX + randi(max(1,maxX-minX+1)) - 1;
            if X>0
                lin = sub2ind([N,N], iList, jList);
                w   = A(lin);
                [~,ord] = sort(w,'ascend');
                top = ord(end-X+1:end);
                for m=1:numel(top)
                    u=iList(top(m)); v=jList(top(m));
                    A(u,v)=A(u,v)*edge_reduction;
                    A(v,u)=A(v,u)*edge_reduction;
                end
            end
        end

        % 1.2) dynamic aging: come in fase 1
        [eI,eJ] = find(triu(A,1));
        E        = numel(eI);
        if E>=2
            maxE = ceil(0.8*E);
            if mod(maxE,2)==1, maxE = maxE-1; end
            m   = 2*randi(maxE/2);
            sel = randperm(E,m);
            for k=1:2:m
                u=eI(sel(k)); v=eJ(sel(k+1));
                A(u,v)=A(u,v)*edge_reduction;
                A(v,u)=A(v,u)*edge_reduction;
            end
        end

        % 2) stochastic reduction of growth a
        a_eff = a_eff * (rand()*0.1 + 0.9);

        % 3) integrazione Forward Euler FK
        L      = diag(sum(A,2)) - A;
        c_old  = c2(i-1, :)';
        dcdt   = -diffusion*(L*c_old) + a_eff*(c_old.*(1-c_old));
        c_new  = c_old + dt2 * dcdt;

        t2(i)   = t2(i-1) + dt2;
        c2(i,:) = c_new';
    end

    % concatena (evita duplicato t1(end))
    t_sol = [t1;     t2(2:end)];
    c_sol = [c1;     c2(2:end,:)];
end
