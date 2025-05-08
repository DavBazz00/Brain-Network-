function [t_sol, p_sol, pt_sol] = HeterodimerInfection_treatment(...
    A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, ...
    dt1, dt2, t_switch, t_end)
% HETERODIMERINFECTION_TREATMENT
%   Simula il modello eterodimero in due fasi:
%     Fase 1: infezione con passo dt1 fino a t = t_switch
%     Fase 2: trattamento con passo dt2 da t = t_switch a t = t_end
%
% INPUT:
%   A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction
%   dt1, dt2, t_switch, t_end
%
% OUTPUT:
%   t_sol       : vettore dei tempi da 0 a t_end (colonna)
%   p_sol, pt_sol: matrici (numel(t_sol)×N) delle concentrazioni

    N = size(A,1);

    % numero di step di ciascuna fase
    num_steps1 = round(t_switch    / dt1);
    num_steps2 = round((t_end - t_switch) / dt2);

    % ------------------------
    % FASE 1: propagazione
    % Step 1: evoluzione fino a prima di t_switch
    [t1a, p1a, pt1a] = HeterodimerInfection( ...
        A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt1, num_steps1);
    t1a = t1a(:);
    
    % Step 2: eventuale passo aggiuntivo per arrivare a t_switch
    t_partial = t1a(end);
    dt_extra = t_switch - t_partial;
    
    if dt_extra > 1e-8
        % iniziali per passo extra
        p_init = p1a(end, :)';
        pt_init = pt1a(end, :)';
        L = diag(sum(A,2)) - A;
        dpdt = -diffusion_coeff*(L*p_init) + k0 - k1*p_init - k12*(p_init.*pt_init);
        dptdt = -diffusion_coeff*(L*pt_init) - ktilde1*pt_init + k12*(p_init.*pt_init);
        p_extra = p_init + dt_extra * dpdt;
        pt_extra = pt_init + dt_extra * dptdt;
    
        % aggiorna output
        t1 = [t1a; t_switch];
        p1 = [p1a; p_extra'];
        pt1 = [pt1a; pt_extra'];
    else
        % già perfettamente allineato
        t1 = t1a;
        p1 = p1a;
        pt1 = pt1a;
    end


    % iniziali per fase 2
    p_init  = p1(end, :)';
    pt_init = pt1(end, :)';

    % ------------------------
    % FASE 2.1: trattamento
    p2  = zeros(num_steps2+1, N);
    pt2 = zeros(num_steps2+1, N);
    t2  = zeros(num_steps2+1, 1);

    p2(1,:)  = p_init';
    pt2(1,:) = pt_init';
    t2(1)    = t1(end);

    k12_eff = k12;

    for step = 2:(num_steps2+1)
        % 1) dynamic edge variation (treatment)
        [iList,jList] = find(triu(A,1));
        E = numel(iList);
            if E>0
                minX = floor(0.2*E);
                maxX = ceil (0.5*E);
                % draw X ∈ {minX,…,maxX}
                X = minX + randi(max(1,maxX-minX+1))-1;

                if X>0
                    linIdx = sub2ind([N,N], iList, jList);
                    w     = A(linIdx);
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

        % 2) stochastic conversion reduction (treatment)
        k12_eff = k12_eff * (rand()*0.1 + 0.9);

        % 3) integrazione Forward Euler
        L      = diag(sum(A,2)) - A;
        p_old  = p2(step-1, :)';
        pt_old = pt2(step-1, :)';
        dpdt   = -diffusion_coeff*(L*p_old) + k0 - k1*p_old - k12_eff*(p_old.*pt_old);
        dptdt  = -diffusion_coeff*(L*pt_old) - ktilde1*pt_old + k12_eff*(p_old.*pt_old);
        p_new  = p_old  + dt2 * dpdt;
        pt_new = pt_old + dt2 * dptdt;

        p2(step, :)  = p_new';
        pt2(step, :) = pt_new';
        t2(step)     = t2(step-1) + dt2;
    end

    % ------------------------
    % concatenazione (evita duplicare t1(end))
    t_sol  = [t1;      t2(2:end)];      % entrambi vettori colonna
    p_sol  = [p1;      p2(2:end, :)];
    pt_sol = [pt1;     pt2(2:end, :)];
end
