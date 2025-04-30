% HeterodimerInfection_dynamic.m
function [t_sol, p_sol, pt_sol] = HeterodimerInfection_treatment_dynamic(...
    A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt, num_steps, p_init, pt_init)
% HETERODIMERINFECTION_DYNAMIC  Come HeterodimerInfection, ma ad ogni passo:
%   – riduce casualmente X archi (tra 20% e 90% degli archi) moltiplicando i
%     loro pesi per edge_reduction;
%   – moltiplica k12 per un fattore random ∈ [0.9,1.0].
%
% INPUT aggiuntivi:
%   edge_reduction – fattore di riduzione dei pesi degli archi;
%   p_init, pt_init – (opzionali) condizioni iniziali per p e pt.
%
% OUTPUT:
%   t_sol  – vettore dei tempi (0:dt:num_steps*dt)
%   p_sol, pt_sol – matrici (numel(t_sol)×N) di concentrazioni.

    N = size(A,1);

    % --- condizioni iniziali ---
    if nargin < 11 || isempty(p_init)
        p0 = ones(N,1) * (k0 / k1);
    else
        p0 = p_init;
    end
    if nargin < 12 || isempty(pt_init)
        pt0 = zeros(N,1);
        infected_mask = strcmp(CoordTable{:,4}, 'Entorhinal');
        pt0(infected_mask) = 0.1;
    else
        pt0 = pt_init;
    end

    t_sol = linspace(0, dt*num_steps, num_steps+1);
    p_sol  = zeros(numel(t_sol), N);
    pt_sol = zeros(numel(t_sol), N);
    p_sol(1,:)  = p0';
    pt_sol(1,:) = pt0';

    disp('Starting dynamic Heterodimer simulation...');
    try
        for step = 2:numel(t_sol)
            %% --- dynamic edge reduction ---
            [iList,jList] = find(triu(A,1));
            E = numel(iList);
            if E>0
                minX = floor(0.2*E);
                maxX = ceil (0.9*E);
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

            %% --- stochastic conversion reduction ---
            k12 = k12 * (rand()*0.1 + 0.9);

            %% --- integrazione (Forward Euler) ---
            L = diag(sum(A,2)) - A;
            p_prev  = p_sol(step-1, :)';
            pt_prev = pt_sol(step-1, :)';

            dpdt  = -diffusion_coeff * L * p_prev ...
                    + k0 - k1 * p_prev - k12 * (p_prev .* pt_prev);
            dptdt = -diffusion_coeff * L * pt_prev ...
                    - ktilde1 * pt_prev + k12 * (p_prev .* pt_prev);

            p_new  = p_prev  + dt * dpdt;
            pt_new = pt_prev + dt * dptdt;

            p_sol(step, :)  = p_new';
            pt_sol(step, :) = pt_new';
        end
        disp('Dynamic Heterodimer simulation completed.');
    catch ME
        rethrow(ME);
    end
end
