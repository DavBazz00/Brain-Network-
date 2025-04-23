function [t, pt_sol] = HeterodimerInfection_combined( ...
    A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, ...
    dt_aging, dt_treat, T_end, t_switch, aging_red, cure_red)
% HETERODIMERINFECTION_COMBINED Simula aging dinamico + trattamento annuale
%
% INPUTS:
%   A               – matrice di adiacenza iniziale
%   CoordTable      – tabella con info sui nodi (colonna 4: regione)
%   k0,k1,ktilde1,k12– parametri del modello heterodimer :contentReference[oaicite:0]{index=0}&#8203;:contentReference[oaicite:1]{index=1}
%   diffusion_coeff – coefficiente di diffusione
%   dt_aging        – intervallo di aging (anni)
%   dt_treat        – intervallo di cura (anni)
%   T_end           – durata simulazione (anni)
%   t_switch        – anno inizio trattamento
%   aging_red       – fattore moltiplicativo per aging (0–1)
%   cure_red        – fattore moltiplicativo per cura (0–1)
%
% OUTPUT:
%   t       – vettore tempi [0:dt_aging:T_end]
%   pt_sol  – matrice misfolded-protein (length(t)×N)

% timeline
t = 0:dt_aging:T_end;
num_steps = length(t);
N = size(A,1);

% preallocazione
p_sol  = zeros(num_steps, N);
pt_sol = zeros(num_steps, N);

% condizione iniziale
p0  = ones(N,1)*(k0/k1);
pt0 = zeros(N,1);
infected_mask = strcmp(CoordTable{:,4}, 'Entorhinal');
pt0(infected_mask) = 0.1;
p_sol(1,:)  = p0';
pt_sol(1,:) = pt0';

% stato di lavoro della rete
A_curr    = A;
next_cure = t_switch;

% funzioni ODE heterodimer
f_p  = @(L,p,pt)  -diffusion_coeff*L*p + k0 - k1*p - k12*(p.*pt);
f_pt = @(L,p,pt)  -diffusion_coeff*L*pt - ktilde1*pt + k12*(p.*pt);

% loop di integrazione con aging+cura
for i = 2:num_steps
    t_now = t(i);

    % --- 1) Aging dinamico (ogni dt_aging) ---
    [ei, ej] = find(triu(A_curr,1));
    M = numel(ei);
    Kmax = ceil(0.8 * M);
    if Kmax >= 2
        if mod(Kmax,2)==1, Kmax=Kmax-1; end
        K = 2 * randi([1, Kmax/2]);
        sel = randperm(M, K);
        for k = 1:2:K
            u = ei(sel(k)); v = ej(sel(k+1));
            A_curr(u,v) = A_curr(u,v) * aging_red;
            A_curr(v,u) = A_curr(v,u) * aging_red;
        end
    end

    % --- 2) Trattamento annuale (ogni dt_treat), da t_switch in poi ---
    if t_now >= next_cure
        A_curr = A_curr * cure_red;
        next_cure = next_cure + dt_treat;
    end

    % --- 3) Forward Euler sul modello heterodimer ---
    L = diag(sum(A_curr,2)) - A_curr;
    p_prev  = p_sol(i-1, :)';
    pt_prev = pt_sol(i-1,: )';
    p_new  = p_prev  + dt_aging * f_p(L, p_prev,  pt_prev);
    pt_new = pt_prev + dt_aging * f_pt(L, p_prev, pt_prev);
    p_sol(i,:)  = p_new';
    pt_sol(i,:) = pt_new';
end
end
