function [t_sol, p_sol, pt_sol] = HeterodimerInfection_dynamic(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt, num_steps)
% HETERODIMERINFECTION_DYNAMIC Simulates the heterodimer infection model with dynamic edge variation.
%
%   [t_sol, p_sol, pt_sol] = heterodimerInfection_dynamic(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt, num_steps)
%
%   INPUTS:
%       A               - Adjacency matrix.
%       CoordTable      - Table containing node data, including region labels in column 4.
%       k0              - Production rate of healthy protein.
%       k1              - Clearance rate of healthy protein.
%       ktilde1         - Clearance rate of misfolded protein.
%       k12             - Conversion rate from healthy to misfolded protein.
%       diffusion_coeff - Diffusion coefficient.
%       edge_reduction  - Factor for reducing edge weights during dynamic update.
%       dt              - Time step size (in years).
%       num_steps       - Number of time steps for the simulation.
%
%   OUTPUTS:
%       t_sol   - Time vector for the simulation.
%       p_sol   - Matrix of healthy protein concentrations over time.
%                 Each row corresponds to a time point.
%       pt_sol  - Matrix of misfolded protein concentrations over time.
%                 Each row corresponds to a time point.
%
%   This function initializes the healthy protein state as p0 = k0/k1 and seeds the
%   misfolded protein (pt0) in the "Entorhinal" region (as specified in CoordTable).
%   It then simulates the heterodimer infection model using a forward Euler integration,
%   while dynamically modifying the network edges at each time step.
%
% Example:
%   [t_sol, p_sol, pt_sol, A_final] = heterodimerInfection_dynamic(A, CoordTable, 1.0, 0.5, 0.5, 0.5, 0.05, 0.7, 0.4, 100);

    % Determine the number of nodes
    N = size(A, 1);

    %% Initial Conditions
    % Healthy protein initial condition: p0 = k0/k1 (healthy state)
    p0 = ones(N, 1) * (k0 / k1);
    
    % Misfolded protein initial condition: initially zero everywhere
    pt0 = zeros(N, 1);
    % Seed misfolded protein in the "Entorhinal" region (assumed stored in column 4)
    infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
    pt0(infected_mask) = 0.1;
    
    %% Time Parameters and Preallocation
    tmax = dt * num_steps;             % Total simulation time
    t_sol = linspace(0, tmax, num_steps+1); % Time vector for simulation
    
    % Preallocate solution matrices for healthy (p) and misfolded (pt) proteins
    p_sol = zeros(num_steps+1, N);
    pt_sol = zeros(num_steps+1, N);
    p_sol(1, :) = p0';
    pt_sol(1, :) = pt0';
    
    %% Define ODE Functions for the Heterodimer Model
    % dp/dt = -diffusion_coeff * L_current * p + k0 - k1 * p - k12 * (p .* pt)
    f_p = @(L_current, p, pt) -diffusion_coeff * L_current * p + k0 - k1 * p - k12 * (p .* pt);
    % dpt/dt = -diffusion_coeff * L_current * pt - ktilde1 * pt + k12 * (p .* pt)
    f_pt = @(L_current, p, pt) -diffusion_coeff * L_current * pt - ktilde1 * pt + k12 * (p .* pt);
    
    %% Simulation with Dynamic Edge Variation
    disp('Starting simulation with dynamic edge variation...');
    try
        for i = 2:length(t_sol)
            % --- Dynamic Update of Edges ---
            % Get upper triangular indices (each edge only once)
            [edgeI, edgeJ] = find(triu(A, 1));
            numEdgesTotal = length(edgeI);
            maxEdgesToModify = ceil(0.8 * numEdgesTotal);
            
            % Ensure we select an even number of edges
            if maxEdgesToModify < 2
                numEdgesToModify = 0;
            else
                if mod(maxEdgesToModify, 2) == 1
                    maxEven = maxEdgesToModify - 1;
                else
                    maxEven = maxEdgesToModify;
                end
                numEdgesToModify = 2 * randi([1, maxEven/2]);
            end
            
            if numEdgesToModify > 0
                selectedIndices = randperm(numEdgesTotal, numEdgesToModify);
                % Update selected edges in pairs
                for k = 1:2:numEdgesToModify
                    idx1 = edgeI(selectedIndices(k));
                    idx2 = edgeJ(selectedIndices(k+1));
                    A(idx1, idx2) = A(idx1, idx2) * edge_reduction;
                    A(idx2, idx1) = A(idx2, idx1) * edge_reduction;
                end
            end
            
            % Update Laplacian with the modified A
            L_current = diag(sum(A, 2)) - A;
            
            % Extract the concentrations from the previous time step
            p_old = p_sol(i-1, :)';
            pt_old = pt_sol(i-1, :)';
            
            % Compute the derivatives using the defined ODE functions
            dpdt = f_p(L_current, p_old, pt_old);
            dptdt = f_pt(L_current, p_old, pt_old);
            
            % Update the concentrations using Forward Euler integration
            p_new = p_old + dt * dpdt;
            pt_new = pt_old + dt * dptdt;
            
            % Store the updated concentrations
            p_sol(i, :) = p_new';
            pt_sol(i, :) = pt_new';
        end
        disp('Simulation completed successfully.');
    catch ME
        rethrow(ME);
    end
end
