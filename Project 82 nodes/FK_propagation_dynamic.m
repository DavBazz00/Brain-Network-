function [t_sol, c_sol] = FK_propagation_dynamic(A, CoordTable, diffusion, a, dt, num_steps, edge_reduction)
%FK_PROPAGATION_DYNAMIC Simulates infection propagation with dynamic edge variation.
%
%   [t_sol, c_sol, A_final, CoordTable] = FK_propagation_dynamic(A, CoordTable, a, dt, num_steps, edge_reduction)
%
%   INPUTS:
%       A             - Adjacency matrix (already loaded, with the first row removed).
%       CoordTable    - Table containing node coordinates and region information.
%                       It is assumed that the region labels are stored in the 4th column.
%       diffusion  - Parameter for slowing the propagation%       a             - Logistic growth rate parameter.
%       dt            - Time step size (in years).
%       num_steps     - Number of time steps for the simulation.
%       edge_reduction- Factor (between 0 and 1) for reducing edge weights during dynamic updates.
%
%   OUTPUTS:
%       t_sol         - Time vector for the simulation.
%       c_sol         - Matrix of infection concentrations over time (each row is a time point).
%       A_final       - The final (modified) adjacency matrix after dynamic updates.
%       CoordTable    - The input coordinate table updated with new node labels.
%
%   This function performs the following steps:
%
%     1. Updates the CoordTable by generating new labels for node pairs.
%     2. Computes the degree matrix and Laplacian.
%     3. Sets initial conditions by seeding nodes in the 'Entorhinal' region with a 
%        concentration of 0.1.
%     4. Runs a simulation with dynamic edge modifications:
%          - At each time step a random even number of edges (up to 20% of the total)
%            is selected and modified (edge weights are reduced by edge_reduction).
%          - The Laplacian is updated and a forward Euler update is applied using the
%            reaction–diffusion model:
%
%               dc/dt = - L_current * c * 5e-4 + a * c .* (1-c)
%
%   Example call from your main file:
%
%       [t_sol, c_sol, A_final, CoordTable] = FK_propagation_dynamic(A, CoordTable, 0.5, 0.4, 100, 0.7);
%

    %% --- Prepare Simulation Parameters ---

    % Time parameters
    tmax = dt * num_steps;
    t_sol = linspace(0, tmax, num_steps + 1);
    N = size(A, 1);
    
    % --- Initial Conditions ---
    c0 = zeros(N, 1);
    % Identify nodes in the 'Entorhinal' region (assumes region info is in column 4)
    infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
    c0(infected_mask) = 0.1;
    
    % Preallocate solution matrix (each row is a time point)
    c_sol = zeros(num_steps + 1, N);
    c_sol(1, :) = c0';
    
    % Define the ODE function for reaction–diffusion update
    % Here, 5e-4 acts as a scaling (diffusion) coefficient
    f = @(L_current, c) - L_current * c * diffusion + a * c .* (1 - c);
    
    %% --- Dynamic Simulation with Edge Variation ---
    disp('Starting simulation with dynamic edge variation...');
    try
        for i = 2:length(t_sol)
            % --- Dynamic Update of Edges (Ensure an Even Number of Modifications) ---
            [edgeI, edgeJ] = find(triu(A, 1));  % Get upper triangular indices (each edge only once)
            numEdgesTotal = length(edgeI);
            maxEdgesToModify = ceil(0.8 * numEdgesTotal); % Give the max number of edges that can be modified for each time
            
            % Ensure we select an even number of edges
            if maxEdgesToModify < 2
                numEdgesToModify = 0;  % Not enough edges to modify in pairs
            else
                if mod(maxEdgesToModify, 2) == 1
                    maxEven = maxEdgesToModify - 1;
                else
                    maxEven = maxEdgesToModify;
                end
                % Generate a random even number between 2 and maxEven
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
            
            % --- Update Laplacian with the Modified A ---
            L_current = diag(sum(A, 2)) - A;
            
            % --- Reaction-Diffusion Update (Forward Euler) ---
            temp = c_sol(i-1, :)';            % Current state as a column vector
            temp = temp + dt * f(L_current, temp); % Euler update
            c_sol(i, :) = temp';              % Save updated state (as a row vector)
        end
        disp('Simulation completed successfully.');
    catch ME
        rethrow(ME);
    end
    
end
