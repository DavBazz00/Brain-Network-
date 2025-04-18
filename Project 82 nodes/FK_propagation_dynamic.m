function [t_sol, c_sol] = FK_propagation_dynamic(A, CoordTable, diffusion, a, dt, num_steps, edge_reduction, c_init)
% FK_PROPAGATION_DYNAMIC Simulates infection propagation with dynamic edge variation.
%
%   [t_sol, c_sol] = FK_propagation_dynamic(A, CoordTable, diffusion, a, dt, num_steps, edge_reduction, c_init)
%
%   INPUTS:
%       A             - Adjacency matrix.
%       CoordTable    - Table containing node coordinates and region information.
%                       Assumes region labels are in the 4th column.
%       diffusion     - Diffusion parameter.
%       a             - Logistic growth rate parameter.
%       dt            - Time step size (in years).
%       num_steps     - Number of simulation steps.
%       edge_reduction- Factor (between 0 and 1) to reduce edge weights dynamically.
%       c_init        - (Optional) Initial concentration (column vector). If not provided or empty,
%                       the standard initial condition is used:
%                           c0 = zeros(N,1);
%                           infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
%                           c0(infected_mask) = 0.1;
%
%   OUTPUTS:
%       t_sol         - Time vector for the simulation.
%       c_sol         - Matrix of infection concentrations over time (each row corresponds to a time point).
%
%   This function performs a dynamic edge update at each time step and integrates the
%   reaction-diffusion model:
%
%         dc/dt = - L*c*diffusion + a*c.*(1-c)
%
%   using a forward Euler integration.
    
    % Number of nodes
    N = size(A, 1);
    
    % --- Initial Conditions: Check if c_init is provided ---
    if nargin < 8 || isempty(c_init)
        % Use the default initial condition: zero everywhere except for nodes
        % in the 'Entorhinal' region, which are seeded with a concentration of 0.1.
        c0 = zeros(N, 1);
        infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
        c0(infected_mask) = 0.1;
    else
        c0 = c_init;
    end

    % Create time vector for the simulation
    t_sol = linspace(0, dt * num_steps, num_steps + 1);
    
    % Preallocate solution matrix and assign the initial condition
    c_sol = zeros(num_steps + 1, N);
    c_sol(1, :) = c0';
    
    % Define the reaction-diffusion update function:
    % dc/dt = - L_current * c * diffusion + a * c .* (1-c)
    f = @(L_current, c) - L_current * c * diffusion + a * c .* (1 - c);
    
    disp('Starting simulation with dynamic edge variation...');
    try
        for i = 2:length(t_sol)
            % --- Dynamic Update of Edges (ensure an even number of modifications) ---
            [edgeI, edgeJ] = find(triu(A, 1));  % Get upper-triangular indices to avoid duplicate edge modifications
            numEdgesTotal = length(edgeI);
            maxEdgesToModify = ceil(0.8 * numEdgesTotal); % Maximum number of edges that can be modified per time step

            if maxEdgesToModify < 2
                numEdgesToModify = 0;
            else
                if mod(maxEdgesToModify, 2) == 1
                    maxEven = maxEdgesToModify - 1;
                else
                    maxEven = maxEdgesToModify;
                end
                % Generate a random even number between 2 and maxEven
                numEdgesToModify = 2 * randi([1, maxEven/2]);
            end

            % Modify selected edges in pairs if any are chosen
            if numEdgesToModify > 0
                selectedIndices = randperm(numEdgesTotal, numEdgesToModify);
                for k = 1:2:numEdgesToModify
                    idx1 = edgeI(selectedIndices(k));
                    idx2 = edgeJ(selectedIndices(k + 1));
                    A(idx1, idx2) = A(idx1, idx2) * edge_reduction;
                    A(idx2, idx1) = A(idx2, idx1) * edge_reduction;
                end
            end
            
            % --- Update Laplacian for the current A ---
            L_current = diag(sum(A, 2)) - A;
            
            % --- Forward Euler Reaction-Diffusion Update ---
            temp = c_sol(i - 1, :)';  % Convert current state to column vector
            temp = temp + dt * f(L_current, temp);
            c_sol(i, :) = temp';      % Store updated state (back as a row vector)
        end
        disp('Dynamic simulation completed successfully.');
    catch ME
        rethrow(ME);
    end
end
