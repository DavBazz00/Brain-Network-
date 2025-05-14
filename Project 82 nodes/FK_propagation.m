function [t_sol, c_sol] = FK_propagation(A, CoordTable, diffusion, a, dt, num_steps)
%SIMULATEINFECTION Simulates infection propagation using a reaction-diffusion model.
%
%   [t_sol, c_sol] = simulateInfection(A, CoordTable, a, dt, num_steps)
%
%   INPUTS:
%       A          - Adjacency matrix (after any necessary preprocessing).
%       CoordTable - Table containing node coordinates and region labels.
%                    The region information is assumed to be in the 4th column.
%       diffusion  - Parameter for slowing the propagation
%       a          - Logistic growth rate parameter.
%       dt         - Time step size (in years).
%       num_steps  - Number of time steps for the simulation.
%
%   OUTPUTS:
%       t_sol - Time vector for the simulation.
%       c_sol - Matrix of infection concentrations over time.
%               Each row corresponds to a time point.
%
%   This function sets the initial condition by seeding the nodes in the 
%   'Entorhinal' region with an infection concentration of 0.1, computes 
%   the Laplacian matrix from A, and performs an implicit (forward Euler)
%   integration using the reaction-diffusion ODE:
%
%       dc/dt = -L*c*diffusion + a*c.*(1-c)
%
%   where L is the Laplacian matrix.

    % Number of nodes
    N = size(A, 1);
    
    % Initial Conditions: set all concentrations to zero
    c0 = zeros(N, 1);
    
    % Identify nodes in the 'Entorhinal' region 
    infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
    % Set initial infection concentration for these nodes
    c0(infected_mask) = 0.1;
    
    % Compute the degree matrix and Laplacian matrix
    D = diag(sum(A, 2));
    L = D - A;
    
    % Total simulation time
    tmax = dt * num_steps;
    % Time vector for the simulation
    t_sol = linspace(0, tmax, num_steps + 1);
    
    % Preallocate solution matrix
    c_sol = zeros(length(t_sol), N);
    c_sol(1, :) = c0;
    
    % Define the ODE function (reaction-diffusion model)
    % dc/dt = - L*c*diffusion + a*c.*(1-c)
    f = @(t, c) -L*c*diffusion + a*c.*(1-c);
    
    disp('Starting explicit time integration...');
    try
        % Perform implicit time integration (here using a forward Euler step)
        for i = 2:length(t_sol)
            t = t_sol(i);
            temp = c_sol(i-1, :)';    % Convert previous state to column vector
            temp = temp + dt * f(t, temp); % Compute the update
            c_sol(i, :) = temp';      % Convert back to row vector and save
        end
        disp('Implicit time integration completed successfully.');
    catch ME
        rethrow(ME);
    end
end
