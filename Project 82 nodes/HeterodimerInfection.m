function [t_sol, p_sol, pt_sol] = HeterodimerInfection(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps)
%HETERODIMERINFECTION Simulates the heterodimer infection propagation model.
%
%   [t_sol, p_sol, pt_sol] = heterodimerInfection(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps)
%
%   INPUTS:
%       A              - Adjacency matrix.
%       CoordTable     - Table containing node data, including region labels in column 4.
%       k0             - Production rate of healthy protein.
%       k1             - Clearance rate of healthy protein.
%       ktilde1       - Clearance rate of misfolded protein.
%       k12            - Conversion rate from healthy to misfolded protein.
%       diffusion_coeff- Diffusion coefficient.
%       dt             - Time step size (in years).
%       num_steps      - Number of time steps for the simulation.
%
%   OUTPUTS:
%       t_sol  - Time vector for the simulation.
%       p_sol  - Matrix containing the healthy protein concentration over time.
%                Each row corresponds to a time point.
%       pt_sol - Matrix containing the misfolded protein concentration over time.
%                Each row corresponds to a time point.
%
%   This function sets the initial conditions based on a healthy state for
%   protein p and seeds the misfolded protein pt in the "Entorhinal" region.
%   It then simulates the heterodimer infection propagation using a forward
%   Euler integration of the following ODEs:
%
%       dp/dt  = -diffusion_coeff * L * p + k0 - k1 * p - k12 * (p .* pt)
%       dpt/dt = -diffusion_coeff * L * pt - ktilde1 * pt + k12 * (p .* pt)
%
%   where L is the Laplacian computed from the adjacency matrix A.
%

    %% Compute the Laplacian matrix from A
    N = size(A, 1);
    D = diag(sum(A, 2));
    L = D - A;
    
    %% Set Initial Conditions
    % Healthy protein initial condition (healthy state)
    p0 = ones(N, 1) * (k0 / k1);
    
    % Misfolded protein initial condition: zero everywhere
    pt0 = zeros(N, 1);
    % Seed misfolded protein in the "Entorhinal" region (assumed stored in column 4)
    infected_mask = strcmp(CoordTable{:, 4}, 'Entorhinal');
    pt0(infected_mask) = 0.1;
    
    %% Time parameters and preallocation
    tmax = dt * num_steps;
    t_sol = linspace(0, tmax, num_steps + 1);
    
    p_sol = zeros(num_steps + 1, N);
    pt_sol = zeros(num_steps + 1, N);
    p_sol(1, :) = p0';
    pt_sol(1, :) = pt0';
    
    %% Define ODE functions for the heterodimer model
    % dp/dt = -diffusion_coeff * L * p + k0 - k1 * p - k12 * (p .* pt)
    f_p = @(p, pt) -diffusion_coeff * L * p + k0 - k1 * p - k12 * (p .* pt);
    % dpt/dt = -diffusion_coeff * L * pt - ktilde1 * pt + k12 * (p .* pt)
    f_pt = @(p, pt) -diffusion_coeff * L * pt - ktilde1 * pt + k12 * (p .* pt);
    
    %% Simulate Heterodimer Infection Propagation using Forward Euler
    disp('Starting Heterodimer simulation...');
    for i = 2:length(t_sol)
        % Extract previous time step concentrations
        p_old = p_sol(i-1, :)';
        pt_old = pt_sol(i-1, :)';
        
        % Compute derivatives
        dpdt = f_p(p_old, pt_old);
        dptdt = f_pt(p_old, pt_old);
        
        % Update concentrations using Euler integration
        p_new = p_old + dt * dpdt;
        pt_new = pt_old + dt * dptdt;
        
        % Store updated concentrations
        p_sol(i, :) = p_new';
        pt_sol(i, :) = pt_new';
    end
    disp('Heterodimer simulation completed successfully.');
    
end
