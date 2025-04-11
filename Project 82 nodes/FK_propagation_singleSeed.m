function [t_sol, c_sol] = FK_propagation_singleSeed(A, CoordTable, diffusion, a, dt, num_steps, seedNode)
%FK_PROPAGATION  Fisher-Kolmogorov infection propagation
%   with optional single-node seeding.

    % 1) Number of nodes
    N = size(A, 1);

    % 2) Initial Conditions: all zero
    c0 = zeros(N, 1);

    % If 'seedNode' is provided, set that node to 0.1
    % If 'seedNode' is missing, default to Entorhinal region
    if nargin < 7 || isempty(seedNode)
        infected_mask = strcmp(CoordTable{:,4}, 'Entorhinal');
        c0(infected_mask) = 0.1;
    else
        c0(seedNode) = 0.1; 
    end

    % 3) Laplacian as before
    D = diag(sum(A,2));
    L = D - A;

    % 4) Time vector
    t_sol = linspace(0, dt*num_steps, num_steps + 1);

    % 5) Allocate
    c_sol = zeros(num_steps+1, N);
    c_sol(1,:) = c0;

    % 6) ODE forward-Euler:
    f = @(c) -diffusion*(L*c) + a*(c.*(1-c));

    for i = 1:num_steps
        c_now = c_sol(i,:)';
        c_next = c_now + dt*f(c_now);

        % optional clamp to [0..1]
        c_next(c_next < 0) = 0;
        c_next(c_next > 1) = 1;

        c_sol(i+1,:) = c_next;
    end
end
