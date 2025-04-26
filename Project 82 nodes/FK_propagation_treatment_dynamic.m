function [t_sol, c_sol] = FK_propagation_treatment_dynamic(A, CoordTable, diffusion, a, dt, num_steps, edge_reduction, c_init)
% FK_PROPAGATION_DYNAMIC  As before, but also applies a random reduction
%                        to the growth rate a at each time‐step.

    N = size(A,1);

    % --- initial condition ---
    if nargin<8 || isempty(c_init)
        c0 = zeros(N,1);
        infected = strcmp(CoordTable{:,4}, 'Entorhinal');
        c0(infected) = 0.1;
    else
        c0 = c_init;
    end

    t_sol = linspace(0,dt*num_steps,num_steps+1);
    c_sol = zeros(num_steps+1, N);
    c_sol(1,:) = c0';

    disp('Starting simulation with stochastic edge & growth reduction...');
    try
        for step = 2:numel(t_sol)
            %% --- dynamic edge reduction (same as before) ---
            [iList,jList] = find(triu(A,1));
            E = numel(iList);
            if E>0
                minX = floor(0.2*E);
                maxX = ceil (0.9*E);
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

            %% --- stochastic growth reduction ---
            % uniform [0,1] factor; change to rand*(0.9-0.8)+0.8 if you
            % want Uniform[0.8,0.9]
            a = a * (rand()*0.1+0.9);

            %% --- build Laplacian and do Euler step ---
            L = diag(sum(A,2)) - A;
            c_prev = c_sol(step-1, :)';
            dc     = -L*c_prev*diffusion + a*c_prev.*(1-c_prev);
            c_next = c_prev + dt*dc;
            c_sol(step,:) = c_next';
        end
        disp('Dynamic simulation completed.');
    catch ME
        rethrow(ME);
    end
end
