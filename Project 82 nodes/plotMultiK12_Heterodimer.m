function plotMultiK12_Heterodimer(A, CoordTable, k0, k1, ktilde1, diffusion_coeff, dt, num_steps)
%PLOTMULTIK12_HETERODIMER Runs heterodimer simulations with varying k12 and plots the 
% average misfolded protein concentration (pt_sol) for each simulation.
%
%   plotMultiK12_Heterodimer(A, CoordTable, k0, k1, ktilde1, diffusion_coeff, dt, num_steps)
%
%   INPUTS:
%       A              - Adjacency matrix.
%       CoordTable     - Table containing node data (with region labels in column 4).
%       k0             - Production rate of healthy protein.
%       k1             - Clearance rate of healthy protein.
%       ktilde1       - Clearance rate of misfolded protein.
%       diffusion_coeff- Diffusion coefficient.
%       dt             - Time step size (years).
%       num_steps      - Number of time steps for the simulation.
%
%   The function varies k12 from 0 to 0.5 (in steps of 0.01). For each k12, it
%   calls HeterodimerInfection to run the simulation, computes the average misfolded
%   protein concentration over time, and plots the result using a color determined by 
%   the current k12 value. Bold (thicker) lines are used for k12 = 0.00, 0.10, 0.20, 
%   0.30, 0.40, and 0.50. A legend is added for these special values.
%
% Example call from your main script:
%   plotMultiK12_Heterodimer(A, CoordTable, 1.0, 0.5, 0.5, 0.05, 0.4, 100);

    % Define k12 values from 0 to 0.5 in steps of 0.01
    k12_values = 0.25:0.01:0.5;
    % Special k12 values for which we want a legend
    special_k12 = [0.25, 0.3, 0.4, 0.5];
    special_handles = [];
    special_labels = {};
    
    % Create a new figure for the plot
    figure;
    hold on;
    
    % Loop over all k12 values
    for k12 = k12_values
        % Run the heterodimer simulation for the current k12 value.
        % HeterodimerInfection is assumed to be defined on the MATLAB path.
        [t_sol, ~, pt_sol] = HeterodimerInfection(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps);
        
        % Compute the average misfolded protein concentration over all nodes
        avg_pt = mean(pt_sol, 2);
        
        % Determine the color and line width based on the current k12 value using a switch
        [thisColor, thisLineWidth] = getColorAndLineWidth(k12);
        
        % Plot the average concentration vs. time for this simulation and get the handle
        h = plot(t_sol, avg_pt, 'Color', thisColor, 'LineWidth', thisLineWidth);
        
        % If the current k12 is one of the special values, store its handle and label
        if ismember(k12, special_k12)
            special_handles(end+1) = h; %#ok<AGROW>
            special_labels{end+1} = sprintf('k_{12} = %.2f', k12); %#ok<AGROW>
        end
    end
    
    % Add legend for the special k12 values
    legend(special_handles, special_labels, 'Location', 'best');
    
    % Add labels, title, and grid to the plot
    xlabel('Time [yr]');
    ylabel('Average Misfolded Protein Concentration');
    title('Plot average misfolded concentration with different values of k_{12}');
    grid on;
    hold off;
end

%% Subfunction: getColorAndLineWidth
function [col, lw] = getColorAndLineWidth(k12)
%GETCOLORANDLINEWIDTH Returns the plot color and line width for a given k12 value
%   using a switch statement based on floor(k12*10).
%
%   [col, lw] = getColorAndLineWidth(k12)
%
%   Color assignment:
%       Case 2:  Green       (k12 in [0.25, 0.29], bold if exactly 0.20)
%       Case 3:  Yellow      (k12 in [0.3, 0.39], bold if exactly 0.30)
%       Case 4:  Orange      (k12 in [0.4, 0.49], bold if exactly 0.40)
%       Case 5:  Red         (k12 = 0.50, always bold)

    % Default line width
    lw = 1;
    switch floor(k12 * 10)
        case 2
            col = [0 1 0]; % green
            if k12 == 0.25
                lw = 2;
            end
        case 3
            col = [1 1 0]; % yellow
            if k12 == 0.3
                lw = 2;
            end
        case 4
            col = [1 0.5 0]; % orange
            if k12 == 0.4
                lw = 2;
            end
        case 5
            col = [1 0 0]; % red
            lw = 2;
        otherwise
            col = [0 0 0]; % default to black (should not occur)
    end
end
