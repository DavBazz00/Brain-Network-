function plotMultiAlpha_FK(A, CoordTable, diffusion, dt, num_steps)
%PLOTMULTIALPHA Runs simulations with varying α and plots the average infection 
% concentration for each simulation. A legend is added for the six special 
% α values: 0.00, 0.10, 0.20, 0.30, 0.40, and 0.50.
%
%   plotMultiAlpha(A, CoordTable, diffusion, dt, num_steps)
%
%   INPUTS:
%       A          - Adjacency matrix (preprocessed as needed).
%       CoordTable - Table containing node coordinates and region labels (with
%                    region info in the 4th column).
%       diffusion  - Diffusion coefficient used in FK_propagation.
%       dt         - Time step size (years).
%       num_steps  - Number of time steps for the simulation.
%
%   This function loops over α values from 0 to 0.5 (in steps of 0.01). For each α,
%   it calls FK_propagation to run the simulation, computes the average infection
%   concentration over time, and plots the result using a color determined by the 
%   current α value (using a switch in the helper function). Bold (thicker) lines are 
%   used for α values 0.00, 0.10, 0.20, 0.30, 0.40, and 0.50, and these are also added 
%   to the legend.
%

    % Define α values from 0 to 0.5 in steps of 0.01
    alpha_values = 0:0.01:0.5;
    % Special α values for which we want a legend
    special_alphas = [0, 0.1, 0.2, 0.3, 0.4, 0.5];
    special_handles = [];
    special_labels = {};
    
    % Create a new figure for the plot
    figure;
    hold on;

    % Loop over all α values
    for alpha = alpha_values
        % Call FK_propagation with the current α as the logistic growth parameter
        [t_sol, c_sol] = FK_propagation(A, CoordTable, diffusion, alpha, dt, num_steps);

        % Compute average concentration over all nodes at each time point
        avg_conc = mean(c_sol, 2);

        % Determine color and line width based on the current α using switch
        [thisColor, thisLineWidth] = getColorAndLineWidth(alpha);

        % Plot the average concentration vs. time for this simulation and get the handle
        h = plot(t_sol, avg_conc, 'Color', thisColor, 'LineWidth', thisLineWidth);
        
        % If the current alpha is one of the special values, store its handle and label
        if ismember(alpha, special_alphas)
            special_handles(end+1) = h; %#ok<AGROW>
            special_labels{end+1} = sprintf('\\alpha = %.2f', alpha); %#ok<AGROW>
        end
    end

    % Add legend for the special α values
    legend(special_handles, special_labels, 'Location', 'best');
    
    % Add labels, title, and grid to the plot
    xlabel('Time [yr]');
    ylabel('Average Infection Concentration');
    title('Plot average concentration with different values of \alpha');
    grid on;
    hold off;
end

%% Subfunction: getColorAndLineWidth
function [col, lw] = getColorAndLineWidth(alpha)
%GETCOLORANDLINEWIDTH Returns the plot color and line width for a given α value
%   using a switch statement on floor(alpha*10).
%
%   [col, lw] = getColorAndLineWidth(alpha)
%
%   This function uses a switch statement on floor(alpha*10) to select the
%   appropriate color. Bold lines (line width 2) are used for the special 
%   α values (0.00, 0.10, 0.20, 0.30, 0.40, 0.50).

    % Default line width
    lw = 1;
    % Determine the case using floor(alpha*10)
    switch floor(alpha*10)
        case 0
            col = [0 0 1]; % blue
            if alpha == 0
                lw = 2;
            end
        case 1
            col = [0 0.7 1]; % light-blue
            if alpha == 0.1
                lw = 2;
            end
        case 2
            col = [0 1 0]; % green
            if alpha == 0.2
                lw = 2;
            end
        case 3
            col = [1 1 0]; % yellow
            if alpha == 0.3
                lw = 2;
            end
        case 4
            col = [1 0.5 0]; % orange
            if alpha == 0.4
                lw = 2;
            end
        case 5
            col = [1 0 0]; % red
            lw = 2;
        otherwise
            col = [0 0 0]; % default to black (should not occur)
    end
end
