function plotSpecificNodePairs(t_sol, c_sol, CoordTable, plotTitle)
%PLOTSPECIFICNODEPAIRS Plots the average infection concentration for specific node pairs.
%
%   plotSpecificNodePairs(t_sol, c_sol, CoordTable, plotTitle)
%
%   INPUTS:
%       t_sol     - Time vector for the simulation.
%       c_sol     - Matrix of infection concentrations over time.
%       CoordTable- Table containing node data, including the 'Label' column.
%       plotTitle - (Optional) Title for the plot.
%
%   The function computes the average concentration over time for the following 
%   node pairs:
%       Temporal:  '22L', '22R' (plotted in green)
%       Frontal:   '10L', '10R' (plotted in red)
%       Parietal:  '5L',  '5R'  (plotted in orange)
%       Occipital: '19L', '19R' (plotted in blue)
%
%   Example:
%       plotSpecificNodePairs(t_sol, c_sol, CoordTable, 'Average Concentration Over Time');
%

    % If no title is provided, set a default
    if nargin < 4 || isempty(plotTitle)
        plotTitle = 'Average Concentration Over Time for Specific Node Pairs';
    end

    %% Define the node labels for each region of interest
    temp_labels    = {'22L', '22R'};   % Temporal
    frontal_labels = {'10L', '10R'};   % Frontal
    parietal_labels= {'5L',  '5R'};    % Parietal
    occipital_labels = {'19L','19R'};   % Occipital

    %% Get indices for each pair using the 'Label' column in CoordTable
    temp_idx    = ismember(CoordTable.Label, temp_labels);
    frontal_idx = ismember(CoordTable.Label, frontal_labels);
    parietal_idx= ismember(CoordTable.Label, parietal_labels);
    occipital_idx = ismember(CoordTable.Label, occipital_labels);

    %% Compute the average concentration over time for each region
    avg_temp     = mean(c_sol(:, temp_idx), 2);
    avg_frontal  = mean(c_sol(:, frontal_idx), 2);
    avg_parietal = mean(c_sol(:, parietal_idx), 2);
    avg_occipital= mean(c_sol(:, occipital_idx), 2);

    %% Plot the averages with the specified colors
    figure;
    hold on;
    plot(t_sol, avg_temp,     'LineWidth', 2, 'Color', [0 1 0]);   % Green for Temporal
    plot(t_sol, avg_frontal,  'LineWidth', 2, 'Color', [1 0 0]);   % Red for Frontal
    plot(t_sol, avg_parietal, 'LineWidth', 2, 'Color', [1 0.5 0]); % Orange for Parietal
    plot(t_sol, avg_occipital,'LineWidth', 2, 'Color', [0 0 1]);   % Blue for Occipital

    xlabel('Time (years)');
    ylabel('Average Misfolded Protein Concentration');
    title(plotTitle);
    legend({'Temporal (22L, 22R)', 'Frontal (10L, 10R)', 'Parietal (5L, 5R)', 'Occipital (19L, 19R)'}, 'Location', 'best');
    grid on;
    hold off;
end
