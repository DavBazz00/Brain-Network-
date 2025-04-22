function plotSpecificNodePairs(t_sol, c_sol, CoordTable, plotTitle)
%PLOTSPECIFICNODEPAIRS Plots the average infection concentration for specific node groups.
%
%   plotSpecificNodePairs(t_sol, c_sol, CoordTable, plotTitle)
%
%   INPUTS:
%       t_sol     - Time vector for the simulation.
%       c_sol     - Matrix of infection concentrations over time.
%       CoordTable- Table containing node data, including the 'Label' column.
%       plotTitle - (Optional) Title for the plot.
%
%   This function computes the average concentration over time for the following 
%   node groups:
%       Temporal:  {'21L', '21R','30L','30R','22L','22R'}
%       Frontal:   {'8L','8R','10L','10R','47L','47R'}
%       Parietal:  {'5L','5R','7L','7R', '40L', '40R'}
%       Occipital: {'18L','18R','19L','19R'}

    if nargin < 4 || isempty(plotTitle)
        plotTitle = 'Average Concentration Over Time for Specific Node Groups';
    end

    % Define node labels for each region
    temp_labels    = {'21L', '21R','30L','30R','22L','22R'};   % Temporal
    frontal_labels = {'8L','8R','10L','10R','47L','47R'};   % Frontal
    parietal_labels= {'5L','5R','7L','7R', '40L', '40R'};    % Parietal
    occipital_labels = {'18L','18R','19L','19R'};   % Occipital

    % Get logical indices for each group
    temp_idx     = ismember(CoordTable.Label, temp_labels);
    frontal_idx  = ismember(CoordTable.Label, frontal_labels);
    parietal_idx = ismember(CoordTable.Label, parietal_labels);
    occipital_idx= ismember(CoordTable.Label, occipital_labels);

    % Compute average concentrations over time
    avg_temp     = mean(c_sol(:, temp_idx), 2);
    avg_frontal  = mean(c_sol(:, frontal_idx), 2);
    avg_parietal = mean(c_sol(:, parietal_idx), 2);
    avg_occipital= mean(c_sol(:, occipital_idx), 2);

    % Plot
    figure;
    hold on;
    plot(t_sol, avg_temp,     'LineWidth', 2, 'Color', [0 1 0]);   % Green for Temporal
    plot(t_sol, avg_frontal,  'LineWidth', 2, 'Color', [1 0 0]);   % Red for Frontal
    plot(t_sol, avg_parietal, 'LineWidth', 2, 'Color', [1 0.5 0]); % Orange for Parietal
    plot(t_sol, avg_occipital,'LineWidth', 2, 'Color', [0 0 1]);   % Blue for Occipital

    xlabel('Time (years)');
    ylabel('Average Misfolded Protein Concentration');
    title(plotTitle);
    legend({ ...
        'Temporal (21L, 21R, 22L, 22R, 30L, 30R)', ...
        'Frontal (8L, 8R, 10L, 10R, 47L, 47R)', ...
        'Parietal (5L, 5R, 7L, 7R, 40L, 40R)', ...
        'Occipital (18L, 18R, 19L, 19R)'}, 'Location', 'best');
    grid on;
    hold off;
end
