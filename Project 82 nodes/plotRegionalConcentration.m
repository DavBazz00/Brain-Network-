function plotRegionalConcentration(t_sol, c_sol, CoordTable, plot_title)
%PLOTREGIONALCONCENTRATION Plots the average infection concentration by brain region.
%
%   plotRegionalConcentration(t_sol, c_sol, CoordTable)
%
%   INPUTS:
%       t_sol      - Time vector for the simulation.
%       c_sol      - Matrix of infection concentrations over time.
%       CoordTable - Table containing node data. It is assumed that the 
%                    region information is stored in the 4th column.
%       plot_title - Choose title for the plot
%
%   This function defines four regions of interest (Temporal, Frontal,
%   Parietal, Occipital) along with associated colors. For each region, it
%   computes the average infection concentration over time and plots the
%   results on the same figure.

    % Define regions of interest and colors
    region_names = {'Temporal', 'Frontal', 'Parietal', 'Occipital'};
    region_colors = [0 1 0;    % Green for Temporal
                     1 0 0;    % Red for Frontal
                     1 0.5 0;  % Orange for Parietal
                     0 0 1];   % Blue for Occipital

    num_groups = length(region_names);

    % Extract region labels from CoordTable (assumed to be in the 4th column)
    all_regions = CoordTable{:, 4};

    % Initialize figure
    figure;
    hold on;

    % Loop through each region and compute average concentration over time
    for i = 1:num_groups
        region = region_names{i};

        % Logical index for nodes in this region
        idx = strcmp(all_regions, region);

        % Extract the relevant columns from c_sol and compute the mean
        avg_conc_region = mean(c_sol(:, idx), 2);

        % Plot the average concentration for this region
        plot(t_sol, avg_conc_region, 'LineWidth', 2, 'Color', region_colors(i, :));
    end

    % Add labels, title, legend, and grid
    xlabel('Time');
    ylabel('Average Misfolded Protein Concentration');
    title(plot_title);
    legend(region_names, 'Location', 'northeast');
    grid on;
    hold off;
end
