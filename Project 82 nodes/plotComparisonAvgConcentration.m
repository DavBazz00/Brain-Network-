function plotComparisonAvgConcentration(t_sol1, c_sol1, t_sol2, c_sol2, plot_title)
%PLOTCOMPARISONAVGINFECTIONCONCENTRATION Plots the average infection concentration 
%   from two simulations on a single figure.
%
%   plotComparisonAvgConcentration(t_sol1, c_sol1, t_sol2, c_sol2, plot_title)
%
%   INPUTS:
%       t_sol1     - Time vector for the first simulation (before dynamic variation)
%       c_sol1     - Concentration matrix for the first simulation
%       t_sol2     - Time vector for the second simulation (after dynamic variation)
%       c_sol2     - Concentration matrix for the second simulation
%       plot_title - Title for the plot
%
%   The function computes the average concentration (over all nodes) at each time point 
%   for both simulations and then plots them. The first simulation is plotted in blue and 
%   the second in orange.

    % Compute the average infection concentration for both simulations
    avg_conc1 = mean(c_sol1, 2);
    avg_conc2 = mean(c_sol2, 2);

    % Create a new figure
    figure;
    hold on;
    
    % Plot first simulation (blue)
    plot(t_sol1, avg_conc1, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]);
    
    % Plot second simulation (orange)
    plot(t_sol2, avg_conc2, 'LineWidth', 2, 'Color', [0.8500, 0.3250, 0.0980]);
    
    % Label axes and add title
    xlabel('Time');
    ylabel('Average Infection Concentration');
    title(plot_title);
    
    % Add a legend to distinguish the two curves
    legend({'Before Dynamic Variation', 'After Dynamic Variation'}, 'Location', 'best');
    
    grid on;
    hold off;
end
