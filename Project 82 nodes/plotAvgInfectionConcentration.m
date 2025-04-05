function plotAvgInfectionConcentration(t_sol, c_sol, plot_title)
%PLOTAVGINFECTIONCONCENTRATION Plots the average infection concentration over time.
%
%   plotAvgInfectionConcentration(t_sol, c_sol)
%
%   INPUTS:
%       t_sol - Time vector for the simulation.
%       c_sol - Matrix of infection concentrations over time.
%       plot_title - Choose title for the plot
%
%   This function computes the average infection concentration at each time
%   point (by averaging over all nodes) and creates a plot.

    % Compute the average infection concentration at each time point
    avg_conc = mean(c_sol, 2);

    % Create a new figure for the plot
    figure;
    plot(t_sol, avg_conc, 'LineWidth', 2, 'Color', [0, 0.4470, 0.7410]); % MATLAB default blue
    xlabel('Time');
    ylabel('Average Infection Concentration');
    title(plot_title);
    grid on;
    hold off;
end
