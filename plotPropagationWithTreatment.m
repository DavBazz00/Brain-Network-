function plotPropagationWithTreatment(t_baseline, c_baseline, t_treatment, c_treatment, t_switch, plot_title)
% PLOTPROPAGATIONWITHTREATMENT Plots the average infection concentration 
% for both the baseline (untreated) simulation and the treatment scenario.
%
%   - The baseline evolution is plotted in blue for all time values.
%   - The treatment evolution is plotted in orange only for times >= t_switch.
%   - A vertical dashed black line is added at t = t_switch.
%
% INPUTS:
%   t_baseline  - Time vector from the full baseline simulation.
%   c_baseline  - Concentration matrix from the baseline simulation.
%   t_treatment - Time vector from the treatment simulation (includes 
%                 both pre- and post-switch segments).
%   c_treatment - Concentration matrix from the treatment simulation.
%   t_switch    - The switching time (in years) when treatment starts.
%   plot_title  - Title of the plot.
%
% The function computes the average concentration at each time step and
% plots the two scenarios accordingly.

    % Compute average concentration (mean over nodes) for baseline
    avg_baseline = mean(c_baseline, 2);
    
    % Compute average concentration for treatment
    avg_treatment = mean(c_treatment, 2);
    
    % Find the first index in t_treatment where time >= t_switch.
    idx = find(t_treatment >= t_switch, 1, 'first');
    
    % Extract the treatment data from t_switch onward.
    treatment_time = t_treatment(idx:end);
    treatment_avg  = avg_treatment(idx:end);
    
    % Create figure and plot the full baseline evolution in blue.
    figure;
    plot(t_baseline, avg_baseline, 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
    hold on;
    
    % Plot the treatment evolution only for times >= t_switch in orange.
    plot(treatment_time, treatment_avg, 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980]);
    
    % Add a vertical dashed black line at t = t_switch.
    xline(t_switch, '--k', 'LineWidth', 2);
    
    % Add labels, title, legend, and grid.
    xlabel('Time (years)');
    ylabel('Average Infection Concentration');
    title(plot_title);
    legend('Baseline (untreated)', 'Treatment (post-switch)', 'Treatment start', 'Location', 'best');
    grid on;
    hold off;
end
