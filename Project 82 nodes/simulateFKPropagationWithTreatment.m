function [t_baseline, c_baseline, t_treatment, c_treatment] = simulateFKPropagationWithTreatment(A, CoordTable, diffusion, a, dt1,dt2, num_steps, t_switch, edge_reduction)
% SIMULATEFKPROPAGATIONWITHTREATMENT Simulates two scenarios:
%   1. Baseline (untreated) propagation using FK_propagation over the full time.
%   2. A treatment scenario where:
%        - Propagation is simulated untreated from t = 0 to t_switch,
%        - Then, from t = t_switch onward, dynamic edge variation is used
%          (using FK_propagation_dynamic with an initial condition).
%
% INPUTS:
%   A             - Adjacency matrix.
%   CoordTable    - Table with node coordinates and region labels.
%   diffusion     - Diffusion parameter.
%   a             - Logistic growth rate parameter.
%   dt            - Time step size (in years).
%   num_steps     - Total number of simulation steps.
%   t_switch      - Time at which the treatment starts (in years).
%   edge_reduction- Edge weight reduction factor for the dynamic treatment.
%
% OUTPUTS:
%   t_baseline    - Time vector for the full baseline (untreated) simulation.
%   c_baseline    - Concentration matrix from the baseline simulation.
%   t_treatment   - Combined time vector for the treatment scenario.
%   c_treatment   - Combined concentration matrix for the treatment scenario.
%
% Note: The switch time must be less than the total simulation time.

    % Determine the simulation step at which to switch (must be < num_steps)
    n_switch = round(t_switch / dt2);
    if n_switch >= num_steps
        error('Switch time must be less than the total simulation time.');
    end

    % --- Baseline Simulation (untreated) ---
    [t_baseline, c_baseline] = FK_propagation(A, CoordTable, diffusion, a, dt1, num_steps/dt1);

    % --- Pre-Treatment Simulation: from t = 0 to t = t_switch ---
    [t_pre, c_pre] = FK_propagation(A, CoordTable, diffusion, a, dt1, n_switch/dt1);

    % --- Post-Treatment Simulation: from t = t_switch to t_end ---
    % Use the last state of the pre-treatment phase as the initial condition for dynamic propagation.
    num_steps_post = num_steps - n_switch;
    c_init = c_pre(end, :)';  % Convert last state from row to column
    [t_post, c_post] = FK_propagation_treatment_dynamic(A, CoordTable, diffusion, a, dt2, num_steps_post, edge_reduction, c_init);

    % Adjust the post-treatment time vector so it starts at t_switch
    t_post_adjusted = t_post + t_pre(end);

    % --- Combine the Pre- and Post-Treatment Segments ---
    % Remove the duplicate time point at the switch for continuity.
    t_treatment = [t_pre, t_post_adjusted(2:end)];
    c_treatment = [c_pre; c_post(2:end, :)];
end
