function [t_baseline, p_baseline, pt_baseline, t_treatment, p_treatment, pt_treatment] = simulateHeterodimerPropagationWithTreatment(...
    A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt, num_steps, t_switch, edge_reduction)
% It Simulates two propagation scenarios for the heterodimer model:
%   1. Baseline (untreated) propagation using HeterodimerInfection_dynamic over the full time.
%   2. A treatment scenario where:
%        - Propagation is simulated untreated from t = 0 to t_switch,
%        - Then, from t = t_switch onward, dynamic edge variation is applied
%          (using HeterodimerInfection_dynamic with an initial condition).
%
% INPUTS:
%   A              - Adjacency matrix.
%   CoordTable     - Table with node coordinates and region labels.
%   k0, k1, ktilde1, k12
%                  - Model parameters for the heterodimer infection model.
%   diffusion_coeff- Diffusion coefficient.
%   dt             - Time step size (in years).
%   num_steps      - Total number of simulation steps.
%   t_switch       - Time at which the treatment starts (in years).
%   edge_reduction - Edge weight reduction factor for the treatment phase.
%
% OUTPUTS:
%   t_baseline     - Time vector from the full baseline (untreated) simulation.
%   p_baseline     - Matrix of healthy protein concentrations (each row is a time point)
%                    from the baseline simulation.
%   pt_baseline    - Matrix of misfolded protein concentrations from the baseline simulation.
%   t_treatment    - Combined time vector for the treatment scenario.
%   p_treatment    - Combined healthy protein concentration matrix for the treatment scenario.
%   pt_treatment   - Combined misfolded protein concentration matrix for the treatment scenario.
%
% Note: The switch time must be less than the total simulation time.

    % Determine the simulation step corresponding to t_switch (must be < num_steps)
    n_switch = round(t_switch / dt);
    if n_switch >= num_steps
        error('Switch time must be less than the total simulation time.');
    end
    num_steps_post = num_steps - n_switch;
    
    % --- Baseline Simulation (untreated) ---
    % Run the full simulation with no treatment (edge_reduction set to 1)
    [t_baseline, p_baseline, pt_baseline] = HeterodimerInfection_dynamic(...
        A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, 1, dt, num_steps);
    
    % --- Pre-Treatment Simulation: from t = 0 to t = t_switch ---
    [t_pre, p_pre, pt_pre] = HeterodimerInfection_dynamic(...
        A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, 1, dt, n_switch);
    
    % --- Post-Treatment Simulation: from t = t_switch to t_end ---
    % Use the final state of the pre-treatment phase as the initial condition.
    % Convert the last state (a row vector) into a column vector.
    p0_init = p_pre(end, :)';
    pt0_init = pt_pre(end, :)';
    [t_post, p_post, pt_post] = HeterodimerInfection_dynamic(...
        A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt, num_steps_post, p0_init, pt0_init);
    
    % Adjust the post-treatment time vector so that it starts at t_switch.
    t_post_adjusted = t_post + t_pre(end);
    
    % --- Combine the Pre- and Post-Treatment Segments ---
    % Remove the duplicate time point at t_switch for continuity.
    t_treatment = [t_pre, t_post_adjusted(2:end)];
    p_treatment = [p_pre; p_post(2:end, :)];
    pt_treatment = [pt_pre; pt_post(2:end, :)];
end
