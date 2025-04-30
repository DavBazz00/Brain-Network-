% simulateHeterodimerPropagationWithTreatment.m
function [t_baseline, p_baseline, pt_baseline, t_treatment, p_treatment, pt_treatment] = ...
    simulateHeterodimerPropagationWithTreatment(...
        A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt1, dt2, num_steps, t_switch, edge_reduction)
% SIMULATEHETERODIMERPROPAGATIONWITHTREATMENT
%   1. Baseline (no treatment) con HeterodimerInfection
%   2. Pre-trattamento fino a t_switch con HeterodimerInfection
%   3. Post-trattamento da t_switch a fine con HeterodimerInfection_dynamic

    % numero di passi per baseline e pre-trat.
    n_steps_baseline = round(num_steps / dt1);
    n_switch     = round(t_switch / dt2);
    if t_switch >= num_steps
        error('Switch time must be less than total simulation time.');
    end
    % passi post-tratt.
    num_steps_post = num_steps - n_switch;

    % --- 1) Baseline ---
    [t_baseline, p_baseline, pt_baseline] = ...
        HeterodimerInfection(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt1, n_steps_baseline);

    % --- 2) Pre-Trattamento ---
    [t_pre, p_pre, pt_pre] = ...
        HeterodimerInfection(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt1, n_steps_pre);

    % --- 3) Post-Trattamento ---
    p_init  = p_pre(end-dt1, :)';
    pt_init = pt_pre(end-dt1, :)';
    [t_post, p_post, pt_post] = ...
        HeterodimerInfection_treatment_dynamic(...
            A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt2, num_steps_post, p_init, pt_init);

    % riallineamento asse temporale
    t_post_shifted = t_post + t_pre(end);

    % concatenazione (eliminando punto doppio a t_switch)
    t_treatment = [t_pre,    t_post_shifted(2:end)];
    p_treatment = [p_pre;   p_post(2:end, :)];
    pt_treatment= [pt_pre;  pt_post(2:end, :)];
end
