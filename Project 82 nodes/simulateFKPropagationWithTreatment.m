function [t_baseline, c_baseline, t_treatment, c_treatment] = simulateFKPropagationWithTreatment(...
    A, CoordTable, diffusion, a, dt1, dt2, num_steps, t_switch, edge_reduction)
% SIMULATEFKPROPAGATIONWITHTREATMENT
%   1) Baseline (untreated) con FK_propagation
%   2) Pre-trattamento fino a t_switch con FK_propagation
%   3) Post-trattamento da t_switch a fine con FK_propagation_treatment_dynamic
%
% INPUT:
%   A             – matrice di adiacenza
%   CoordTable    – tabella con coordinate e regioni (colonna 4 = regione)
%   diffusion     – coefficiente di diffusione
%   a             – tasso di crescita logistica
%   dt1           – passo temporale per baseline e pre-trattamento
%   dt2           – passo temporale per post-trattamento
%   num_steps     – tempo totale simulazione (in anni)
%   t_switch      – tempo di inizio trattamento (in anni)
%   edge_reduction– fattore di riduzione pesi per FK_propagation_treatment_dynamic
%
% OUTPUT:
%   t_baseline    – vettore temporale per la simulazione baseline
%   c_baseline    – matrice concentrazioni per la simulazione baseline
%   t_treatment   – vettore temporale combinato (pre + post treatment)
%   c_treatment   – matrice concentrazioni combinata (pre + post treatment)

    %--- calcolo numero passi ---
    % passi per la fase baseline (tutto a dt1)
    n_steps_baseline = round(num_steps / dt1);
    % passi di pre-trattamento (dt1)
    n_switch = round(t_switch / dt2);
    if n_switch >= n_steps_baseline
        error('t_switch deve essere minore del tempo totale num_steps.');
    end
    % passi post-trattamento (dt2)
    % remaining_time = num_steps - n_switch * dt1;
    num_steps_post = num_steps - n_switch;

    %--- 1) Baseline Simulation ---
    [t_baseline, c_baseline] = FK_propagation(...
        A, CoordTable, diffusion, a, dt1, n_steps_baseline);

    %--- 2) Pre-Treatment Simulation ---
    [t_pre, c_pre] = FK_propagation(...
        A, CoordTable, diffusion, a, dt1, n_switch);

    %--- 3) Post-Treatment Simulation ---
    % condizione iniziale = ultimo stato di pre-trattamento
    c_init = c_pre(end, :)';
    [t_post, c_post] = FK_propagation_treatment_dynamic(...
        A, CoordTable, diffusion, a, dt2, num_steps_post, edge_reduction, c_init);

    % riallineo l’asse dei tempi in modo che t_post(1) = t_pre(end)
    t_post = t_post + t_pre(end);

    %--- concatenazione (tolgo il punto duplicato a t_switch) ---
    t_treatment = [t_pre,       t_post(2:end)];
    c_treatment = [c_pre;       c_post(2:end, :)];
end
