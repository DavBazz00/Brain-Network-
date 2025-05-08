function plotFourScenariosHeterodimer(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt_aging, dt_treat, t_switch, t_end)
% plotFourScenariosHeterodimer  Compute and plot four heterodimer scenarios:
%   baseline, treatment-only, aging+treatment, and pure-aging restart.
%
% Usage:
%   plotFourScenariosHeterodimer(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, ...
%       edge_reduction, dt_aging, dt_treat, t_switch, t_end)

% Number of steps for baseline
num_steps_base = round(t_end / dt_aging);

% 11.1 Static baseline
[t_bs, ~, pt_bs] = HeterodimerInfection(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, dt_aging, num_steps_base);
avg_bs = mean(pt_bs, 2);

% 11.3 Treatment only
[t_tr, ~, pt_tr] = HeterodimerInfection_treatment(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt_aging, dt_treat, t_switch, t_end);
avg_tr = mean(pt_tr, 2);

% 11.4 Aging + Treatment
[t_ct, ~, pt_ct] = HeterodimerInfection_treatment_aging(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt_aging, dt_treat, t_switch, t_end);
avg_ct = mean(pt_ct, 2);

% 11.5 Pure-aging restart from t_switch
idx_sw   = find(t_ct >= t_switch, 1);
p_init   = pt_ct(idx_sw, :)';
pt_init  = pt_ct(idx_sw, :)';  % seed misfolded restart
dt_pa    = dt_aging;
num_steps_pa = round((t_end - t_switch) / dt_pa);
[t_pa2, ~, pt_pa2] = HeterodimerInfection_dynamic(A, CoordTable, k0, k1, ktilde1, k12, diffusion_coeff, edge_reduction, dt_pa, num_steps_pa, [], pt_init);
t_pa     = t_switch + t_pa2;
avg_pa   = mean(pt_pa2, 2);

% 11.6 Prepare piecewise segments
idx_tr   = find(t_tr >= t_switch, 1);
t_tr_pre = t_tr(1:idx_tr);      avg_tr_pre  = avg_tr(1:idx_tr);
t_tr_post= t_tr(idx_tr:end);    avg_tr_post = avg_tr(idx_tr:end);

idx_ct   = find(t_ct >= t_switch, 1);
t_ct_pre = t_ct(1:idx_ct);      avg_ct_pre  = avg_ct(1:idx_ct);
t_ct_post= t_ct(idx_ct:end);    avg_ct_post = avg_ct(idx_ct:end);

% 11.7 Plot
figure; hold on;
  % Baseline (blue)
  h1 = plot(t_bs,      avg_bs,     '-b', 'LineWidth',2);
  % Pure-aging restart (green)
  h2 = plot(t_pa,      avg_pa,     '-g', 'LineWidth',2);
  % Treatment-only: blue pre-switch, red post-switch
  plot(t_tr_pre,       avg_tr_pre, '-b', 'LineWidth',2);
  h3 = plot(t_tr_post,  avg_tr_post,'-r', 'LineWidth',2);
  % Aging+Treatment: green pre-switch, magenta post-switch
  plot(t_ct_pre,       avg_ct_pre, '-g', 'LineWidth',2);
  h4 = plot(t_ct_post,  avg_ct_post,'-m', 'LineWidth',2);
  % switch marker
  xline(t_switch,     '--k','LineWidth',2);
hold off;

xlabel('Time (years)');
ylabel('⟨misfolded protein⟩');
title('Comparison of Four Heterodimer Scenarios');
legend([h1,h2,h3,h4], 'Baseline','Aging','Treatment only','Aging+Treatment','Location','northwest');

xlim([0, t_end]);
grid on;
end
