function plotFourScenariosFK(A, CoordTable, diffusion, a, edge_reduction, dt_aging, dt_treat, t_switch, t_end)
% plotFourScenariosFK  Compute and plot four FK scenarios:
%   baseline, treatment-only, aging+treatment, and pure-aging restart.
%
% Usage:
%   plotFourScenariosFK(A, CoordTable, diffusion, a, edge_reduction, ...
%       dt_aging, dt_treat, t_switch, t_end)

% Number of steps for baseline
num_steps_base = round(t_end / dt_aging);

% 1) Static baseline
[t_bs, c_bs] = FK_propagation(A, CoordTable, diffusion, a, dt_aging, num_steps_base);
avg_bs = mean(c_bs, 2);

% 2) Treatment only
[t_tr, c_tr] = FK_propagation_treatment(A, CoordTable, diffusion, a, edge_reduction, dt_aging, dt_treat, t_switch, t_end);
avg_tr = mean(c_tr, 2);

% 3) Aging + Treatment
[t_ct, c_ct] = FK_propagation_treatment_aging(A, CoordTable, diffusion, a, edge_reduction, dt_aging, dt_treat, t_switch, t_end);
avg_ct = mean(c_ct, 2);

% 4) Pure-aging restart from t_switch
idx_sw = find(t_ct >= t_switch, 1);
c_init = c_ct(idx_sw, :)';
num_steps_pa = round((t_end - t_switch) / dt_aging);
[t_pa2, c_pa2] = FK_propagation_dynamic(A, CoordTable, diffusion, a, dt_aging, num_steps_pa, edge_reduction, c_init);
t_pa = t_switch + t_pa2;
avg_pa = mean(c_pa2, 2);

% Prepare piecewise segments for treatment-only and aging+treatment
idx_tr = find(t_tr >= t_switch, 1);
t_tr_pre = t_tr(1:idx_tr);      avg_tr_pre  = avg_tr(1:idx_tr);
t_tr_post = t_tr(idx_tr:end);   avg_tr_post = avg_tr(idx_tr:end);

idx_ct = find(t_ct >= t_switch, 1);
t_ct_pre = t_ct(1:idx_ct);      avg_ct_pre  = avg_ct(1:idx_ct);
t_ct_post = t_ct(idx_ct:end);   avg_ct_post = avg_ct(idx_ct:end);

% Plot all four with piecewise colors
figure; hold on;
  % Baseline (blue)
  h1 = plot(t_bs,    avg_bs,    '-b','LineWidth',2);
  % Pure-aging restart (green)
  h2 = plot(t_pa,    avg_pa,    '-g','LineWidth',2);
  % Treatment-only: blue before, red after
  plot(t_tr_pre,     avg_tr_pre,'-b','LineWidth',2);
  h3 = plot(t_tr_post, avg_tr_post,'-r','LineWidth',2);
  % Aging+Treatment: green before, magenta after
  plot(t_ct_pre,     avg_ct_pre,'-g','LineWidth',2);
  h4 = plot(t_ct_post, avg_ct_post,'-m','LineWidth',2);
  % switch marker
  xline(t_switch,   '--k','LineWidth',2);
hold off;

xlabel('Time (years)');
ylabel('⟨infection concentration⟩');
title('Comparison of Four FK Scenarios');
legend([h1,h2,h3,h4], 'Baseline','Aging','Treatment only','Aging+Treatment','Location','northwest');

xlim([0, t_end]);
grid on;
end
