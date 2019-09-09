clc
clear

load workspace_RC_WF_R245fa_high_T_100_200_C_high_P_1_100_bar
load workspace_RC_economizer_WF_R245fa_high_T_100_200_C_high_P_1_100_bar

% ------------------------- 3 dimensional graphs -------------------------

figure ('Name', 'Net work output (J / kg of working fluid) vs. expander inlet temperature vs. Pump outlet pressure');
surf (high_P_list_kPa, high_T_list_K, ORC_net_W_output_Jperkg, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', 'red');
hold on;
surf (high_P_list_kPa, high_T_list_K, ORC_economizer_net_W_output_Jperkg, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', 'green');
ylabel ('Expander inlet temperature (K)')
xlabel ('Pump outlet pressure (kPa)')
zlabel ('Net work output (J / kg of working fluid)')
legend ({'Without economizer', 'With economizer'});

figure ('Name', 'Overall efficiency vs. expander inlet temperature vs. Pump outlet pressure');
surf (high_P_list_kPa, high_T_list_K, ORC_overall_efficiency, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', 'red');
hold on;
surf (high_P_list_kPa, high_T_list_K, ORC_economizer_overall_efficiency, 'EdgeColor', 'none', 'FaceAlpha', 0.5, 'FaceColor', 'green');
ylabel ('Expander inlet temperature (K)')
xlabel ('Pump outlet pressure (kPa)')
zlabel ('Overall efficiency')
legend ({'Without economizer', 'With economizer'});

% ------------------------- 2 dimensional graphs -------------------------

figure ('Name', 'Net work output (J / kg) vs. pump outlet pressure, at expander inlet temperature of 200 ?C');
plot (high_P_list_kPa, ORC_net_W_output_Jperkg (101,:), 'Color', 'red', 'LineWidth', 2)
hold on;
plot (high_P_list_kPa, ORC_economizer_net_W_output_Jperkg (101,:), 'Color', 'green', 'LineWidth', 2)
xlabel ('Pump outlet pressure (kPa)')
ylabel ('Net work output (J / kg of working fluid)')
legend ({'Without economizer', 'With economizer'}, 'Location', 'southeast');

figure ('Name', 'Overall efficiency vs. pump outlet pressure, at expander inlet temperature of 200 ?C');
plot (high_P_list_kPa, ORC_overall_efficiency (101,:), 'Color', 'red', 'LineWidth', 2)
hold on;
plot (high_P_list_kPa, ORC_economizer_overall_efficiency (101,:), 'Color', 'green', 'LineWidth', 2)
xlabel ('Pump outlet pressure (kPa)')
ylabel ('Overall efficiency')
legend ({'Without economizer', 'With economizer'}, 'Location', 'southeast');
