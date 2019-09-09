clc
clear
 
load workspace_RC_economizer_different_high_T_(100_200_C)_high_P_(1_100_bar)_3_C_mixture_WF;

figure ('Name', 'Net work output (J / kg) vs. pump outlet pressure, at expander inlet temperature of 200 ?C');
for WF = 1 : n_WF
    plot (high_P_list_kPa, net_W_output_matrix_Jperkg (51,:,WF), 'Color', surface_18_colors (:, WF+1), 'LineWidth', 1)
    hold on
end
xlabel ('Pump outlet pressure (kPa)')
ylabel ('Net work output (J / kg of working fluid)')
legend_object = legend (mass_composition_array);
title (legend_object, legend_title);

figure ('Name', 'Overall efficiency vs. pump outlet pressure, at expander inlet temperature of 200 ?C');
for WF = 1 : n_WF
    plot (high_P_list_kPa, overall_efficiency_matrix (51,:,WF), 'Color', surface_18_colors (:, WF+1), 'LineWidth', 1)
    hold on
end
xlabel ('Pump outlet pressure (kPa)')
ylabel ('Overall efficiency')
legend_object = legend (mass_composition_array);
title (legend_object, legend_title);