clc;
clear;

load workspace_RC_economizer_high_T_100_500_C_high_P_1_100_bar_pure_WF_73.mat

selected_WF_list = {"WATER", "R245FA", "TOLUENE", "ETHANOL", "ACETONE", "BENZENE", "MLINOLEN"};
n_selected_WF = numel (selected_WF_list);

selected_WF_index_list = zeros (1, n_selected_WF);

selected_WF_net_W_output_matrix_Jperkg = zeros (n_high_T, n_high_P, n_selected_WF);
selected_WF_efficiency_matrix = zeros (n_high_T, n_high_P, n_selected_WF);

for WF_select = 1 : n_selected_WF
    for WF_all = 1 : n_WF
        if WF_list{WF_all} == selected_WF_list{WF_select}
            selected_WF_net_W_output_matrix_Jperkg (:,:,WF_select) = net_W_output_matrix_Jperkg (:,:,WF_all);
            selected_WF_efficiency_matrix (:,:,WF_select) = overall_efficiency_matrix (:,:,WF_all);
        end
    end
end

% color = zeros (3, n_selected_WF);
% for WF_i = 1 : n_WF
%     for c = 1 : 3
%         color (c, WF_i) = rand;
%     end
% end

color_list = {'blue', 'red', 'green', 'yellow', 'magenta', 'cyan', 'black'};

% -------------- 3 dimensional graphs ----------------------

figure ('Name', 'Net work output (J / kg) vs. expander inlet temperature vs. expander inlet pressure');
for WF_i = 1 : n_selected_WF
    surf (high_P_list_kPa, high_T_list_K, selected_WF_net_W_output_matrix_Jperkg(:,:,WF_i), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', color_list {WF_i});
    hold on;
end
xlabel ('Pump outlet pressure (kPa)')
ylabel ('Expander inlet temperature (K)')
zlabel ('Net work output (J / kg)')
legend (selected_WF_list)

figure ('Name', 'Overall efficiency vs. expander inlet temperature vs. expander inlet pressure');
for WF_i = 1 : n_selected_WF
    surf (high_P_list_kPa, high_T_list_K, selected_WF_efficiency_matrix(:,:,WF_i), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', color_list {WF_i});
    hold on;
end
xlabel ('Pump outlet pressure (kPa)')
ylabel ('Expander inlet temperature (K)')
zlabel ('Overall efficiency')
legend (selected_WF_list)

% -------------- 2 dimensional graphs ----------------------

highest_T_index = numel (high_T_list_K);
% high_T_list_K (highest_T_index)

figure ('Name', 'Net work output (J / kg) vs. expander inlet pressure, at highest expander inlet temperature.');
for WF_i = 1 : n_selected_WF
    plot (high_P_list_kPa, selected_WF_net_W_output_matrix_Jperkg(highest_T_index,:,WF_i), 'Color', color_list {WF_i}, 'LineWidth', 1.5)
    hold on;
end
xlabel ('Pump outlet pressure (kPa)')
ylabel ('Net work output (J / kg)')
legend (selected_WF_list)

figure ('Name', 'Overall efficiency vs. expander inlet pressure, at highest expander inlet temperature.');
for WF_i = 1 : n_selected_WF
    plot (high_P_list_kPa, selected_WF_efficiency_matrix(highest_T_index,:,WF_i), 'Color', color_list {WF_i}, 'LineWidth', 1.5)
    hold on;
end
xlabel ('Pump outlet pressure (kPa)')
ylabel ('Overall efficiency')
legend (selected_WF_list)
