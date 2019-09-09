
clc;
clear;

% -------------------------- PRE-PROCESSING -------------------------

ambient_T_K = 30 + 273.15;
pump_isentropic_efficiency = 0.85;
expander_isentropic_efficiency = 0.65;
economizer_effectiveness = 0.8;

working_fluid = 'R245fa'

high_T_lower_limit_K = 100 + 273.15;
high_T_upper_limit_K = 200 + 273.15;
n_high_T = 101;
high_T_list_K = linspace (high_T_lower_limit_K, high_T_upper_limit_K, n_high_T);     % Chosen high T, in K (at heater exit or expander inlet).

high_P_lower_limit_kPa = 100;
high_P_upper_limit_kPa = 10000;
n_high_P = 100;
high_P_list_kPa = linspace (high_P_lower_limit_kPa, high_P_upper_limit_kPa, n_high_P); % Chosen high P, in kPa (at pump exit or economizer inlet).

ORC_net_W_output_Jperkg = zeros (n_high_T, n_high_P);
ORC_overall_efficiency = zeros (n_high_T, n_high_P);

% -------------------------- COMPUTING OR SOLVING -------------------------

for T = 1 : n_high_T
    for P = 1 : n_high_P
        try
            [ORC_net_W_output_Jperkg(T,P), ORC_overall_efficiency(T,P)] = function_RC_input_high_T_high_P_pure_WF (expander_isentropic_efficiency, pump_isentropic_efficiency, ambient_T_K, high_T_list_K(T), high_P_list_kPa(P), working_fluid);
        catch
            ORC_net_W_output_Jperkg (T,P) = NaN;
            ORC_overall_efficiency (T,P) = NaN;
            % ORC_type (T,P) = NaN;
            % economizer (T,P) = NaN;                
        end
    end
end

% -------------------------- POST-PROCESSING -------------------------

% Finding maxima of work output and efficiency.

[max_net_W_output_Jperkg,  max_net_W_output_index] = max(ORC_net_W_output_Jperkg(:))
[max_net_W_output_T_index, max_net_W_output_P_index] = ind2sub(size(ORC_net_W_output_Jperkg), max_net_W_output_index);
max_net_W_output_T_K = high_T_list_K (max_net_W_output_T_index)
max_net_W_output_P_kPa = high_P_list_kPa (max_net_W_output_P_index)

[max_overall_efficiency,  max_overall_efficiency_index] = max(ORC_overall_efficiency(:))
[max_overall_efficiency_T_index, max_overall_efficiency_P_index] = ind2sub(size(ORC_overall_efficiency), max_overall_efficiency_index);
max_overall_efficiency_T_K = high_T_list_K (max_overall_efficiency_T_index)
max_overall_efficiency_P_kPa = high_P_list_kPa (max_overall_efficiency_P_index)

% 3-dimensional graphs.

figure ('Name', 'Net work output (J / kg) vs. expander inlet temperature vs. expander inlet pressure');
surf (high_P_list_kPa, high_T_list_K, ORC_net_W_output_Jperkg, 'EdgeColor', 'none') %, 'FaceAlpha', 0.5);
ylabel ('Expander inlet temperature (K)')
xlabel ('Pump outlet pressure (kPa)')
zlabel ('Net work output (J / kg)')
colorbar;
view (0, 90);

figure ('Name', 'Overall efficiency vs. expander inlet temperature vs. expander inlet pressure');
surf (high_P_list_kPa, high_T_list_K, ORC_overall_efficiency, 'EdgeColor', 'none') %, 'FaceAlpha', 0.5);
ylabel ('Expander inlet temperature (K)')
xlabel ('Pump outlet pressure (kPa)')
zlabel ('Overall efficiency')
colorbar;
view (0, 90);

% 345.961 s (5 minutes, 45.961 seconds) for R-245fa WF, 101 T, 100 P (10100 cases).