clc;
clear;

substance1 = 'pentane'
substance2 = 'propane'
substance3 = 'butane'
substance4 = 'hexane'
% substance4 = 'toluene'
mass_composition = [0.3, 0.35, 0.15, 0.2] % Mass fractions of substances.
n_mass_composition = numel (mass_composition);
sum = 0;
for a = 1 : n_mass_composition
    sum = sum + mass_composition (a);
end
if sum ~= 1
    error ('Sum of elements in mass_composition vector should equal 1')
end

expander_isentropic_efficiency = 0.65;
if (expander_isentropic_efficiency < 0) || (expander_isentropic_efficiency > 1)
    error ('expander_isentropic_efficiency must be in range [0, 1].')
end
pump_isentropic_efficiency = 0.85;
if (pump_isentropic_efficiency < 0) || (pump_isentropic_efficiency > 1)
    error ('pump_isentropic_efficiency must be in range [0, 1].')
end
economizer_effectiveness = 0.8;
if (economizer_effectiveness < 0) || (economizer_effectiveness > 1)
    error ('economizer_effectiveness must be in range [0, 1].')
end

ambient_T_K = 30 + 273.15;

high_T_lower_limit_K = 100 + 273.15;
high_T_upper_limit_K = 500 + 273.15;
n_high_T = 15;
high_T_K = linspace (high_T_lower_limit_K, high_T_upper_limit_K, n_high_T);     % Chosen high T, in K (at heater exit or expander inlet).

high_P_lower_limit_kPa = 200;
high_P_upper_limit_kPa = 10000;
n_high_P = 15;
high_P_kPa = linspace (high_P_lower_limit_kPa, high_P_upper_limit_kPa, n_high_P); % Chosen high P, in kPa (at pump exit or economizer inlet).

net_W_output_Jperkg = zeros (n_high_T, n_high_P);
overall_efficiency = zeros (n_high_T, n_high_P);
% ORC_type = zeros (n_high_T, n_high_P);
% economizer = zeros (n_high_T, n_high_P);

for T = 1 : n_high_T
    for P = 1 : n_high_P
        try
            [net_W_output_Jperkg(T,P), overall_efficiency(T,P)] = function_RC_economizer_input_high_T_high_P_4_C_mixture_WF (expander_isentropic_efficiency, pump_isentropic_efficiency, ambient_T_K, economizer_effectiveness, high_T_K(T), high_P_kPa(P), mass_composition, substance1, substance2, substance3, substance4);
        catch
            net_W_output_Jperkg (T,P) = NaN;
            overall_efficiency (T,P) = NaN;
            % ORC_type (T,P) = NaN;
            % economizer (T,P) = NaN;
        end
    end
end

% 3 dimensional graphs

figure ('Name', 'Net work output vs. expander inlet temperature vs. high pressure.');
surf (high_P_kPa, high_T_K, net_W_output_Jperkg);
ylabel ('Expander inlet temperature (K)')
xlabel ('Pump outlet pressure (kPa)')
zlabel ('Net work output (J / kg)')
colorbar

figure ('Name', 'Overall efficiency vs. expander inlet temperature vs. high pressure.');
surf (high_P_kPa, high_T_K, overall_efficiency);
ylabel ('Expander inlet temperature (K)')
xlabel ('Pump outlet pressure (kPa)')
zlabel ('Overall efficiency')
colorbar

% Displaying maxima of net W output, overall efficiciency, and their T, P.

[max_net_W_output_Jperkg, max_net_W_output_index] = max (net_W_output_Jperkg(:))
[max_net_W_output_row, max_net_W_output_column] = ind2sub (size (net_W_output_Jperkg), max_net_W_output_index);
max_net_W_output_T_K = high_T_K (max_net_W_output_row)
max_net_W_output_P_kPa = high_P_kPa (max_net_W_output_column)

[max_overall_efficiency, max_overall_efficiency_index] = max (overall_efficiency(:))
[max_overall_efficiency_row, max_overall_efficiency_column] = ind2sub (size (overall_efficiency), max_overall_efficiency_index);
max_overall_efficiency_T_K = high_T_K (max_overall_efficiency_row)
max_overall_efficiency_P_kPa = high_P_kPa (max_overall_efficiency_column)
