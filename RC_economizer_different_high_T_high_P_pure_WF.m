
clc;
clear;

% ----------------------- PRE-PROCESSING -----------------------

ambient_T_K = 30 + 273.15;
pump_isentropic_efficiency = 0.85;
expander_isentropic_efficiency = 0.7;
economizer_effectiveness = 0.8;

WF_list = function_RC_pure_WF_list ();
n_WF = numel (WF_list)

high_T_lower_limit_K = 100 + 273.15;
high_T_upper_limit_K = 200 + 273.15;
n_high_T = 51;
high_T_list_K = linspace (high_T_lower_limit_K, high_T_upper_limit_K, n_high_T);     % Chosen high T, in K (at heater exit or expander inlet).

high_P_lower_limit_kPa = 100;
high_P_upper_limit_kPa = 10000;
n_high_P = 100;
high_P_list_kPa = linspace (high_P_lower_limit_kPa, high_P_upper_limit_kPa, n_high_P); % Chosen high P, in kPa (at pump exit or economizer inlet).

net_W_output_matrix_Jperkg = zeros (n_high_T, n_high_P, n_WF);
overall_efficiency_matrix = zeros (n_high_T, n_high_P, n_WF);

current_WF_n = 0; 

% ----------------------- COMPUTING OR SOLVING -----------------------

for WF = 1 : n_WF
    current_WF_n = current_WF_n + 1
    WF_list{WF}
    for T = 1 : n_high_T
        for P = 1 : n_high_P
            try
                [net_W_output_matrix_Jperkg(T,P,WF), overall_efficiency_matrix(T,P,WF)] = function_RC_economizer_input_high_T_high_P_pure_WF (expander_isentropic_efficiency, pump_isentropic_efficiency, ambient_T_K, economizer_effectiveness, high_T_list_K(T), high_P_list_kPa(P), WF_list{WF});
            catch
                net_W_output_matrix_Jperkg (T,P,WF) = NaN;
                overall_efficiency_matrix (T,P,WF) = NaN;
                % ORC_type (T,P,WF) = NaN;
                % economizer (T,P,WF) = NaN;                
            end
        end
    end
end

% ----------------------- POST-PROCESSING -----------------------

% Finding global maxima (for all working fluids considered together) of net work output and overall efficiency.

[max_net_W_output_Jperkg,  max_net_W_output_index] = max(net_W_output_matrix_Jperkg(:))
[max_net_W_output_T_index, max_net_W_output_P_index, max_net_W_output_WF_index] = ind2sub(size(net_W_output_matrix_Jperkg), max_net_W_output_index);
max_net_W_output_T_K = high_T_list_K (max_net_W_output_T_index)
max_net_W_output_P_kPa = high_P_list_kPa (max_net_W_output_P_index)
max_net_W_output_WF = WF_list (max_net_W_output_WF_index)

[max_overall_efficiency,  max_overall_efficiency_index] = max(overall_efficiency_matrix(:))
[max_overall_efficiency_T_index, max_overall_efficiency_P_index, max_overall_efficiency_WF_index] = ind2sub(size(overall_efficiency_matrix), max_overall_efficiency_index);
max_overall_efficiency_T_K = high_T_list_K (max_overall_efficiency_T_index)
max_overall_efficiency_P_kPa = high_P_list_kPa (max_overall_efficiency_P_index)
max_overall_efficiency_WF = WF_list (max_overall_efficiency_WF_index)

disp ('------------------------------------------------------------------')

% Finding local maxima (for each working fluid considered separately) of net work output and overall efficiency.

highest_net_W_output_list_Jperkg = zeros (1,n_WF);
highest_net_W_output_T_list_K = zeros (1,n_WF);
highest_net_W_output_P_list_kPa = zeros (1,n_WF);

highest_overall_efficiency_list = zeros (1,n_WF);
highest_overall_efficiency_T_list_K = zeros (1,n_WF);
highest_overall_efficiency_P_list_kPa = zeros (1,n_WF);

temp = zeros (n_high_T, n_high_P);
final_matrix = cell (n_WF,7);

for WF = 1 : n_WF
    
    final_matrix {WF,1} = WF_list {WF};
    
    temp = net_W_output_matrix_Jperkg(:,:,WF);
    [highest_net_W_output_list_Jperkg(1,WF),  highest_net_W_output_index] = max(temp(:));
    final_matrix {WF,2} = highest_net_W_output_list_Jperkg(1,WF);
    [highest_net_W_output_T_index, highest_net_W_output_P_index] = ind2sub(size(temp), highest_net_W_output_index);
    highest_net_W_output_T_list_K(1,WF) = high_T_list_K (highest_net_W_output_T_index);
    final_matrix {WF,3} = highest_net_W_output_T_list_K(1,WF);
    highest_net_W_output_P_list_kPa(1,WF) = high_P_list_kPa (highest_net_W_output_P_index);
    final_matrix {WF,4} = highest_net_W_output_P_list_kPa(1,WF);

    temp = overall_efficiency_matrix(:,:,WF);
    [highest_overall_efficiency_list(1,WF),  highest_overall_efficiency_index] = max(temp(:));
    final_matrix {WF,5} = highest_overall_efficiency_list(1,WF);
    [highest_overall_efficiency_T_index, highest_overall_efficiency_P_index] = ind2sub(size(temp), highest_overall_efficiency_index);
    highest_overall_efficiency_T_list_K(1,WF) = high_T_list_K (highest_overall_efficiency_T_index);
    final_matrix {WF,6} = highest_overall_efficiency_T_list_K(1,WF);
    highest_overall_efficiency_P_list_kPa(1,WF) = high_P_list_kPa (highest_overall_efficiency_P_index);
    final_matrix {WF,7} = highest_overall_efficiency_P_list_kPa(1,WF);
    
end

disp ('Working fluid,    max W output (J/kg),    T (K),    P (kPa),    max efficiency,    T (K),    P (kPa)')
disp (final_matrix)

% 3 dimensional graphs.

surface_color = zeros (3, n_WF);
for WF = 1 : n_WF
    for c = 1 : 3
        surface_color (c, WF) = rand;
    end
end

figure ('Name', 'Net work output (J / kg) vs. expander inlet temperature vs. expander inlet pressure');
for WF = 1 : n_WF
    surf (high_P_list_kPa, high_T_list_K, net_W_output_matrix_Jperkg(:,:,WF), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', surface_color (:, WF));
    hold on;
end
ylabel ('Expander inlet temperature (K)')
xlabel ('High pressure (kPa)')
zlabel ('Net work output (J / kg)')
legend (WF_list)

figure ('Name', 'Overall efficiency vs. expander inlet temperature vs. expander inlet pressure');
for WF = 1 : n_WF
    surf (high_P_list_kPa, high_T_list_K, overall_efficiency_matrix(:,:,WF), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', surface_color (:, WF));
    hold on;
end
ylabel ('Expander inlet temperature (K)')
xlabel ('High pressure (kPa)')
zlabel ('Overall efficiency')
legend (WF_list)

% Total run time 18602.028 seconds (5.16723 hours) on Intel Core i5 6500 system with RAM 24 GB and Windows 10 OS (IISc CST lab computer).