% about 14 minutes 45 seconds for 6 WF, 40 T, 40 P (9600 cases).
% 2592.253 s (43 minutes and 12.253 seconds) for 15 WF, 21 T, 34 P (10710 cases).

clc;
clear;

% PRE-PROCESSING.

expander_isentropic_efficiency = 0.65;
pump_isentropic_efficiency = 0.85;
economizer_effectiveness = 0.8;
ambient_T_K = 30 + 273;

substance1 = 'toluene'
substance2 = 'benzene'
substance3 = 'cyclohex'

step = 0.25;
n_WF = 0;
for c1 = 0 : step : 1
    for c2 = 0 : step : 1 - c1
        c3 = 1 - c1 - c2;
        n_WF = n_WF + 1;
    end
end
n_WF
mass_composition = zeros (n_WF, 3);
n = 0;
for c1 = 0 : step : 1
    for c2 = 0 : step : 1 - c1
        c3 = 1 - c1 - c2;
        n = n + 1;
        mass_composition (n,:) = [c1, c2, c3];
    end
end

high_T_lower_limit_K = 100 + 273;
high_T_upper_limit_K = 200 + 273;
n_high_T = 51;
high_T_list_K = linspace (high_T_lower_limit_K, high_T_upper_limit_K, n_high_T);     % Chosen high T, in K (at heater exit or expander inlet).

high_P_lower_limit_kPa = 100;
high_P_upper_limit_kPa = 10000;
n_high_P = 100;
high_P_list_kPa = linspace (high_P_lower_limit_kPa, high_P_upper_limit_kPa, n_high_P); % Chosen high P, in kPa (at pump exit or economizer inlet).

net_W_output_matrix_Jperkg = zeros (n_high_T, n_high_P, n_WF);
overall_efficiency_matrix = zeros (n_high_T, n_high_P, n_WF);
% ORC_type = zeros (n_high_T, n_high_P, n_WF);
% economizer = zeros (n_high_T, n_high_P, n_WF);

% COMPUTING.

current_WF_n = 0;

for WF = 1 : n_WF
    current_WF_n = current_WF_n + 1
    mass_composition(WF,:)
    for T = 1 : n_high_T
        for P = 1 : n_high_P
            % try
                [net_W_output_matrix_Jperkg(T,P,WF), overall_efficiency_matrix(T,P,WF)] = function_RC_economizer_input_high_T_high_P_3_C_mixture_WF (expander_isentropic_efficiency, pump_isentropic_efficiency, ambient_T_K, economizer_effectiveness, high_T_list_K(T), high_P_list_kPa(P), mass_composition(WF,:), substance1, substance2, substance3);
            % catch
                % net_W_output_Jperkg (T,P,WF) = NaN;
                % overall_efficiency (T,P,WF) = NaN;
                % ORC_type (T,P,WF) = NaN;
                % economizer (T,P,WF) = NaN;                
            % end
        end
    end
end

% POST-PROCESSING.

mass_composition_array = cell (n_WF, 1);
for WF = 1 : n_WF
    mass_composition_array{WF,1} = num2str (mass_composition(WF,:));
end

% Displaying global maxima (considering all woking fluids together) of work output and efficicency.

[max_net_W_output_Jperkg, max_net_W_output_index] = max (net_W_output_matrix_Jperkg(:))
[max_net_W_output_row, max_net_W_output_column, max_net_W_output_page] = ind2sub (size (net_W_output_matrix_Jperkg), max_net_W_output_index);
max_net_W_output_T_K = high_T_list_K (max_net_W_output_row)
max_net_W_output_P_kPa = high_P_list_kPa (max_net_W_output_column)
max_net_W_output_mass_composition = mass_composition (max_net_W_output_page,:)

[max_overall_efficiency, max_overall_efficiency_index] = max (overall_efficiency_matrix(:))
[max_overall_efficiency_row, max_overall_efficiency_column, max_overall_efficiency_page] = ind2sub (size (overall_efficiency_matrix), max_overall_efficiency_index);
max_overall_efficiency_T_K = high_T_list_K (max_overall_efficiency_row)
max_overall_efficiency_P_kPa = high_P_list_kPa (max_overall_efficiency_column)
max_overall_efficiency_mass_composition = mass_composition (max_overall_efficiency_page,:)

disp ('------------------------------------------------------------------')

% Displaying local maxima (considering each woking fluid seperately) of work output and efficicency.

highest_net_W_output_list_Jperkg = zeros (1,n_WF);
highest_net_W_output_T_list_K = zeros (1,n_WF);
highest_net_W_output_P_list_kPa = zeros (1,n_WF);

highest_overall_efficiency_list = zeros (1,n_WF);
highest_overall_efficiency_T_list_K = zeros (1,n_WF);
highest_overall_efficiency_P_list_kPa = zeros (1,n_WF);

temp = zeros (n_high_T, n_high_P);
final_matrix = cell (n_WF,7);

for WF = 1 : n_WF
    
    final_matrix {WF,1} = mass_composition_array {WF,1};
    
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

% -------------- 3 dimensional graphs ----------------------

legend_title = cat (2, substance1, ' : ', substance2, ' : ', substance3);

% surface_color = zeros (3, n_WF);
% for WF = 1 : n_WF
%     for c = 1 : 3
%         surface_color (c, WF) = rand;
%     end
% end

surface_18_colors = zeros (3, 18);
column = 0;
for c1 = [0, 1]
    for c2 = [0, 0.5, 1]
        for c3 = [0, 0.5, 1]
            column = column + 1;
            surface_18_colors (1, column) = c1;
            surface_18_colors (2, column) = c2;
            surface_18_colors (3, column) = c3;
        end
    end
end
% Color 1 is black and color 18 is white. These may be avoided for clarity in figures.

figure ('Name', 'Net work output (J / kg) vs. expander inlet temperature vs. expander inlet pressure');
for WF = 1 : n_WF
    surf (high_P_list_kPa, high_T_list_K, net_W_output_matrix_Jperkg(:,:,WF), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', surface_18_colors (:, WF+1));
    % WF+1 to avoid the first color (black).
    hold on;
end
ylabel ('Expander inlet temperature (K)')
xlabel ('Pump outlet pressure (kPa)')
zlabel ('Net work output (J / kg)')
legend_object = legend (mass_composition_array);
title (legend_object, legend_title);

figure ('Name', 'Overall efficiency vs. expander inlet temperature vs. expander inlet pressure');
for WF = 1 : n_WF
    surf (high_P_list_kPa, high_T_list_K, overall_efficiency_matrix(:,:,WF), 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'FaceColor', surface_18_colors (:, WF+1));
    % WF+1 to avoid the first color (black).
    hold on;
end
ylabel ('Expander inlet temperature (K)')
xlabel ('Pump outlet pressure (kPa)')
zlabel ('Overall efficiency')
legend_object = legend (mass_composition_array);
title (legend_object, legend_title);