% For plotting temperature vs. specific entropy and pressure vs. specific enthalpy diagrams of selected fluids.

clc;
clear;

selected_pure_fluids_list = {'Water', 'R245fa', 'Toluene', 'Ethanol', 'Acetone', 'Benzene', 'Mlinolen'};
n_F = numel(selected_pure_fluids_list);

% color_matrix = zeros (3, n_F);
% for F_i = 1 : n_F
%     for c = 1 : 3
%         color_matrix (c, F_i) = rand;
%     end
% end

colors_7_list = {'blue', 'red', 'green', 'yellow', 'magenta', 'cyan', 'black'};

% ------------- Temperature vs. specific entropy diagram. -------------

figure ('Name', 'Temperature vs. specific entropy');

for F_i = 1 : n_F
    
    F = selected_pure_fluids_list {F_i};
    critical_T_K = refpropm ('T', 'C', 0, ' ', 0, F);
    T_max_integer_K = round (critical_T_K) + 1;
    T_list_K = 0 : T_max_integer_K;
    s_liquid_list_JperKkg = zeros (1, T_max_integer_K + 1);
    s_vapour_list_JperKkg = zeros (1, T_max_integer_K + 1);

    for T_K = 0 : T_max_integer_K
        try
            s_liquid_list_JperKkg (T_K + 1) = refpropm ('S', 'T', T_K, 'Q', 0, F);
        catch
            s_liquid_list_JperKkg (T_K + 1) = NaN;
        end
        try
            s_vapour_list_JperKkg (T_K + 1) = refpropm ('S', 'T', T_K, 'Q', 1, F);
        catch        
            s_vapour_list_JperKkg (T_K + 1) = NaN;
        end
    end
    
    s_total_list_JperKkg = cat (2, s_liquid_list_JperKkg, s_vapour_list_JperKkg);
    T_double_list_K = cat (2, T_list_K, T_list_K);
    plot (s_total_list_JperKkg, T_double_list_K, 'LineWidth', 1, 'Color', colors_7_list {F_i});
    hold on
    
end

xlabel ('Specific entropy (J/K-kg)')
ylabel ('Temperature (K)')
legend (selected_pure_fluids_list)

% ------------- Pressure vs. specific enthalpy diagram. -------------

figure ('Name', 'Pressure vs. specific enthalpy');

for F_i = 1 : n_F

    F = selected_pure_fluids_list {F_i};
    P_critical_kPa = refpropm ('P', 'C', 0, ' ', 0, F);
    P_max_integer_kPa = round (P_critical_kPa) + 1;
    P_list_kPa = 0 : P_max_integer_kPa;
    h_liquid_list_Jperkg = zeros (1, P_max_integer_kPa + 1);
    h_vapour_list_Jperkg = zeros (1, P_max_integer_kPa + 1);

    for P_kPa = 0 : P_max_integer_kPa
        try
            h_liquid_list_Jperkg (P_kPa + 1) = refpropm ('H', 'P', P_kPa, 'Q', 0, F);
        catch
            h_liquid_list_Jperkg (P_kPa + 1) = NaN;
        end
        try
            h_vapour_list_Jperkg (P_kPa + 1) = refpropm ('H', 'P', P_kPa, 'Q', 1, F);
        catch        
            h_vapour_list_Jperkg (P_kPa + 1) = NaN;
        end
    end
    
    h_total_list_Jperkg = cat (2, h_liquid_list_Jperkg, h_vapour_list_Jperkg);
    P_double_list_kPa = cat (2, P_list_kPa, P_list_kPa);
    plot (h_total_list_Jperkg, P_double_list_kPa, 'LineWidth', 1, 'Color', colors_7_list {F_i});
    hold on

end
xlabel ('Specific enthalpy (J/kg)')
ylabel ('Pressure (kPa)')
legend (selected_pure_fluids_list)
