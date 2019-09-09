% Rankine Cycle without economizer, using a pure substance as working fluid.
% Input parameters: working fluid mixture, expander inlet temperature, expander inlet pressure.
% Output parameters: net work output, overall efficiency.

% State 1: Liquid receiver exit or pump inlet.
% State 2: Pump exit or heater inlet.
% State 3: Heater exit or expander inlet.
% State 4: Expander exit or condenser inlet.
% State 5: Condenser exit or liquid receiver inlet.

clc;
clear;

WF = 'R245fa'   % Working fluid for Rankine cycle.
[critical_T_K, critical_P_kPa] = refpropm('TP', 'C', 0, ' ', 0, WF);
T_high_limit_K = critical_T_K - 10;
mass_flow_rate_kgpers = 0.35
% expander_power_output_W = 100000
% auxilliary_load_W = 30000;
expander_isentropic_efficiency = 0.65;
pump_isentropic_efficiency = 0.85;
ambient_T_K = 30 + 273.15;
T_condensation_K = ambient_T_K + 15;

T1_K = T_condensation_K;

P1_kPa = refpropm('P', 'T', T1_K, 'Q', 0, WF);
P2_kPa = 637.01   % Chosen high pressure (at pump exit or heater inlet).
P3_kPa = 0.99 * P2_kPa;  % 1% pressure drop in boiler.
P4_kPa = 1.01 * P1_kPa; % 1% pressure drop in condenser.
T3_K = 140 + 273.15    % Heater exit or expander inlet temperature.
T_dew_point_low_P_K = refpropm('T', 'P', P4_kPa, 'Q', 1, WF);

if (P2_kPa >= P1_kPa + 100) && (T_high_limit_K >= T_dew_point_low_P_K + 10)      % Feasible conditions.
    
    if P2_kPa < critical_P_kPa        % Sub-critical ORC.

        RC_type = 'sub-critical'
        T_dew_point_high_P_K = refpropm('T', 'P', P3_kPa, 'Q', 1, WF); % valid only for sub-critical ORC.
        T_bubble_point_high_P_K = refpropm('T', 'P', P2_kPa, 'Q', 0, WF); % valid only for sub-critical ORC.
        T_evaporation_K = (T_dew_point_high_P_K + T_bubble_point_high_P_K) / 2;
        % Finding maximum s for 0.8 Q saturated phase, between dew point T (low P), and dew point T (high P) or high T limit (whichever is less).
        n = 100;
        if T_dew_point_high_P_K <= T_high_limit_K
            T_vector_K = linspace (T_dew_point_low_P_K, T_dew_point_high_P_K, n);
        else
            T_vector_K = linspace (T_dew_point_low_P_K, T_high_limit_K, n);
        end
        saturated_s_vector_JperKkg = zeros (1, n);
        for m = 1 : n
            saturated_s_vector_JperKkg (m) = refpropm('S', 'T', T_vector_K (m), 'Q', 0.8, WF);
        end
        max_saturated_s_JperKkg = max (saturated_s_vector_JperKkg);
        
        % Finding T corresponding to high P and saturated maximum s between dew point T (low P) and dew point T (high P).
        high_P_max_s_T_K = refpropm('T', 'P', P3_kPa, 'S', max_saturated_s_JperKkg, WF);
        
    else        % Trans-critical ORC.
        
        RC_type = 'trans-critical'
        
        T_dew_point_low_P_K = refpropm('T', 'P', P4_kPa, 'Q', 1, WF);
        
        % Finding maximum s for 0.8 Q saturated phase between dew point T (low P) and high T limit.
        n = 100;
        T_vector_K = linspace (T_dew_point_low_P_K, T_high_limit_K, n);
        saturated_s_vector_JperKkg = zeros (1, n);
        for m = 1 : n
            saturated_s_vector_JperKkg (m) = refpropm('S', 'T', T_vector_K (m), 'Q', 0.8, WF);
        end
        max_saturated_s_JperKkg = max (saturated_s_vector_JperKkg);
        
        % Finding T corresponding to high P and saturated maximum s between dew point T (low P) and critical T.
        high_P_max_s_T_K = refpropm('T', 'P', P3_kPa, 'S', max_saturated_s_JperKkg, WF);
    
    end
        
    if T3_K >= high_P_max_s_T_K
        
        % Expander calculations.
        [h3_Jperkg, s3_JperKkg] = refpropm('HS', 'T', T3_K, 'P', P3_kPa, WF);
        s4prime_JperKkg = s3_JperKkg;
        [h4prime_Jperkg, T4prime_K] = refpropm('HT', 'P', P4_kPa, 'S', s4prime_JperKkg, WF);
        delta_h_expander_isentropic_Jperkg = h4prime_Jperkg - h3_Jperkg;
        delta_h_expander_Jperkg = delta_h_expander_isentropic_Jperkg * expander_isentropic_efficiency;
        h4_Jperkg = h3_Jperkg + delta_h_expander_Jperkg;
        [T4_K, s4_JperKkg] = refpropm('TS', 'P', P4_kPa, 'H', h4_Jperkg, WF);
        expander_specific_work_output_Jperkg = h3_Jperkg - h4_Jperkg;

        % Pump calculations.
        [h1_Jperkg, s1_JperKkg] = refpropm('HS', 'T', T1_K, 'Q', 0, WF);
        s2prime_JperKkg = s1_JperKkg;
        [h2prime_Jperkg, T2prime_K] = refpropm('HT', 'P', P2_kPa, 'S', s2prime_JperKkg, WF);
        delta_h_pump_isentropic_Jperkg = h2prime_Jperkg - h1_Jperkg;
        delta_h_pump_Jperkg = delta_h_pump_isentropic_Jperkg / pump_isentropic_efficiency;
        h2_Jperkg = h1_Jperkg + delta_h_pump_Jperkg;
        [T2_K, s2_JperKkg] = refpropm('TS', 'P', P2_kPa, 'H', h2_Jperkg, WF);
        pump_specific_work_input_Jperkg = h2_Jperkg - h1_Jperkg;

        % Calculating net work output, heat input, and overall efficiciency.
        net_specific_work_output_Jperkg = expander_specific_work_output_Jperkg - pump_specific_work_input_Jperkg
        net_power_output_W = mass_flow_rate_kgpers * net_specific_work_output_Jperkg
        % mass_flow_rate_kgpers = expander_power_output_W / expander_specific_work_output_Jperkg
        heat_input_specific_Jperkg = h3_Jperkg - h2_Jperkg;
        heat_input_power_W = mass_flow_rate_kgpers * heat_input_specific_Jperkg
        pump_input_power_W = mass_flow_rate_kgpers * pump_specific_work_input_Jperkg
        % net_power_output_W = net_specific_work_output_Jperkg * WF_mass_flow_rate_kgpers
        % other_parasitic_power_W = auxilliary_load_W - pump_input_power_W
        overall_efficiency = net_specific_work_output_Jperkg / heat_input_specific_Jperkg
          
    else
        error (['To avoid too much condensation in expander, selected expander inlet temperature (T3) must be greater than ', num2str(high_P_max_s_T_K), 'K, which is fluid temperature corresponding to expander inlet pressure (P3 = ', num2str(P3_kPa), ' kPa) and maximum specific entropy for 0.8 quality (', num2str(max_saturated_s_JperKkg), ' J/Kkg).'])
    end       
    
else           
    if P2_kPa < P1_kPa + 100     % Too little P difference for useful ORC.
        error (['Selected high pressure (P2) should be at least 100 kPa more than fluid condensation pressure for ', num2str(T1_K), ' K which is ', num2str(P1_kPa), ' kPa.'])
    end
    if T_high_limit_K < T_dew_point_low_P_K + 10                           % Too little T difference for analysis.
        error(['There should be at least 20 K temperature difference between critical temperature (currently ', num2str(critical_T_K), ' K) and dew point temperature at low pressure (currently ', num2str(T_dew_point_low_P_K), ' K).'])
    end
end
