function [net_work_output_Jperkg, overall_efficiency] = function_RC_input_high_T_high_P_pure_WF (expander_isentropic_efficiency, pump_isentropic_efficiency, ambient_T_K, T3_K, P2_kPa, working_fluid)

% Organic Rankine Cycle without economizer, using a pure substance as
% working fluid.
% Input parameters: working fluid mixture, expander inlet temperature,
% expander inlet pressure, expander isentropic efficiency, pump isentropic
% efficiency, ambient temperature.
% Output parameters: net work output, overall efficiency.

% State 1: Liquid receiver exit or pump inlet.
% State 2: Pump exit or heater inlet.
% State 3: Heater exit or expander inlet.
% State 4: Expander exit or condenser inlet.
% State 5: Condenser exit or liquid receiver inlet.

% try
    
[critical_T_K, critical_P_kPa] = refpropm ('TP', 'C', 0, ' ', 0, working_fluid);
high_T_limit_K = critical_T_K - 10;

% Fixing T1. (Condensation temperature.)
T1_K = ambient_T_K + 15;

P1_kPa = refpropm('P', 'T', T1_K, 'Q', 0, working_fluid);
P3_kPa = 0.99 * P2_kPa; % 1% P drop in heater.
P4_kPa = 1.01 * P1_kPa; % 1% P drop in condenser.
dew_point_T_low_P_K = refpropm ('T', 'P', P4_kPa, 'Q', 1, working_fluid);

if (P2_kPa >= P1_kPa + 100) && (high_T_limit_K >= dew_point_T_low_P_K + 10)      % Feasible conditions.
    
           
    if P2_kPa < critical_P_kPa          % Sub-critical ORC.
        % ORC_type = 'sub-critical';
        % bubble_point_T_high_P_K = refpropm('T', 'P', P2_kPa, 'Q', 0, working_fluid); % valid only for sub-critical ORC.
        dew_point_T_high_P_K = refpropm ('T', 'P', P3_kPa, 'Q', 1, working_fluid); % valid only for sub-critical ORC.
        if dew_point_T_high_P_K <= high_T_limit_K
            n = dew_point_T_high_P_K - dew_point_T_low_P_K;
            n = round (n);
            T_vector_K = linspace (dew_point_T_low_P_K, dew_point_T_high_P_K, n);
        else
            n = high_T_limit_K - dew_point_T_low_P_K;
            n = round (n);
            T_vector_K = linspace (dew_point_T_low_P_K, high_T_limit_K, n);
        end
        
    else                                % Trans-critical ORC.
        % ORC_type = 'trans-critical';
        n = high_T_limit_K - dew_point_T_low_P_K;
        n = round (n);
        T_vector_K = linspace (dew_point_T_low_P_K, high_T_limit_K, n);
    end
    
    % Finding maximum s for 0.8 Q saturated phase, in region of interest.
    saturated_s_vector_JperKkg = zeros (1, n);
    for m = 1 : n
        saturated_s_vector_JperKkg (m) = refpropm('S', 'T', T_vector_K (m), 'Q', 0.8, working_fluid);
    end
    max_saturated_s_JperKkg = max (saturated_s_vector_JperKkg);
        
    % Finding T corresponding to high P and 0.8 Q maximum s in region of interest.
    high_P_max_s_T_K = refpropm ('T', 'P', P3_kPa, 'S', max_saturated_s_JperKkg, working_fluid);
        
    if T3_K >= high_P_max_s_T_K

        % Calculating expander work output (between states 3 and 4).
        [h3_Jperkg, s3_JperKkg] = refpropm ('HS', 'T', T3_K, 'P', P3_kPa, working_fluid);
        s4prime_JperKkg = s3_JperKkg;
        h4prime_Jperkg = refpropm ('H', 'P', P4_kPa, 'S', s4prime_JperKkg, working_fluid);
        delta_h_expander_isentropic_Jperkg = h4prime_Jperkg - h3_Jperkg;
        delta_h_expander_Jperkg = delta_h_expander_isentropic_Jperkg * expander_isentropic_efficiency;
        h4_Jperkg = h3_Jperkg + delta_h_expander_Jperkg;
        expander_work_output_Jperkg = h3_Jperkg - h4_Jperkg;

        % Calculating pump work input (between states 1 and 2).
        [h1_Jperkg, s1_JperKkg] = refpropm('HS', 'T', T1_K, 'Q', 0, working_fluid);
        s2prime_JperKkg = s1_JperKkg;
        h2prime_Jperkg = refpropm('H', 'P', P2_kPa, 'S', s2prime_JperKkg, working_fluid);
        delta_h_pump_isentropic_Jperkg = h2prime_Jperkg - h1_Jperkg;
        delta_h_pump_Jperkg = delta_h_pump_isentropic_Jperkg / pump_isentropic_efficiency;
        h2_Jperkg = h1_Jperkg + delta_h_pump_Jperkg;
        pump_work_input_Jperkg = h2_Jperkg - h1_Jperkg;
              
        % Calculating net work output, heat input, and overall efficiciency.
        net_work_output_Jperkg = expander_work_output_Jperkg - pump_work_input_Jperkg;
        heat_input_Jperkg = h3_Jperkg - h2_Jperkg;
        overall_efficiency = net_work_output_Jperkg / heat_input_Jperkg;
            
    else
        error (['To avoid too much condensation in expander, selected expander inlet temperature (T3) must be at least ', num2str(high_P_max_s_T_K), 'K, which is fluid temperature corresponding to expander inlet pressure (P4 = ', num2str(P3_kPa), ' kPa) and maximum specific entropy for 0.8 quality (', num2str(max_saturated_s_JperKkg), ' J/Kkg) in region of interest.'])
    end
    
else           
    if P2_kPa < P1_kPa + 100     % Too little P difference for useful ORC.
        error (['Selected high pressure (P2) should be at least 100 kPa more than ', num2str(P1_kPa), ' kPa which is fluid condensation pressure for ', num2str(T1_K), ' K.'])
    end
    if high_T_limit_K < dew_point_T_low_P_K + 10                           % Too little T difference for analysis.
        error(['There should be at least 20 K temperature difference between critical temperature (currently ', num2str(critical_T_K), ' K) and dew point temperature at low pressure (currently ', num2str(dew_point_T_low_P_K), ' K).'])
    end
end

% catch
    % net_work_output_Jperkg = NaN;
    % overall_efficiency = NaN;
    % ORC_type = NaN;
    % economizer = NaN;
% end

end
