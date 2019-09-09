function [net_work_output_Jperkg, overall_efficiency] = function_RC_economizer_input_high_T_high_P_3_C_mixture_WF (expander_isentropic_efficiency, pump_isentropic_efficiency, ambient_T_K, economizer_effectiveness, T4_K, P2_kPa, mass_composition, substance1, substance2, substance3)

% Rankine Cycle with economizer, using a mixture as working fluid
% (can have a pure substance like element or compound as working fluid, as 
% a special case).
% Input parameters: working fluid mixture, expander inlet temperature,
% expander inlet pressure, expander isentropic efficiency, pump isentropic
% efficiency, ambient temperature.
% Output parameters: net work output, overall efficiency.

% State 1: Liquid receiver exit or pump inlet.
% State 2: Pump exit or economizer (high pressure stream) inlet.
% State 3: Economizer (high pressure stream) exit or heater inlet.
% State 4: Heater exit or expander inlet.
% State 5: Expander exit or economizer (low pressure stream) inlet.
% State 6: Economizer (low pressure stream) exit or condenser inlet.
% State 7: Condenser exit or liquid receiver inlet.

try
    
[critical_T_K, critical_P_kPa] = refpropm ('TP', 'C', 0, ' ', 0, substance1, substance2, substance3, mass_composition);
high_T_limit_K = critical_T_K - 10;

ambient_T_condensation_P_kPa = refpropm ('P', 'T', ambient_T_K, 'Q', 0, substance1, substance2, substance3, mass_composition);
dew_point_T_initial_condensation_P_K = refpropm ('T', 'P', ambient_T_condensation_P_kPa, 'Q', 1, substance1, substance2, substance3, mass_composition);
initial_T_glide_K = dew_point_T_initial_condensation_P_K - ambient_T_K;

% Fixing T1. (Bubble point T for low P, and lowest T of WF. At pump inlet.)
if initial_T_glide_K > 20
    T1_K = ambient_T_K + 5;
elseif (initial_T_glide_K <= 20) && (initial_T_glide_K > 10)
    T1_K = ambient_T_K + 10;
else
    T1_K = ambient_T_K + 15;
end

P1_kPa = refpropm('P', 'T', T1_K, 'Q', 0, substance1, substance2, substance3, mass_composition);
% P3_kPa = 0.99 * P2_kPa; % 1% P drop in economizer (high pressure stream).
P4_kPa = 0.98 * P2_kPa; % 1% P drop in heater.
P5_kPa = 1.02 * P1_kPa; % 1% P drop in economizer (low pressure stream).
P6_kPa = 1.01 * P1_kPa; % 1% P drop in condenser.
dew_point_T_low_P_K = refpropm ('T', 'P', P6_kPa, 'Q', 1, substance1, substance2, substance3, mass_composition);
% T_glide_low_P_K = dew_point_T_low_P_K - T1_K;

if (P2_kPa >= P1_kPa + 100) && (high_T_limit_K >= dew_point_T_low_P_K + 10)      % Feasible conditions.
    
    if P2_kPa < critical_P_kPa          % Sub-critical ORC.
        RC_type = 'sub-critical';
        % bubble_point_T_high_P_K = refpropm('T', 'P', P2_kPa, 'Q', 0, substance1, substance2, mass_composition); % valid only for sub-critical ORC.
        dew_point_T_high_P_K = refpropm ('T', 'P', P4_kPa, 'Q', 1, substance1, substance2, substance3, mass_composition); % valid only for sub-critical ORC.
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
        RC_type = 'trans-critical';
        n = high_T_limit_K - dew_point_T_low_P_K;
        n = round (n);
        T_vector_K = linspace (dew_point_T_low_P_K, high_T_limit_K, n);
    end
    
    % Finding maximum s for 0.8 Q saturated phase, in region of interest.
    saturated_s_vector_JperKkg = zeros (1, n);
    for m = 1 : n
        saturated_s_vector_JperKkg (m) = refpropm('S', 'T', T_vector_K (m), 'Q', 0.8, substance1, substance2, substance3, mass_composition);
    end
    max_saturated_s_JperKkg = max (saturated_s_vector_JperKkg);
        
    % Finding T corresponding to high P and 0.8 Q maximum s in region of interest.
    high_P_max_s_T_K = refpropm ('T', 'P', P4_kPa, 'S', max_saturated_s_JperKkg, substance1, substance2, substance3, mass_composition);
        
    if T4_K >= high_P_max_s_T_K

        % Calculating expander work output (between states 4 and 5).
        [h4_Jperkg, s4_JperKkg] = refpropm ('HS', 'T', T4_K, 'P', P4_kPa, substance1, substance2, substance3, mass_composition);
        s5prime_JperKkg = s4_JperKkg;
        h5prime_Jperkg = refpropm ('H', 'P', P5_kPa, 'S', s5prime_JperKkg, substance1, substance2, substance3, mass_composition);
        delta_h_expander_isentropic_Jperkg = h5prime_Jperkg - h4_Jperkg;
        delta_h_expander_Jperkg = delta_h_expander_isentropic_Jperkg * expander_isentropic_efficiency;
        h5_Jperkg = h4_Jperkg + delta_h_expander_Jperkg;
        T5_K = refpropm('T', 'P', P5_kPa, 'H', h5_Jperkg, substance1, substance2, substance3, mass_composition);
        expander_work_output_Jperkg = h4_Jperkg - h5_Jperkg;

        % Calculating pump work input (between states 1 and 2).
        [h1_Jperkg, s1_JperKkg] = refpropm('HS', 'T', T1_K, 'Q', 0, substance1, substance2, substance3, mass_composition);
        s2prime_JperKkg = s1_JperKkg;
        h2prime_Jperkg = refpropm('H', 'P', P2_kPa, 'S', s2prime_JperKkg, substance1, substance2, substance3, mass_composition);
        delta_h_pump_isentropic_Jperkg = h2prime_Jperkg - h1_Jperkg;
        delta_h_pump_Jperkg = delta_h_pump_isentropic_Jperkg / pump_isentropic_efficiency;
        h2_Jperkg = h1_Jperkg + delta_h_pump_Jperkg;
        T2_K = refpropm('T', 'P', P2_kPa, 'H', h2_Jperkg, substance1, substance2, substance3, mass_composition);
        pump_work_input_Jperkg = h2_Jperkg - h1_Jperkg;
        
        % Economizer calculations.
        if T5_K >= T2_K + 20 % Sufficient T difference. Economizer feasible.
            economizer = 'yes';
            Cp_5_JperKkg = refpropm ('C', 'P', P5_kPa, 'H', h5_Jperkg, substance1, substance2, substance3, mass_composition);
            Cp_2_JperKkg = refpropm ('C', 'P', P2_kPa, 'H', h2_Jperkg, substance1, substance2, substance3, mass_composition);
            Cp_min_JperKkg = min (Cp_5_JperKkg, Cp_2_JperKkg);
            q_max_Jperkg = Cp_min_JperKkg * (T5_K - T2_K);
            q_low_P_Jperkg = economizer_effectiveness * q_max_Jperkg;
            q_high_P_Jperkg = 0.95 * q_low_P_Jperkg;                       % Assuming that, in economizer, 95 % of heat lost by low pressure stream is gained by high pressure stream.
            h3_Jperkg = h2_Jperkg + q_high_P_Jperkg;
            % T6_K = T2_K + 10;
            % h6_Jperkg = refpropm('H', 'T', T6_K, 'P', P6_kPa, substance1, substance2, mass_composition);
            % h3_Jperkg = h2_Jperkg + 0.9 * (h5_Jperkg - h6_Jperkg);         
            % T3_K = refpropm('T', 'H', h3_Jperkg, 'P', P3_kPa, substance1, substance2, mass_composition);
        else          % Insufficient T difference. Economizer not feasible.
            economizer = 'no';
            % T6_K = T5_K;
            % h6_Jperkg = h5_Jperkg;
            h3_Jperkg = h2_Jperkg;
            % T3_K = T2_K;
        end

        % Calculating net work output, heat input, and overall efficiciency.
        net_work_output_Jperkg = expander_work_output_Jperkg - pump_work_input_Jperkg;
        heat_input_Jperkg = h4_Jperkg - h3_Jperkg;
        overall_efficiency = net_work_output_Jperkg / heat_input_Jperkg;
            
    else
        error (['To avoid too much condensation in expander, selected expander inlet temperature (T3) must be at least ', num2str(high_P_max_s_T_K), 'K, which is fluid temperature corresponding to expander inlet pressure (P4 = ', num2str(P4_kPa), ' kPa) and maximum specific entropy for 0.8 quality (', num2str(max_saturated_s_JperKkg), ' J/Kkg) in region of interest.'])
    end
    
else           
    if P2_kPa < P1_kPa + 100     % Too little P difference for useful ORC.
        error (['Selected high pressure (P2) should be at least 100 kPa more than ', num2str(P1_kPa), ' kPa which is fluid condensation pressure for ', num2str(T1_K), ' K.'])
    end
    if high_T_limit_K < dew_point_T_low_P_K + 10                           % Too little T difference for analysis.
        error(['There should be at least 20 K temperature difference between critical temperature (currently ', num2str(critical_T_K), ' K) and dew point temperature at low pressure (currently ', num2str(dew_point_T_low_P_K), ' K).'])
    end
end

catch
    net_work_output_Jperkg = NaN;
    overall_efficiency = NaN;
    RC_type = NaN;
    economizer = NaN;
end

end
