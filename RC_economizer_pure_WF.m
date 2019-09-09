%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% refpropm  Thermophysical properties of pure substances and mixtures.
%   Calling sequence for pure substances:
%      result=refpropm(prop_req, spec1, value1, spec2, value2, substance1)
%
%   Calling predefined mixtures:
%      result=refpropm(prop_req, spec1, value1, spec2, value2, mixture1)
%
%   Calling user defined mixtures:
%      result=refpropm(prop_req, spec1, value1, spec2, value2,
%                                           substance1, substance2, ..., x)
%
%   where
%       prop_req    character string showing the requested properties
%                   Each property is represented by one character:
%                           0   Refprop DLL version number
%                           A   Speed of sound [m/s]
%                           B   Volumetric expansivity (beta) [1/K]
%                           C   Cp [J/(kg K)]
%                           D   Density [kg/m^3]
%                           F   Fugacity [kPa] (returned as an array)
%                           G   Gross heating value [J/kg]
%                           H   Enthalpy [J/kg]
%                           I   Surface tension [N/m]
%                           J   Isenthalpic Joule-Thompson coeff [K/kPa]
%                           K   Ratio of specific heats (Cp/Cv) [-]
%                           L   Thermal conductivity [W/(m K)]
%                           M   Molar mass [g/mol]
%                           N   Net heating value [J/kg]
%                           O   Cv [J/(kg K)]
%                           P   Pressure [kPa]
%                           Q   Quality (vapor fraction) (kg/kg)
%                           S   Entropy [J/(kg K)]
%                           T   Temperature [K]
%                           U   Internal energy [J/kg]
%                           V   Dynamic viscosity [Pa*s]
%                           X   Liquid phase & gas phase comp.(mass frac.)
%                           Y   Heat of Vaporization [J/kg]
%                           Z   Compressibility factor
%                           $   Kinematic viscosity [cm^2/s]
%                           %   Thermal diffusivity [cm^2/s]
%                           ^   Prandtl number [-]
%                           )   Adiabatic bulk modulus [kPa]
%                           |   Isothermal bulk modulus [kPa]
%                           =   Isothermal compressibility [1/kPa]
%                           ~   Cstar [-]
%                           `   Throat mass flux [kg/(m^2 s)]
%                           +   Liquid density of equilibrium phase
%                           -   Vapor density of equilibrium phase
%
%                           E   dP/dT (along the saturation line) [kPa/K]
%                           #   dP/dT     (constant rho) [kPa/K]
%                           R   d(rho)/dP (constant T)   [kg/m^3/kPa]
%                           W   d(rho)/dT (constant p)   [kg/(m^3 K)]
%                           !   dH/d(rho) (constant T)   [(J/kg)/(kg/m^3)]
%                           &   dH/d(rho) (constant P)   [(J/kg)/(kg/m^3)]
%                           (   dH/dT     (constant P)   [J/(kg K)]
%                           @   dH/dT     (constant rho) [J/(kg K)]
%                           *   dH/dP     (constant T)   [J/(kg kPa)]
%
%       spec1           first input character:  T, P, H, D, C, R, or M
%                         T, P, H, D:  see above
%                         C:  properties at the critical point
%                         R:  properties at the triple point
%                         M:  properties at Tmax and Pmax
%                            (Note: if a fluid's lower limit is higher
%                             than the triple point, the lower limit will
%                             be returned)
%
%       value1          first input value
%
%       spec2           second input character:  P, D, H, S, U or Q
%
%       value2          second input value
%
%       substance1      file name of the pure fluid (or the first
%                       component of the mixture)
%
%       mixture1        file name of the predefined fluid mixture
%                       with the extension ".mix" included
%
%       substance2,substance3,...substanceN
%                       name of the other substances in the
%                       mixture. Up to 20 substances can be handled.
%                       Valid substance names are equal to the file names
%                       in the C:\Program Files\REFPROP\fluids\' directory.
%
%       x               vector with mass fractions of the substances
%                       in the mixture.
%
%   Examples:
%   1) P = refpropm('P','T',373.15,'Q',0,'water') gives
%      Vapor pressure of water at 373.15 K in [kPa]
%
%   2) [S Cp] = refpropm('SC','T',373.15,'Q',1,'water') gives
%      Entropy and Cp of saturated steam at 373.15 K
%
%   3) D = refpropm('D','T',323.15,'P',1e2,'water','ammonia',[0.9 0.1])
%      Density of a 10% ammonia/water solution at 100 kPa and 323.15 K.
%
%   4) [x y] = refpropm('X','P',5e2,'Q',0.4,'R134a','R32',[0.8, 0.2])
%      Temperature as well as gas and liquid compositions for a mixture
%      of two refrigerants at a certain pressure and quality.
%      Note that, when 'X' is requested, two variables must be sent, the
%      first contains the liquid phase composition and the second
%      the vapor phase composition.
%
%   5) T=refpropm('T','C',0,' ',0,'water')
%      Critical temperature
%
%   6) T=refpropm('T','M',0,' ',0,'r410a.mix')
%      Maximum temperature that can be used to call properties.
%      Shows how to call a predefined mixture.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Rankine Cycle with economizer, using a pure substance as working fluid.
% Inputs: working fluid mixture, expander inlet temperature, pump outlet
% pressure, expander isentropic efficiency, pump isentropic efficiency, 
% ambient temperature, economizer effectiveness.
% Outputs: net work output, overall efficiency, working fluid's mass flow rate.

% State 1: Liquid receiver exit or pump inlet.
% State 2: Pump exit or economizer (high pressure stream) inlet.
% State 3: Economizer (high pressure stream) exit or heater inlet.
% State 4: Heater exit or expander inlet.
% State 5: Expander exit or economizer (low pressure stream) inlet.
% State 6: Economizer (low pressure stream) exit or condenser inlet.
% State 7: Condenser exit or liquid receiver inlet.

clc;
clear;

WF = 'water'  % Working fluid of Rankine cycle.
[critical_T_K, critical_P_kPa] = refpropm ('TP', 'C', 0, ' ', 0, WF);
high_T_limit_K = critical_T_K - 10;
expander_power_output_W = 100000
expander_isentropic_efficiency = 0.65
pump_isentropic_efficiency = 0.85
ambient_T_K = 30 + 273.15;
T_bubble_point_low_P_K = ambient_T_K + 15;
economizer_effectiveness = 0.8

T1_K = T_bubble_point_low_P_K;
T4_K = 500 + 273.15 % Chosen high T, in K (at heater exit or expander inlet).
P1_kPa = refpropm('P', 'T', T1_K, 'Q', 0, WF);
P2_kPa = 9200 % Chosen high P, in kPa (at pump exit or economizer inlet).
P3_kPa = 0.99 * P2_kPa; % 1% P drop in economizer (high pressure stream).
P4_kPa = 0.98 * P2_kPa; % 1% P drop in heater.
P5_kPa = 1.02 * P1_kPa; % 1% P drop in economizer (low pressure stream).
P6_kPa = 1.01 * P1_kPa; % 1% P drop in condenser.
T_dew_point_low_P_K = refpropm ('T', 'P', P6_kPa, 'Q', 1, WF);
% T_glide_low_P_K = dew_point_T_low_P_K - T1_K;

if (P2_kPa >= P1_kPa + 100) && (high_T_limit_K >= T_bubble_point_low_P_K + 10)      % Feasible conditions.
    
    if P2_kPa < critical_P_kPa          % Sub-critical ORC.
        RC_type = 'sub-critical'
        T_bubble_point_high_P_K = refpropm('T', 'P', P3_kPa, 'Q', 0, WF); % valid only for sub-critical ORC.
        T_dew_point_high_P_K = refpropm ('T', 'P', P4_kPa, 'Q', 1, WF); % valid only for sub-critical ORC.
        T_superheat_K = T4_K - T_bubble_point_high_P_K;
        if T_dew_point_high_P_K <= high_T_limit_K
            n = T_dew_point_high_P_K - T_dew_point_low_P_K;
            n = round (n);
            T_vector_K = linspace (T_dew_point_low_P_K, T_dew_point_high_P_K, n);
        else
            n = high_T_limit_K - T_dew_point_low_P_K;
            n = round (n);
            T_vector_K = linspace (T_dew_point_low_P_K, high_T_limit_K, n);
        end
        
    else                                % Trans-critical ORC.
        RC_type = 'trans-critical'
        n = high_T_limit_K - T_dew_point_low_P_K;
        n = round (n);
        T_vector_K = linspace (T_dew_point_low_P_K, high_T_limit_K, n);
    end
    
    % Finding maximum s for 0.8 Q saturated phase, in region of interest.
    saturated_s_vector_JperKkg = zeros (1, n);
    for m = 1 : n
        saturated_s_vector_JperKkg (m) = refpropm('S', 'T', T_vector_K (m), 'Q', 0.8, WF);
    end
    max_saturated_s_JperKkg = max (saturated_s_vector_JperKkg);
        
    % Finding T corresponding to high P and 0.8 Q maximum s in region of interest.
    high_P_max_s_T_K = refpropm ('T', 'P', P4_kPa, 'S', max_saturated_s_JperKkg, WF);
        
    if T4_K >= high_P_max_s_T_K

        % Expander calculations (between states 4 and 5).
        [h4_Jperkg, s4_JperKkg] = refpropm ('HS', 'T', T4_K, 'P', P4_kPa, WF);
        s5prime_JperKkg = s4_JperKkg;
        [h5prime_Jperkg, T5prime_K] = refpropm ('HT', 'P', P5_kPa, 'S', s5prime_JperKkg, WF);
        delta_h_expander_isentropic_Jperkg = h5prime_Jperkg - h4_Jperkg;
        delta_h_expander_Jperkg = delta_h_expander_isentropic_Jperkg * expander_isentropic_efficiency;
        h5_Jperkg = h4_Jperkg + delta_h_expander_Jperkg;
        [T5_K, s5_JperKkg] = refpropm('TS', 'P', P5_kPa, 'H', h5_Jperkg, WF);
        expander_specific_work_output_Jperkg = h4_Jperkg - h5_Jperkg;

        % Pump calculations (between states 1 and 2).
        [h1_Jperkg, s1_JperKkg] = refpropm('HS', 'T', T1_K, 'Q', 0, WF);
        s2prime_JperKkg = s1_JperKkg;
        [h2prime_Jperkg, T2prime_K] = refpropm('HT', 'P', P2_kPa, 'S', s2prime_JperKkg, WF);
        delta_h_pump_isentropic_Jperkg = h2prime_Jperkg - h1_Jperkg;
        delta_h_pump_Jperkg = delta_h_pump_isentropic_Jperkg / pump_isentropic_efficiency;
        h2_Jperkg = h1_Jperkg + delta_h_pump_Jperkg;
        [T2_K, s2_JperKkg] = refpropm('TS', 'P', P2_kPa, 'H', h2_Jperkg, WF);
        pump_specific_work_input_Jperkg = h2_Jperkg - h1_Jperkg;
        
        % Economizer calculations.
        if T5_K >= T2_K + 20 % Sufficient T difference. Economizer feasible.
            economizer = 'yes'
            Cp_5_JperKkg = refpropm ('C', 'P', P5_kPa, 'H', h5_Jperkg, WF);
            Cp_2_JperKkg = refpropm ('C', 'P', P2_kPa, 'H', h2_Jperkg, WF);
            Cp_min_JperKkg = min (Cp_5_JperKkg, Cp_2_JperKkg);
            q_max_Jperkg = Cp_min_JperKkg * (T5_K - T2_K);
            q_low_P_Jperkg = economizer_effectiveness * q_max_Jperkg;
            q_high_P_Jperkg = 0.95 * q_low_P_Jperkg;                       % Assuming that, in economizer, 95 % of heat lost by low pressure stream is gained by high pressure stream.
            h3_Jperkg = h2_Jperkg + q_high_P_Jperkg;
            h6_Jperkg = h5_Jperkg - q_low_P_Jperkg;
            T3_K = refpropm('T', 'H', h3_Jperkg, 'P', P3_kPa, WF);
            T6_K = refpropm('T', 'H', h6_Jperkg, 'P', P6_kPa, WF);
            s3_JperKkg = refpropm('S', 'H', h3_Jperkg, 'P', P3_kPa, WF);
            s6_JperKkg = refpropm('S', 'H', h6_Jperkg, 'P', P6_kPa, WF);
        else          % Insufficient T difference. Economizer not feasible.
            economizer = 'no'
            T6_K = T5_K;
            T3_K = T2_K;
            h6_Jperkg = h5_Jperkg;
            h3_Jperkg = h2_Jperkg;
            s3_JperKkg = s2_JperKkg;
            s6_JperKkg = s5_JperKkg;
        end
        
        % Calculating net work output, heat input, overall efficiency, and working fluid's mass flow rate.
        net_specific_work_output_Jperkg = expander_specific_work_output_Jperkg - pump_specific_work_input_Jperkg
        % net_power_output_W = mass_flow_rate_kgpers * net_specific_work_output_Jperkg
        mass_flow_rate_kgpers = expander_power_output_W / expander_specific_work_output_Jperkg
        heat_input_Jperkg = h4_Jperkg - h3_Jperkg;
        heat_input_W = mass_flow_rate_kgpers * heat_input_Jperkg
        overall_efficiency = net_specific_work_output_Jperkg / heat_input_Jperkg
            
    else
        error (['To avoid too much condensation in expander, selected expander inlet temperature (T4) must be at least ', num2str(high_P_max_s_T_K), 'K, which is fluid temperature corresponding to expander inlet pressure (P4 = ', num2str(P4_kPa), ' kPa) and maximum specific entropy for 0.8 quality (', num2str(max_saturated_s_JperKkg), ' J/Kkg) in region of interest.'])
    end
    
else           
    if P2_kPa < P1_kPa + 100     % Too little P difference for useful RC.
        error (['Selected pump outlet pressure (P2) should be at least 100 kPa more than ', num2str(P1_kPa), ' kPa which is fluid condensation pressure for ', num2str(T1_K), ' K.'])
    end
    if high_T_limit_K < T_dew_point_low_P_K + 10   % Too little T difference for analysis.
        error(['There should be at least 20 K temperature difference between critical temperature (currently ', num2str(critical_T_K), ' K) and dew point temperature at low pressure (currently ', num2str(T_dew_point_low_P_K), ' K).'])
    end
end