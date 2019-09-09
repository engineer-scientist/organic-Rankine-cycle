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
