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

% For a Rankine cycle, given the pure working fluid, condensation
% temperature, expander's built-in volume ratio, and expander inlet
% temperature, this program calculates the expander inlet pressure, 
% expander outlet temperature, pressure ratios, etc., through an iterative 
% process. Expansion process is assumed to be adiabatic.

clc;
clear;

WF = 'R245fa';
T_ambient_K = 30 + 273.15;
T1_K = T_ambient_K + 15;
P1_kPa = refpropm ('P', 'T', T1_K, 'Q', 0, WF);
P4_kPa = 1.01 * P1_kPa;     % Assuming 1 % pressure drop in condenser.
T3_K = 140 + 273.15;        % Expander inlet temperature.
V_4_3_ratio = 2;            % Built-in volume ratio of expander.

first_T4_K = refpropm ('T', 'P', P4_kPa, 'Q', 1, WF);
first_P3_kPa = 2 * P4_kPa;

gamma3 = refpropm ('K', 'T', T3_K, 'P', first_P3_kPa, WF);
gamma4 = refpropm ('K', 'T', first_T4_K, 'P', P4_kPa, WF);
gamma_average = (gamma3 + gamma4) / 2;

second_P3_kPa = P4_kPa * (V_4_3_ratio ^ gamma_average);
second_T4_K = T3_K / (V_4_3_ratio ^ (gamma_average - 1));

iteration = 1;      % First iteration already done above.

while ((abs(first_P3_kPa - second_P3_kPa) > 0.01) || (abs(first_T4_K - second_T4_K) > 0.01))
    
    first_T4_K = second_T4_K;
    first_P3_kPa = second_P3_kPa;
    
    gamma3 = refpropm ('K', 'T', T3_K, 'P', first_P3_kPa, WF);
    gamma4 = refpropm ('K', 'T', first_T4_K, 'P', P4_kPa, WF);
    gamma_average = (gamma3 + gamma4) / 2;      % gamma is adiabatic index (Cp / Cv).

    second_P3_kPa = P4_kPa * (V_4_3_ratio ^ gamma_average);
    second_T4_K = T3_K / (V_4_3_ratio ^ (gamma_average - 1));

    iteration = iteration + 1;
end

iteration
second_P3_kPa
second_T4_K
P2_kPa = second_P3_kPa / 0.99   % Assuming 1 % pressure drop in evaporator.
P_3_4_ratio = second_P3_kPa / P4_kPa
P_2_1_ratio = P2_kPa / P1_kPa
