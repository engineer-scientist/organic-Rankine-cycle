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

% For plotting temperature vs. specific entropy and pressure vs. specific
% enthalpy diagrams of selected fluids.

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
