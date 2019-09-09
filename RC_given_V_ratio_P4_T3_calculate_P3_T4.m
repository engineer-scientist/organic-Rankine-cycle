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
