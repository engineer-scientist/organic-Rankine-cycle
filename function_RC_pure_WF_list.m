function [selected_pure_fluids_list] = function_RC_pure_WF_list ()
ambient_T_high_K = 40 + 273.15;
ambient_T_low_K = 10 + 273.15;
condensation_T_high_K = 10 + ambient_T_high_K;
condensation_T_low_K = 10 + ambient_T_low_K;

% ......................... PURE FLUIDS ..............................

refprop_pure_fluids_list = {'R32', 'H2S', '1butene', 'acetone', 'ammonia', 'argon', 'benzene', 'butane', 'C1CC6', 'C2butene', 'C3CC6', 'C4F10', 'C5F12', 'C12', 'CF3I', 'CO', 'CO2', 'COS', 'cyclohex', 'cyclopen', 'cyclopro', 'D2', 'D2O', 'D4', 'D5', 'D6', 'decane', 'DMC', 'DME', 'ethane', 'ethanol', 'ethylene', 'fluorine', 'helium', 'heptane', 'hexane', 'hydrogen', 'ibutene', 'ihexane', 'ipentane', 'isobutan', 'krypton', 'MD2M', 'MD3M', 'MD4M', 'MDM',  'methane', 'methanol', 'mlinolea', 'mlinolen', 'MM', 'moleate', 'mpalmita', 'mstearat', 'N2O', 'neon', 'neopentn', 'NF3', 'nitrogen', 'nonane', 'octane', 'orthohyd', 'oxygen', 'parahyd', 'pentane', 'propane', 'propylen', 'propyne', 'R11', 'R12', 'R13', 'R14', 'R21', 'R22', 'R23', 'R41', 'R113', 'R114', 'R115', 'R116', 'R123', 'R124', 'R125', 'R134a', 'R141b', 'R142b', 'R143a', 'R152a', 'R161', 'R218', 'R227ea', 'R236ea', 'R236fa', 'R245ca', 'R245fa', 'R365mfc', 'R1234yf', 'R1234ze', 'RC318', 'SF6', 'SO2', 't2butene', 'toluene', 'water', 'xenon'}; % 'R32', 'H2S', };
n_R_pure_F = numel(refprop_pure_fluids_list);
critical_T_pure_K = zeros(1, n_R_pure_F);
triple_T_pure_K = zeros(1, n_R_pure_F);

for j = 1 : n_R_pure_F
    critical_T_pure_K(j) = refpropm('T', 'C', 0, ' ', 0, refprop_pure_fluids_list{j});
    triple_T_pure_K(j) = refpropm('T', 'R', 0, ' ', 0, refprop_pure_fluids_list{j});
end

L = 1;
for k = 1 : n_R_pure_F
    if ((triple_T_pure_K (k) + 10 < condensation_T_low_K) && (condensation_T_high_K < critical_T_pure_K (k) - 10))
        selected_pure_fluids_list{L} = refprop_pure_fluids_list{k};
        % S_pure_F_critical_T_K(L) = critical_T_pure_K(k);
        % S_pure_F_triple_T_K(L) = triple_T_pure_K(k);
        L = L + 1;
    end
end

% ......... MIXTURES (can not calculate triple point)......................
% [Error using refpropm (line 368). Triple point not known for mixtures.]

end

