% For step = 0.05 and 20 components: p = 41,685,448,298. [Run time = 24059.183 s (6 hours, 40 minutes, 59.183 seconds) on Intel Core i7-6700 CPU; 3.4 GHz, 3.4 GHz; 64 GB RAM.]
% For step = 0.1 and 20 components: p = 16,648,863. [Run time = 17.192 s on Intel Core i7-6700 CPU; 3.4 GHz, 3.4 GHz; 64 GB RAM.]
% For step = 0.05 and 10 components: p = 8,060,931. [Run time = 2.595 s on Intel Core i7-6700 CPU; 3.4 GHz, 3.4 GHz; 64 GB RAM.]
% For step = 0.1 and 10 components: p = 76,290. [Run time = 0.077 s on Intel Core i7-6700 CPU; 3.4 GHz, 3.4 GHz; 64 GB RAM.]

p = 0;
step = 0.1;
for c1 = 0 : step : 1
    for c2 = 0 : step : 1 - c1
        for c3 = 0 : step : 1 - c1 - c2
            for c4 = 0 : step : 1 - c1 - c2 - c3
                for c5 = 0 : step : 1 - c1 - c2 - c3 - c4
                    for c6 = 0 : step : 1 - c1 - c2 - c3 - c4 - c5
                        for c7 = 0 : step : 1 - c1 - c2 - c3 - c4 - c5 - c6
                            for c8 = 0 : step : 1 - c1 - c2 - c3 - c4 - c5 - c6 - c7
                                for c9 = 0 : step : 1 - c1 - c2 - c3 - c4 - c5 - c6 - c7 - c8
                                    % c20 = 1 - c1 - c2 - c3 - c4 - c5 - c6 - c7 - c6 - c9 - c10 - c11 - c12 - c13 - c14 - c15 - c16 - c17 - c18 - c19;
                                    p = p + 1;
                                    % mass_composition = [c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c14, c15, c16, c17, c18, c19, c20];
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
disp (['For 10 components varied in steps of ', num2str(step), ', number of possible mixtures is ', num2str(p)])
% disp (mass_composition)