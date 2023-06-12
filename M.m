function M = M(M0,P0,P)
%PO2 dependent tissue consumption
% Michaelis-Menten kinetics
M = M0 * P / (P + P0);
end

