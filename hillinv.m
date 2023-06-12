function P = hillinv(n,P50,S)
%PO2 dependent tissue consumption
% Michaelis-Menten kinetics
P = P50 * (S/(1-S))^(1/n);
end
