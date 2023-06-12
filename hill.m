function S = hill(n,P50,P)
%PO2 dependent tissue consumption
% Michaelis-Menten kinetics
S = (P/P50)^n/(1 + (P/P50)^n);
end
