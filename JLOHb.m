function Hbest = JLOHb(P50)
% JLO Correlation Hb = c1 * P50 + c0
% Shepherd et al. 2019, J Physiol 597:4193
    coeff = [-0.3135 23.636];
    Hbest = polyval(coeff,P50);
end

