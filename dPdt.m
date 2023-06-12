function dPdt = dPdt(tlung, P, P50a, P50v, hilln, QT, DLeff, Palv, C0Hb)
% dPdt is used in the calculation of lung oxygen uptake
%   Assumes constant P50
    P50 = P50v;                                                 % Can change to vary along lung capillary
    DLcgs = DLeff/60;
    dPdt = (DLcgs/QT) * (Palv - P) * ((P50^hilln + P^hilln)^2)/(C0Hb * P50^hilln * hilln * P^(hilln-1));
end

