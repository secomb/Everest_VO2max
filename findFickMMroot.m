function FickMM = findFickMMroot(Pa,P50a,P50v,hilln,QTcgs,C0Hb,M0cgs,P0,Pv0)

FickMM = fzero(@VO2diff,[0 Pv0]);

    function FickMMdiff = VO2diff(Pv)
        
    % Given PaO2 and estimated VO2
    % O2delivery = Q * C0Hb * (S(PaO2) - S(PvO2))
    % VO2 = M0 * Pv/(Pv + P0)

    Sa = S(Pa, P50a, hilln);
    VO2MMcgs = M(M0cgs,P0,Pv);
    VO2Fickcgs = QTcgs * C0Hb * (Sa - S(Pv,P50v,hilln));
    FickMMdiff = VO2MMcgs - VO2Fickcgs;
    end
end

