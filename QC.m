function QC = QC(VO2)
% Calbet et al 2014
    Q_TD = 4.37 + 5.33 * (VO2/1000);
    Q_ICG = 4.43 + 5.22 * (VO2/1000);
    QC = 0.5 * (Q_TD + Q_ICG);
end

