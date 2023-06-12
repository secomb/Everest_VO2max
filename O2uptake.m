function [Pa, Sa] = O2uptake(Pv, P50a, P50v, hilln, QT, DLeff, Palv, C0Hb, lp, qc)
% O2 uptake
    Sv = S(Pv, P50v, hilln);
    tlung = [0:0.01:1];
    [T, Plung] = ode45(@dPdt, tlung, Pv, [], P50a, P50v, hilln, QT, DLeff, Palv, C0Hb);
    Plung = real(Plung);
    if lp > 0
        % qcolor = {'r','k','b'};
        qcolor = flipud(copper(3));
        style = {':','-','--'};
        figure(2);
        % plot(T,P,'linestyle',style{lp},'color',qcolor{qc});
        plot(T,Plung,'linestyle',style{lp},'color',qcolor(qc,:));
        hold on;
    end
    Pa = real(Plung(end));
    Sa = S(Pa, P50a, hilln);
end

