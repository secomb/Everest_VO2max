function [Pv, Sv, Pmt] = O2delivery_Everest(dm, mc, Paref, Pa, P50a, P50v, hilln, QTcgs, VO2cgs, M0cgs, P0, C0Hb, NAcf)

% Tissue parameters

rc = 2.5E-4;    % cm
L = 0.05;       % cm
K = 9.4E-10;    % (cm^2/s).(cm^3_O2/cm^3/torr)
Sh = 2.5;
Kpl = 8.3E-10;          % (cm^2/s).(cm^3_O2/cm^3/torr)
Mt = pi * Kpl * Sh;     % (cm^2/s).(cm^3_O2/cm^3/torr)

m_skel = 10.8;
m_other = 50.0;
m_tot = 70.0;
m_nonconsuming = m_tot - (m_skel + m_other);

mc

if (mc)
    massarray = [m_skel m_other]    
    cptmax = length(massarray)
else
    massarray = [m_skel 0]
    cptmax = 1;
end
massfrac = massarray/sum(massarray)

if (mc)
    VO2rest = 250;                                          % mlO2/min
    VO2cgsrest = VO2rest/60                                % mlO2/s
    VO2cgs_other = VO2cgsrest * massfrac(2)
    VO2cgs_skel = VO2cgs - VO2cgs_other;
    VO2cgsarray = [VO2cgs_skel VO2cgs_other]
else
    VO2cgsarray = [VO2cgs 0]
end
VO2cgs
sum(VO2cgsarray)
VO2frac = VO2cgsarray/sum(VO2cgsarray)

if (mc)                                                 % Multiple compartment mode
    Qrest = QC(VO2rest);
    Qrestcgs = Qrest * 1000/60;
    Qrestcgsarray = Qrestcgs * massfrac;
    Qmcgs = QTcgs - Qrestcgsarray(2);
    Qcgsarray = [Qmcgs Qrestcgsarray(2)]
else
    Qcgsarray = [QTcgs 0]
end
Qcgsarray
flowfrac = Qcgsarray/sum(Qcgsarray)

% NAcfarray = [1288 304]                             % Roy 2019
NAcfarray = [NAcf 468]                               % NAcf = 937 from Brooke 2003 Ref to Schaffartzik
rtarray = 0.1./sqrt(pi.*NAcfarray)
rho = 1.06;                                          % g/cm^3
tissvolarray = 1E3 * massarray/rho                  % cm^3
perfarray = Qcgsarray./tissvolarray;
perfarray(isnan(perfarray)) = 0
Sa = S(Paref, P50a, hilln)
Svarray = Sa - (VO2cgsarray./Qcgsarray) / C0Hb;
Svarray(isnan(Svarray)) = 0
Sv = dot(Svarray,flowfrac)
% Sv = Sa - VO2cgs/(QTcgs * C0Hb)
Pvarray = [P(Svarray(1), P50a, hilln) P(Svarray(2), P50a, hilln)];

% Establish baseline demand M0 such that venous PO2 satisfies FickMM at sea level with normal P50
M0cgsarray = VO2cgsarray.*(Pvarray + P0)./Pvarray;                    % mlO2/s
M0cgsarray(isnan(M0cgsarray)) = 0
M0volcgsarray = M0cgsarray./tissvolarray;
M0volcgsarray(isnan(M0volcgsarray)) = 0
M0BrookeUnits = M0volcgsarray * 60 * 100

Varray = L* perfarray .* (rtarray./rc).^2;
Varray(isnan(Varray)) = 0
Qcaparray = Varray*pi*rc^2;                                              % Flow, cm^3/s
Qcaparray(isnan(Qcaparray)) = 0

switch dm
    case 1 % Fick
        % O2 delivery using Fick
        Sa = S(Pa, P50a, hilln);
        Sv = Sa - VO2cgs/(QTcgs * C0Hb);
        Sv = max(Sv,0.0001);
        Pv = P(Sv, P50v, hilln);
        Pmt = Pv;
    case 2 % FickMM
        Sa = S(Pa, P50a, hilln);
        if (mc)
            % Solve set of 4 nonlinear equations for x = {P1,P2,M1,M2}
            % Initial guesses
            x0 = [hillinv(hilln,P50a,0.25),hillinv(hilln,P50a,0.75),75,25]
            x = findFickMM2cpt(Pa,P50a,hilln,P0,Qcgsarray,C0Hb,VO2cgs,x0)
            Pvarray = x(1:2)
            Svarray = [S(Pvarray(1),P50a,hilln) S(Pvarray(2),P50a,hilln)]
            Sv = dot(Svarray,flowfrac)
            Pv = P(Sv,P50a,hilln)
            Pmt = Pv;
        else
            % Use Fick as initial guess for Pv
            Sa = S(Pa, P50a, hilln);
            Sv = Sa - VO2cgs/(QTcgs * C0Hb);
            Sv = max(Sv,0.0001);
            Pv0 = P(Sv, P50v, hilln);
%            Pv0 = 20;                               % Initial guess for Pv
            Pv0 = Pa                               % Initial guess for Pv
            Pv = findFickMMroot(Pa,P50a,P50v,hilln,QTcgs,C0Hb,M0cgs,P0,Pv0)
            if(isnan(Pv))
                Pa
                P50a
                QTcgs
                C0Hb
                M0cgs
                P0
                Pv0
                error('Pv isnan')
            end
            Sv = S(Pv, P50v, hilln)
            VO2MM = M(M0cgs,P0,Pv)
            VO2Fick = (QTcgs * C0Hb)*(Sa - Sv)
            Pmt = Pv;
        % Given PaO2 and estimated VO2
        end
        % O2delivery = Q * C0Hb * (S(PaO2) - S(PvO2))
        % Since M = M0 * Pv/(Pv+P0), calculate M0 from VO2*(Pv+P0)/Pv

%       Return Sv, Pv, Pmt

    case 3 % Krogh
        % Use P50 = P50(z), interpolated along cap between P50a and P50v
        % Calculate velocity v from Q
        % Calculate demand from VO2

        Pvarray = [0 0];
        Svarray = [0 0];
        Pmtarray = [0 0];
        for cpt = 1:cptmax
            [Pvarray(cpt), Svarray(cpt), Pmtarray(cpt)] = krogh_cpt(Pa, rtarray(cpt), rc, L, M0volcgsarray(cpt), Qcaparray(cpt), perfarray(cpt), P0, Mt, K, P50a, P50v, hilln, C0Hb)
        end


%         % Mb parameters
% 
%         DMb = 1.73E-7;  % cm^2/s
%         CMb = 3.83E-7;  % mol/cm^3
%         Vm = 2.2E4;     % cm^3
%         P50Mb = 3.2;    % torr
% 

%         tissvol = 2.3E3;                           % cm^3
%         V = 2.25e-1;                               % cm/s
%         V = QT/(pi * rc^2)
%         Perf = (pi * rc^2 * V)/(pi * rt^2 * L)     % cm^3/cm^3/s
%         M0 = 40 / (60 * 100);                      % cm^3_O2/cm^3/s


%        Perf = QTmcgs/tissvol;

%       Other compartment Krogh calculation

        % Find demand at rest
%         
% 
%         tissvol = 1E3 * (m_tot - m_skel)/rho;                  % cm^3
%         M0volcgs = M0cgsrest/tissvol;
%         Perf = QT0othercgs/tissvol;
%         V = Perf * L * (rt/rc)^2;
%         M0BrookeUnits = (M0cgsrest/tissvol) * 60 * 100;
%         Qcgs = V * pi * rc^2;                         % Flow, cm^3/s
%         [Pv_other, Sv_other, Pmt_other] = krogh_cpt(Pa, rt, rc, L, M0volcgs, Qcgs, Perf, P0, Mt, K, P50a, P50v, hilln, C0Hb);
% 

%       Return Sv, Pv, Pmt
        Sv = dot(Svarray,flowfrac)
        Pv = P(Sv,P50a,hilln)
        Pmt = dot(Pmtarray,massfrac)
    end
end
