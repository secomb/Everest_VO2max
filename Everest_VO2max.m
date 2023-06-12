
% Everest_VO2max 
% TKR, June 2023
% Simulation of oxygen uptake and utilization with altered hemoglobin-oxygen affinity at altitude

clc
clear all
close all
format compact
format shortG

% Set up cases for Hb and P50
limitHb = true; % Limit Hb to maximum summit Hb at normal P50

% -1: Hb fixed at normal Hb at sea level for all altitudes
%  0: Hb fixed at value corresponding to affinity type for all altitudes
%  1: Hb starts at value corresponding to affinity type and varies with altitude
varyHb = [-1 -1 -1; 0 0 0; 0 1 0; 1 1 1; 0 0 0; 0 1 0; 1 1 1];

%  0: Hb fixed at value corresponding to P50 for all altitudes
%  1: Hb starts at value corresponding to P50 and varies with altitude
varyP50 = [0 0 0; 0 0 0; 0 0 0; 0 0 0; 1 1 1; 1 1 1; 1 1 1];

Patm = PB(0);                   % Atmospheric pressure at sea level
Pw = 47;
FiO2 = 0.21;
RQ = 0.8;
P0 = 10.5;                      % Effective Km for Michaelis-Menten

VO2rest = 250;                                       % mlO2/min
VO2max = 2750;                                       % mlO2/min

VO2 = VO2max;
QT = QC(VO2);                                        % L/min, Calbet et al

summitalt = 8.848;

altarray = 0.001 * [75 5300 6400 7100 8400];         % Grocott altitude values
PBalt = PB(altarray);
PiO2predalt = FiO2 * (PBalt - Pw);

RQalt = RQ * ones(size(altarray))
PaO2measalt = (455 - [324 390 399 409 422]) * 250/(455 - 120) % Grocott Fig 2
SaO2measalt = 0.01 * (455 - [324 339 344 355 383]) * 250/(455 - 120) 

Hbmeasalt = 0.1 * (455 - [259 215 189 203 196]) * 250/(455 - 120) 
CaO2measalt = 1e-3 * (455 - [186 165 147 191 260]) * 250/(455 - 120)

PAO2measalt = PaO2measalt + 5.41
PaCO2measalt = [36.6 20.4 18.2 16.7 13.3];
pHmeasalt = [7.40 7.46 7.51 7.53 7.53];

Pfit = [29.5 19.1 21.0 28.7];
Sfit = 0.01 * [68.1 34.4 43.7 69.7];

% Fit Grocott summit P and S values to find P50 with fixed n = 2.7
hilln = 2.7;
P50ref = 26.3;    
HillF = @(coeff,Pdata)realpow(Pdata,hilln)/(realpow(Pdata,hilln)+realpow(coeff(1),hilln));
coeff0 = [P50ref]; % Initial guess
Fsumsquares = @(coeff)sum((HillF(coeff,Pfit) - Sfit).^2);
opts = optimoptions('fminunc','Algorithm','quasi-newton');
[xunc,ressquared,eflag,outputu] = fminunc(Fsumsquares,coeff0,opts);
fixednP50fit = xunc

plot(Pfit,Sfit,'o');
hold on;
for Ptest = 1:100
    Prange(Ptest) = Ptest;
    Srange(Ptest) = HillF(fixednP50fit,Ptest);
end
plot(Prange,Srange,'-')

sealevelindex = 1;
summitindex = length(altarray);
P50coeff = polyfit([altarray(sealevelindex) altarray(summitindex)],[P50ref fixednP50fit],1)
P50predalt = polyval(P50coeff,altarray)

Hbcoeff = polyfit(altarray,Hbmeasalt,1)
Hbmeasalt
Hbpredalt = polyval(Hbcoeff,altarray,1)

PaO2coeff = polyfit(altarray,PaO2measalt,1);
PaO2measalt

PaCO2coeff = polyfit(altarray,PaCO2measalt,1)
PaCO2measalt
PaCO2predalt = polyval(PaCO2coeff,altarray)

Palvtrunc = PiO2predalt - PaCO2predalt/RQ;
PAO2predalt = PiO2predalt - PaCO2predalt/RQ + PaCO2predalt * FiO2 * ((1-RQ)/RQ)

pHcoeff = polyfit(altarray,pHmeasalt,1)

altsamp = 47;               
minalt = 0;
maxalt = 10;
altrange = linspace(minalt,maxalt,altsamp);
summit = length(altrange);

RQaltrange = linspace(RQ,RQ,altsamp);
pHrange = polyval(pHcoeff,altrange);
PaCO2range = polyval(PaCO2coeff,altrange);
figure
subplot(2,2,1)
plot(altarray,PaCO2measalt,'o')
hold on
plot(altrange,PaCO2range)
xlabel('Altitude (km)');
ylabel('PaCO_2');
subplot(2,2,2)
plot(altarray,pHmeasalt,'o')
hold on
plot(altrange,pHrange)
xlabel('Altitude (km)');
ylabel('pH');
subplot(2,2,3)
PBaltspace = PB(altrange);
plot(altrange,PBaltspace);
xlabel('Altitude (km)');
ylabel('Atmospheric pressure (mmHg)');
subplot(2,2,4)
PiO2range = FiO2 * (PBaltspace - Pw);
PAO2range = PiO2range - PaCO2range./RQaltrange + FiO2 * PaCO2range .* ((1-RQaltrange)./RQaltrange);
plot(altrange,PAO2range);
hold on;
plot(altarray,Palvtrunc,'*')
hold on;
plot(altarray,PAO2predalt,'^')
xlabel('Altitude (km)');
ylabel('PAO_2')

% Wagner 2010 High Alt Med Biol Figure 6
figure
PiO2predalt = FiO2 * (PBalt - Pw);
PAO2predalt = PiO2predalt - PaCO2predalt./RQalt + FiO2 * PaCO2predalt .* ((1-RQalt)./RQalt);
plot(PAO2predalt,PaCO2predalt,'o')
xlabel('Calculated PAO_2');
ylabel('Measured PACO_2')
title('Wagner 2010 High Alt Med Biol Figure 6')

DLeff = 74; 

P50_LS = 15.6;                                  % Dominelli 2020
P50array = [P50_LS P50ref P50ref + (P50ref - P50_LS)]
Hbcalcarray = JLOHb(P50array)
Hbcalcarray(2) = Hbpredalt(sealevelindex)

Hbval = zeros(length(altrange),length(P50array));
P50val = zeros(length(altrange),length(P50array));

PBalt = PB(altrange);
PiO2predalt = FiO2 * (PBalt - Pw);
Hbpredalt = polyval(Hbcoeff,altrange,1);
P50predalt = polyval(P50coeff,altrange);
PaO2predalt = polyval(PaO2coeff,altrange);
PaCO2predalt = polyval(PaCO2coeff,altrange);
Palvtrunc = PiO2predalt - PaCO2predalt/RQ;
PAO2predalt = PiO2predalt - PaCO2predalt/RQ + PaCO2predalt * FiO2 * ((1-RQ)/RQ);

for P50index = 1:3
    for altindex = 1:length(altrange)
        Hbfix(altindex,P50index) = Hbpredalt(sealevelindex);
        Hbcon(altindex,P50index) = Hbcalcarray(P50index);
        Hbvar(altindex,P50index) = Hbcalcarray(P50index) * Hbpredalt(altindex)/Hbpredalt(sealevelindex);
        P50con(altindex,P50index) = P50array(P50index);
        P50var(altindex,P50index) = P50array(P50index) * P50predalt(altindex)/P50predalt(sealevelindex);
    end
    if (limitHb)
        if (Hbcalcarray(P50index) * Hbpredalt(summit)/Hbpredalt(sealevelindex) > Hbpredalt(summit))
            %(Hbpredalt(summit) - Hbpredalt(sealevel))/(altrange(end) - altrange(sealevel))
            newslope = (Hbpredalt(summit) - Hbcalcarray(P50index))/(altrange(end) - 0);
            newintcp = Hbcalcarray(P50index) - newslope * altrange(sealevelindex);
            Hbadjcoeff = [newslope newintcp]
            Hbvar(:,P50index) = polyval(Hbadjcoeff,altrange);
        end
    end
end
    
% Comparison of FickMM to Krogh

% Find demand at VO2max

VO2 = VO2max;                                         % mlO2/min
VO2cgs = VO2/60                                      % mlO2/s
QT = QC(VO2);                                            % L/min
QTcgs = QT * 1000/60                                    % ml/s
% Establish baseline demand M0 such that venous PO2 satisfies FickMM at sea level with normal P50
Paref = PaO2predalt(sealevelindex)
Sa = S(Paref, fixednP50fit, hilln)
C0Hb = 1.34 * Hbpredalt(sealevelindex) * 0.01;                         % mlO2/ml
Sv = Sa - VO2cgs/(QTcgs * C0Hb)
Pv = P(Sv, fixednP50fit, hilln)
M0cgs = VO2cgs*(Pv + P0)/Pv                     % mlO2/s

mc = false;                                     % Use multiple compartment model
PaO2range = 25:5:100;
Pvcfcasedm = [2 3 3 3];
PvcfcaseNacf = [0 1468 1100 700];
maxcase = length(Pvcfcasedm);

for P50case = 2:2
    P50val = P50array(P50case);
    for Pvcfcase = 1:maxcase
        dm = Pvcfcasedm(Pvcfcase)
        NAcf = PvcfcaseNacf(Pvcfcase)
        for i = 1:length(PaO2range)
            Pa = PaO2range(i);
            [PvO2, SvO2, PmtO2] = O2delivery_Everest(dm, mc, Paref, Pa, P50val, P50val, hilln, QTcgs, VO2cgs, M0cgs, P0, C0Hb, NAcf);
            SaO2 = S(Pa, P50val, hilln);
            extrac = (SaO2 - SvO2)/SaO2;
            VO2 = 60 * QTcgs * C0Hb * (SaO2 - SvO2);
            extraccf(Pvcfcase,i) = extrac;
            VO2cf(Pvcfcase, i) = VO2;
            Pvcf(Pvcfcase,i) = PvO2;
            Pmtcf(Pvcfcase,i) = PmtO2;
        end
    end
    figure
    plot(PaO2range,Pvcf)
    hold on;
    xlabel('P_a (mmHg)')
    ylabel('P_v (mmHg)')
    legend({'FickMM','Krogh 1468 mm^{-2}','Krogh 1100 mm^{-2}','Krogh 700 mm^{-2}'},'Location','southeast')
end

for dm = 2:2                                                % Delivery method 1 = Fick, 2 = FickMM, 3 = Krogh
    VO2 = VO2max;                                           % mlO2/min
    VO2cgs = VO2/60;                                        % mlO2/s
    QT = QC(VO2);                                           % L/min
    QTcgs = QT * 1000/60;                                   % ml/s
    % Establish baseline demand M0 such that venous PO2 satisfies FickMM at sea level with normal P50
    Pa = PaO2predalt(sealevelindex);
    Sa = S(Pa, fixednP50fit, hilln);
    C0Hb = 1.34 * Hbpredalt(sealevelindex) * 0.01;                         % mlO2/ml
    Sv = Sa - VO2cgs/(QTcgs * C0Hb);
    Pv = P(Sv, fixednP50fit, hilln);
    
    switch dm
        case 1
            M0cgs = VO2cgs;
        case 2
            % Establish baseline demand M0 such that venous PO2 satisfies FickMM at sea level
            M0cgs = VO2cgs*(Pv + P0)/Pv;                     % mlO2/s
        case 3
            % Calculate M0 using Krogh such that VO2 = VO2max at sea level using P50fit
            NAcf = 937;
            M0cgs = VO2cgs*(Pv + P0)/Pv;                     % mlO2/s Initial estimate
            tol = 1e-4;
            err = 1e+4;
            itmax = 10000;
            it = 0;
            lp = 0;
            qc = 0;
            while ((abs(err) > tol) && (it < itmax))
                M0cgsold = M0cgs;
                [PvO2, SvO2, PmtO2] = O2delivery_Everest(dm, oc, Paref, Pa, fixednP50fit, fixednP50fit, hilln, QTcgs, VO2cgs, M0cgs, P0, C0Hb, NAcf);
                VO2calc = 60 * (QTcgs * C0Hb)*(Sa - SvO2);
                % fprintf('%f <<< %f\n',PvO2,PaO2);
                err = (VO2calc - VO2max)/VO2max;
                M0cgs = M0cgs * (1.0 - err);
                it = it + 1;
            end
    end
    dm
    M0cgs
    %  Plot sample VO2MM and VO2Fick curves
    PaO2 = 100;
    Pvarray = linspace(0,PaO2);
    VO2Fickarray = 60 * QTcgs * C0Hb * (SaO2 - Pvarray.^hilln./(Pvarray.^hilln + fixednP50fit^hilln));
    VO2MMarray = 60 * M0cgs * Pvarray./(Pvarray + P0);
    figure
    hold on;
    plot(Pvarray,VO2Fickarray);
    hold on;
    plot(Pvarray,VO2MMarray,'--');
    xlabel('P_v (mmHg)')
    ylabel('Oxygen Utilization (mlO_2/min)')
    axis([0 100 0 4500])
    box on
    Pv0 = 20;
    Pvroot = findFickMMroot(Pa,fixednP50fit,fixednP50fit,hilln,QTcgs,C0Hb,M0cgs,P0,Pv0)
%     error('stop')
    
    PaO2array = zeros([size(varyHb),length(altrange)]);
    PvO2array = zeros([size(varyHb),length(altrange)]);
    SaO2array = zeros([size(varyHb),length(altrange)]);
    SvO2array = zeros([size(varyHb),length(altrange)]);
    CaO2array = zeros([size(varyHb),length(altrange)]);
    CvO2array = zeros([size(varyHb),length(altrange)]);
    extraction = zeros([size(varyHb),length(altrange)]);
    VO2array = zeros([size(varyHb),length(altrange)]);
    C0Hbarray = zeros([size(varyHb),length(altrange)]);
    Hbpredarray = zeros([size(varyHb),length(altrange)]);
    QTcgsarray = zeros([size(varyHb),length(altrange)]);
    casemax = 7;
    P50casemax = 3;
    
    for casenum = 1:casemax
        for P50index = 1:P50casemax
            for altindex = 1:length(altrange)
                switch varyP50(casenum,P50index)
                    case 0
                        P50a = P50con(altindex,P50index)
                    case 1
                        P50a = P50var(altindex,P50index)
                end
                P50v = P50a;
                switch varyHb(casenum,P50index)
                    case -1
                        Hbpred = Hbfix(altindex,P50index)
                    case 0
                        Hbpred = Hbcon(altindex,P50index)
                    case 1
                        Hbpred = Hbvar(altindex,P50index)
                end
                altindex
                C0Hb = 1.34 * Hbpred * 0.01;                         % mlO2/ml
                Sv = Sa - VO2cgs/(QTcgs * C0Hb);
                Pv = P(Sv, P50ref, hilln);
                PvO2 = Pv;                     % Initial estimate
                Palv = PAO2predalt(altindex);

                tol = 1e-4;
                err = 1e+4;
                itmax = 10000;
                it = 0;
                lp = 0;
                qc = 0;
                while ((err > tol) && (it < itmax))
                    PvO2old = PvO2;
                    [PaO2, SaO2] = O2uptake(PvO2, P50a, P50v, hilln, QTcgs, DLeff, Palv, C0Hb, lp, qc);
                    fprintf('%f >>> %f\n',PvO2,PaO2);
                    [PvO2, SvO2, PmtO2] = O2delivery_Everest(dm, mc, Paref, PaO2, P50a, P50v, hilln, QTcgs, VO2cgs, M0cgs, P0, C0Hb, NAcf);
                    fprintf('%f <<< %f\n',PvO2,PaO2);
                    err = abs(PvO2 - PvO2old)/PvO2;
                    it = it + 1;
                end
                CaO2 = C0Hb * SaO2;
                CvO2 = C0Hb * SvO2;
                extrac = (CaO2 - CvO2)/CaO2;
                VO2calc = 60 * (QTcgs * C0Hb)*(SaO2 - SvO2);
                if (isnan(VO2calc))
                    QTcgs
                    C0Hb
                    SaO2
                    SvO2
                    error('VO2calc isnan')
                end
                SaO2array(casenum, P50index, altindex) = SaO2;
                SvO2array(casenum, P50index, altindex) = SvO2;
                PaO2array(casenum, P50index, altindex) = PaO2;
                PvO2array(casenum, P50index, altindex) = PvO2;
                CaO2array(casenum, P50index, altindex) = CaO2;
                CvO2array(casenum, P50index, altindex) = CvO2;
                VO2array(casenum, P50index, altindex) = VO2calc;
                C0Hbarray(casenum, P50index, altindex) = C0Hb;
                Hbpredarray(casenum, P50index, altindex) = Hbpred;
                QTcgsarray(casenum, P50index, altindex) = QTcgs;
                extraction(casenum, P50index, altindex) = extrac;
                switch dm
                    case 2
                        %  Plot VO2MM and VO2Fick curves
                        Pvarray = linspace(0,PaO2);
                        VO2Fickarray = 60 * QTcgs * C0Hb * (SaO2 - Pvarray.^hilln./(Pvarray.^hilln + P50v^hilln));
                        VO2MMarray = 60 * M0cgs * Pvarray./(Pvarray + P0);
                        figure(10)
                        hold on;
                        plot(Pvarray,VO2Fickarray);
                        hold on;
                        plot(Pvarray,VO2MMarray,'--');
                        hold on;
                        xlabel('PvO_2')
                        ylabel('Oxygen Consumption (mlO_2/min)')
                end
            end
        end
    end
    
    testarray = VO2array - 60*QTcgsarray.*C0Hbarray.*(SaO2array-SvO2array);
    tol = 1e-12;
    zerotest = any(abs(testarray(:)) < tol)                                % B = any(A(:) < 0)
    stylearray = ["-","-","-","-","-","-","-"];
    rgb = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];
    figure
    subplot(2,2,1)
    for casenum = 1:casemax
        for P50index = 1:P50casemax
            set(gca,'ColorOrderIndex',P50index)
            plot(altrange,squeeze(PaO2array(casenum,P50index,:)),stylearray(casenum))
            hold on;            
            xticks([0:maxalt]);
            axis([0 maxalt 0 110]);
            ymax = ylim
            plot([summitalt summitalt],[0 ymax(2)],'--');
            hold on;
        end
    end
    leg = legend({'Low P_{50}','','Normal P_{50}','','High P_{50}'},'Location','northeast','FontSize',8)
    xlabel('Altitude (km)');
    ylabel('Arterial oxygen tension (mmHg)');

    subplot(2,2,2)
    for casenum = 1:casemax
        for P50index = 1:P50casemax
            set(gca,'ColorOrderIndex',P50index)
            plot(altrange,squeeze(100*SaO2array(casenum,P50index,:)),stylearray(casenum))
            axis([0 maxalt 0 100]);
            hold on
            ymax = ylim;
            plot([summitalt summitalt],[0 ymax(2)],'--');
            hold on;
        end
    end
    
    hold on
    xticks([0:maxalt]);
    axis([0 maxalt 0 inf]);
    ymax = ylim
    plot([summitalt summitalt],[0 ymax(2)],'--');
    xlabel('Altitude (km)');
    ylabel('Arterial blood saturation (%)');
    
    subplot(2,2,3)
    for casenum = 1:casemax
        for P50index = 1:P50casemax
            set(gca,'ColorOrderIndex',P50index)
            plot(altrange,squeeze(100*CaO2array(casenum,P50index,:)),stylearray(casenum))
            hold on
        end
    end
    
    hold on
    xticks([0:maxalt]);
    axis([0 maxalt 0 30]);
    ymax = ylim
    plot([summitalt summitalt],[0 ymax(2)],'--');
    
    xlabel('Altitude (km)');
    ylabel('Arterial oxygen content (ml O_2/dL)');
    
    subplot(2,2,4)
    for casenum = 1:casemax
        for P50index = 1:P50casemax
            set(gca,'ColorOrderIndex',P50index)
            plot(altrange,squeeze(100*extraction(casenum,P50index,:)),stylearray(casenum))
            hold on
        end
    end
    xlabel('Altitude (km)');
    ylabel('Systemic oxygen extraction (%)');
    xticks([0:maxalt]);
    axis([0 maxalt 50 100]);
    hold on
    ymax = ylim;
    plot([summitalt summitalt],[0 ymax(2)],'--');

    figure
    for casenum = 1:casemax
        for P50index = 1:P50casemax
            set(gca,'ColorOrderIndex',P50index)
            plot(altrange,squeeze(VO2array(casenum,P50index,:))/1e3,stylearray(casenum))
            hold on;
        end
    end
    xlabel('Altitude (km)');
    ylabel('VO_{2}max (L/min)');
    
  
    for P50index = P50casemax:-1:1
        mindata = min([VO2array(:,P50index,:)],[],1)/1e3;
        maxdata = max([VO2array(:,P50index,:)],[],1)/1e3;
        pln(P50index) = fill([altrange';flipud(altrange')],[squeeze(mindata);flipud(squeeze(maxdata))],rgb(P50index,:))
    end
    pln
    leg = legend(pln,{'Low P_{50}','Normal P_{50}','High P_{50}'},'Location','northeast','FontSize',8)
end
leg.AutoUpdate = 'off';
hold on;
ymax = ylim;
plot([summitalt summitalt],[0 ymax(2)],'--');
figure(20)
hold on;
plot(Pvarray,VO2Fickarray);
hold on;
plot(Pvarray,VO2MMarray,'--');
hold on;
xlabel('P_v (mmHg)')
ylabel('Oxygen Utilization (mlO_2/min)')

rc = 2.5E-4;    % cm
L = 0.05;       % cm
K = 9.4E-10;    % (cm^2/s).(mlO2/cm^3/mmHg)
Sh = 2.5;
Kpl = 8.3E-10;          % (cm^2/s).(mlO2/cm^3/mmHg)
Mt = pi * Kpl * Sh;     % (cm^2/s).(mlO2/cm^3/mmHg)

m_skel = 10.8e3;                % g
rho = 1.06;                     % g/cm^3
V_skel = m_skel/rho;            % cm^3

NAcf = 1468;                    % Capillary density, mm^-2
rt = 0.1/sqrt(pi*NAcf);         % cm

P50ref = 26.3;
hilln = 2.7;

Prefarray = [PaO2array(2,:,1); PaO2array(2,:,end)]
altnum = size(Prefarray,1)
Pvrefarray = [PvO2array(2,:,1); PvO2array(2,:,end)]
Hbpredarray = [Hbpredarray(2,:,1); Hbpredarray(2,:,end)]
VO2calcarray = [VO2array(2,:,1); VO2array(2,:,end)]
VO2rest = 250;         % mlO2/min
term = rt^2 - rc^2;
PbKroghslope = term/(Kpl*Sh) + (term + 2*rt^2*log(rt/rc))/(4*K)
VO2index = 0;
for VO2fac = 0:0.1:15
    VO2index = VO2index + 1;
    VO2 = VO2rest * VO2fac;
    Q = QC(VO2);
    VO2cgs = VO2/60;
    PbKrogh = (VO2cgs/V_skel) * PbKroghslope;
    VO2FickMMarray(VO2index) = VO2;
    % Qarray(VO2index) = Q;
    Qarray(VO2index) = 315.1*60/1000;
    PbKrogharray(VO2index) = PbKrogh;
end
Prefstr = strings(altnum,1);
figure;
for altindex = 1:altnum
    for P50index = 1:length(P50array)
        P50 = P50array(P50index);
        Sa = S(Prefarray(altindex,P50index),P50,hilln);
        for VO2index = 1:length(VO2FickMMarray)
            VO2 = VO2FickMMarray(VO2index);
            VO2cgs = VO2/60;
            Q = Qarray(VO2index);
            Qcgs = Q*1000/60;
            C0Hb = 1.34 * Hbpredarray(altindex,P50index) * 0.01;
            Sv = Sa - VO2cgs/(Qcgs*C0Hb);
            PbFick = P(Sv,P50,hilln);
            Qarray(VO2index) = Q;
            PbFickarray(P50index,VO2index) = PbFick;       
            Prefstr(altindex,P50index) = ['P_a = ' num2str(floor(Prefarray(altindex,P50index))) '; P_{50} = '  num2str(floor(P50))];
        end
        VO2maxarray(P50index) = interp1(PbFickarray(P50index,:)-PbKrogharray,VO2FickMMarray,0,'spline');
    end
    plot(PbKrogharray,VO2FickMMarray/1e3,':')
    hold on;
    if altindex == 1
        plot(PbFickarray,VO2FickMMarray/1e3,'-')
    else
        plot(PbFickarray,VO2FickMMarray/1e3,'--')
    end
    xlabel('P_v (mmHg)')
    ylabel('VO_2 (L/min)');
end
% Prefstr(altindex) = strcat('Pa = ',num2str(floor(Prefarray(altindex))))
nullstr = "";
warning('off','MATLAB:legend:IgnoringExtraEntries');
Prefstr = [[nullstr nullstr]' Prefstr];
legend([reshape(Prefstr',[1 8])]')
legend('AutoUpdate','off')
for altindex = 1:altnum
    plot(Pvrefarray(altindex,:),VO2calcarray(altindex,:)/1e3,'o')
end
plot([0 0],[0 4],'k-')
