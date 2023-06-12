  function [Pvk, Svk, Pmtk] = krogh_cpt(Pa, rt, rc, L, M0cgs, Q, Perf, P0, Mt, K, P50a, P50v, hilln, C0Hb) 

        consump = 0;                               % Initialize consumption
        % Boundary conditions

        Pb0 = Pa;      % torr (Inlet Pb at z = 0)

        % Integration parameters

        nz = 100;
        dz = L/nz;
        zarray = 0:dz:L;

        nr = 47;
        dr = (rt - rc) / nr;
        rarray = rc:dr:rt;

        maxit = 100000;
        Ptol = 1E-6;

        maxqit = 1000;
        qtol = 1E-14;

        iz = 1;
        % Estimate q (oxygen consumption per unit length) for first slice
        q0 = pi*(rt^2 - rc^2) * M0cgs;      % cm^3_O2/cm^3/s
        q(iz) = q0;
        
        % Initialize P(r,z) at z = 0
        Pb(iz) = Pb0;
        Pk(:,iz) = (Pb(iz) - q(iz)/Mt) * ones([nr + 1,1]);

        itreq = 0;
        for iz = 1:nz+1,iz;
            qit = 1;
            while (qit < maxqit)
                qprev = q(iz);
                it = 1;
                while (it < maxit)
                    Pzprev = Pk(:,iz);
                    Pk(1,iz) = Pb(iz) - q(iz)/Mt;
                    for ir = 2:nr
                        ri = rarray(ir);
                        rip = ri + 0.5 * dr;
                        rim = ri - 0.5 * dr;
                        Pk(ir,iz) = (0.5 / ri) * (rip * Pk(ir + 1,iz) + rim * Pk(ir - 1,iz)) - 0.5 * dr^2 * M(M0cgs,P0,Pk(ir,iz))/K;
                    end
                    Pk(nr+1,iz) = Pk(nr,iz) - 0.5 * dr^2 * M(M0cgs,P0,Pk(ir,iz))/K;
                    err(it) = max(abs(Pk(:,iz) - Pzprev));
                    if (err(it) < Ptol)
                        itreq = it;
                        it = maxit;
                    else
                        it = it + 1;
                    end
                end
                itreq;
                % recalculate q
                qtemp = 0;
                for ir = 1:nr
                    qtemp = qtemp + 2 * pi * rarray(ir) * dr * M(M0cgs,P0,Pk(ir,iz));
                end
                q(iz) = qtemp;
                qerr = abs(q(iz) - qprev);
                if (qerr < qtol)
                    qitreq = qit;
                    qit = maxqit;
                else
                    qit = qit + 1;
                end
            end
            consump = consump + q(iz) * dz;             % cm^3_O2/s
            qitreq;
            hold on;
            plot([0 rc rarray],[Pb(iz) Pb(iz) Pk(:,iz)']);
            Pk(:,iz+1) = Pk(:,iz);                        % Estimate of P in next slice
            Spb(iz) = hill(hilln,P50a,Pb(iz));
            Spbnew = Spb(iz) - q(iz) * dz / (Q * C0Hb);
            q(iz+1) = q(iz);                            % Estimate of q in next slice
            Pb(iz+1) = hillinv(hilln,P50a,Spbnew);
        end
        Pvk = Pb(end);
        Svk = hill(hilln,P50a,Pvk);
        Demand = M0cgs;
        Supply = C0Hb * Perf;                         % cm^3_O2/cm^3/s
        Consumption = consump/(pi*(rt^2-rc^2)*L);        % Volume specific consumption cm^3_O2/cm^3/s
        DeficitRatio = 1 - Consumption/Demand;
        Extraction = Consumption/Supply;
        Pk(:,nz+1) = [];                                 % Remove last column of P (not used)
        q(nz+1) = [];
        Pb(nz+1) = [];

        % Calculate mean PO2 in each slice
        nz;
        nr;
        Pksum = 0.0;
        Msum = 0.0;
        for iz = 1:nz+1
            Pkzsum = 0.0;
            Akzsum = 0.0;
            Mzsum = 0.0;
            for ir = 1:nr
                Akzsum = Akzsum + 2*pi * rarray(ir) * dr;
                Pkzsum = Pkzsum + 2*pi * rarray(ir) * dr * Pk(ir,iz);
                Mzsum = Mzsum + 2*pi * rarray(ir) * dr * M(M0cgs,P0,Pk(ir,iz));
            end
            Akz(iz) = Akzsum/(pi*(rt^2-rc^2));
            Pkz(iz) = Pkzsum/(pi*(rt^2-rc^2));
            Mz(iz) = Mzsum/(pi*(rt^2-rc^2));
            Pksum = Pksum + Pkz(iz);
            Msum = Msum + Mz(iz);
        end
        Pmtk = Pksum/nz;
        Mavg = Msum/nz;
        Mavg * L* pi*(rt^2-rc^2);
        consump;

%         figure
%         plot(zarray,Pkz);
%         hold on;
%         plot([0 L],[Pmtk Pmtk]);
%         hold off;
%         xlabel('Cap distance')
%         ylabel('Tissue PO_2 [Mean]')
% 
% 
%         figure
%         plot(zarray,Spb);
%         hold on;
%         xlabel('Cap distance')
%         ylabel('Tissue SO_2')
% 
%         figure
%         plot(zarray,Mz);
%         hold on;
%         plot([0 L],[Mavg Mavg]);
%         hold off;
%         xlabel('z')
%         ylabel('Consumption [Mean]')
% 
%         figure
%         plot(zarray,60 * 100 * q/(pi*(rt^2-rc^2)));
%         hold on;
%         plot([0 L],[60*100*Supply 60*100*Supply]);
%         hold on;
%         plot([0 L],[60*100*Demand 60*100*Demand]);
%         hold on;
%         plot([0 L],[60*100*Consumption 60*100*Consumption]);
%         legend ('Consumption','Supply','Demand','Mean Consumption')
%         
%         figure;
%         plot(zarray,Pb);
%         hold on;
%         plot(zarray,Pkz);
%         hold on;
%         plot(zarray,Pk(1,:));
%         hold on;
%         plot(zarray,Pk(nr,:));
%         legend({'$P_b$','$\overline{P}_t(slice)$','$P_t(wall)$','$P_t(edge)$'},'interpreter','latex')
%         % plot(log(err));
% 
%         figure;
%         imagesc(1e4*[0 L],1e4*[rc rt],Pk)
%         title('Tissue PO_2')
%         xlabel('z')
%         ylabel('r')
% 
%         figure;
%         [C,h] = contour(Pk);
%         clabel(C,h);
%         title('Tissue PO_2')
%         xlabel('z')
%         ylabel('r')
