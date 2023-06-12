function S = S(P, P50, hilln)
% Hill saturation
    S = P^hilln/(P^hilln + P50^hilln);
end

