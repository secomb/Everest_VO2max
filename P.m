function P = P(S, P50, hilln)
% Inverse Hill equation
    % P = P50 * (S/(1-S))^(1/hilln);
    if S > 0
        P = P50 * nthroot((S/(1-S)),hilln);
    else
        P = 0;
    end
end

