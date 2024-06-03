function Cd = InterpCd(A_BEAVS,A_ref)
%% InterpCd()

% slope relationship between a ratio of diameters between a blade and the
% body
% a ratio of 1/3 gives a Cd of 0.23
% assume a linear relationship
% information sourced from textbook provided by Dr. Roberto Albertani
m =  1.6 / (1/3);

% Cd ~= m * d/D;  <-- ratio of diameters
% sqrt( A_BEAVS / A_ref ) = d/D
% therefore,

Cd = m.*sqrt(A_BEAVS ./ A_ref );

end