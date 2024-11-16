%% InterpCd()
function Cd = InterpCd(A_BEAVS,A_ref)
%% SUMMARY
%   This function is used to determine an approximate Cd for BEAVS's blades
%
%% INPUTS
%   A_BEAVS     double, surface area of all BEAVS blades, combined, m2
%   A_ref       double, reference area of the rocket, m2
%
%% OUTPUTS
%   Cd          double, coefficient of drag for BEAVS

% slope relationship between a ratio of diameters between a blade and the
% rocket body
% a ratio of 1/3 gives a Cd of 1.6
% assume a linear relationship
% information sourced from textbook provided by Dr. Roberto Albertani
m =  1.6 / (1/3);

% assume a linear relationship for Cd (it's not)
%            \/ ratio of diameters
% Cd =~ m * d/D
% d/D = sqrt( A_BEAVS / A_ref )
% therefore,

Cd = m.*sqrt(A_BEAVS ./ A_ref );

end