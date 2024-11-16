%% InterpCdInverse()
function A_BEAVS = InterpCdInverse(Cd,A_ref)
%% SUMMARY
%   This function is used to backwards-calculate the area BEAVS needs to
%   produce in order to achieve a certain Cd. This is useful for
%   determining how far the blades need to extend in order to achieve a
%   certain braking force.
%   
%   See InterpCd() as well
%
%% INPUTS
%   Cd          double, desired Cd for BEAVS blades
%   A_ref       double, reference area of the rocket, from openrocket, m2
%
%% OUTPUTS
%   A_BEAVS     double, surface area that BEAVS must achieve between all
%                       it's blades, m2

% slope relationship between a ratio of diameters between a blade and the
% rocket body
% a ratio of 1/3 gives a Cd of 1.6
% assume a linear relationship
% information sourced from textbook provided by Dr. Roberto Albertani
m =  1.6 / (1/3);

% Cd ~= m * d/D;  <-- ratio of diameters
% sqrt( A_BEAVS / A_ref ) = d/D
% therefore,

A_BEAVS = A_ref.*(Cd./m).^2;

end