function [Fb,Fn] = SolveFdrag(ma,mb,a,g,Fd)
%% SolveFdrag()
% INPUTS
    % ma | mass of top half of rocket (kg)
    % mb | mass of bottom half of rocket (kg)
    % a  | net acceleration of the rocket (m/s2)
    % g  | force of gravity (m/s2)
    % Fd | force of drag (N)
% OUTPUTS
    % Fb | force on a SINGULAR BEAVS blade
    % Fn | normal force through coupler

% Solve for known forces
Fnet = (ma+mb)*a;
FgA = ma*g;
FgB = mb*g;

% setup matrix
M = [ 0 -1 Fnet+Fd+FgA;...
     -2  1  Fnet+FgB   ];

% sovle system
F = rref(M);

% extract answers
Fb = F(1,end);
Fn = F(2,end);
end

