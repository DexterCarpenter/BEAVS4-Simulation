%% GetFeedback()
function [Cd,BladeExtn,BladeExtnDesire] = GetFeedbackPred(Time, h0, V0, g, rho, A_ref, m, alpha, BladeExtnCurrent, Cd_rocket, BladeWdth, BladeCnt, BladeExtnRate, BEAVSExtnMAX)
% this function takes the current state of the rocket parameters and
% calculates how much to increase/decrease the blade extension.
% INPUTS ------------------------------------------------------------------
    % Time = [CurrentTime NextTime]
    % h = Current Altitude
    % V = Current Velocity
% OUTPUTS -----------------------------------------------------------------

%% Constants

hTarg = 10000*0.3048; % 10k ft Target Apogee convert to meters

%% Find Desired Blade Extension
% use fzero() to obtain a desired Cd

% if initial altitude is within 90% of target altitude
% or vertical velocity is negative (in descent)
if (h0/hTarg >= 0.90) || (V0 < 0)
    % retract beavs
    BladeExtnDesire = 0;
else
    % write c1 and k as a function of Cd
    k  = @(Cd) rho*Cd*A_ref/2/m*sin(alpha);
    c1 = @(Cd) -1/sqrt(g)/sqrt(k(Cd))*atan(V0*sqrt(k(Cd))/sqrt(g));
    
    % use fzero() to determine desired Cd
    Cd_guess = 1;
    Cd_desire = fzero(@(Cd) h0 - 1/k(Cd)*log(cos(sqrt(g)*sqrt(k(Cd))*c1(Cd))) - hTarg, Cd_guess);
    
    % calculate desired area from desired Cd
    % use idea behind InterpCd()
%     ExtnDesire = (Cd_desire-Cd_rocket)*Aref/Cd_BEAVS/BladeWdth/BladeCnt;
    BladeExtnDesire = A_ref/BladeWdth/BladeCnt*( (Cd_desire-Cd_rocket)/4.8 )^(2/3);
end

%% Adjust Blade Extension

% if desired to be beyond extension
if BladeExtnDesire > BEAVSExtnMAX
    % set to max extension
    BladeExtn = BEAVSExtnMAX;
% check if unable to adjust extension to desired in time step
elseif abs(BladeExtnDesire-BladeExtnCurrent) > BladeExtnRate*(Time(2)-Time(1))
    % adjust extension by maximum amount in time step
    BladeExtn = BladeExtnCurrent + sign(BladeExtnDesire-BladeExtnCurrent)*BladeExtnRate*(Time(2)-Time(1));
% check if desired is below zero (cannot be negative)
elseif BladeExtnDesire < 0
    % set to zero
    BladeExtn = 0;
else
    % set to desired extension
    BladeExtn = BladeExtnDesire;
end

%% Calculate Cd

A_BEAVS = BladeExtn*BladeWdth*BladeCnt;
Cd_BEAVS = InterpCd(A_BEAVS,A_ref);
Cd = Cd_rocket + Cd_BEAVS*(A_BEAVS/A_ref);

end
