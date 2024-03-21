function [Cd,Extn,ExtnDesire] = GetFeedback(Time, h0, V0, g, rho, Aref, m, alpha, ExtnCurrent, Cd_rocket, BladeWdth, BladeCnt, BladeExtnRate, BEAVSExtnMAX)
%% GetFeedback()
% this function takes the current state of the rocket parameters and
% calculates how much to increase/decrease the blade extension.
% INPUTS ------------------------------------------------------------------
    % Time = [CurrentTime NextTime]
    % h = Current Altitude
    % V = Current Velocity
% OUTPUTS -----------------------------------------------------------------

%% Constants

hTarg = 10000*0.3048; % 10k ft Target Apogee covert to meters

%% Desired Blade Extension
% use fzero() to obtain a desired Cd

% if initial altitude is within 99% of target altitude
% or vertical velocity is negative (in descent)
if (h0/hTarg >= 0.90) || (V0 < 0)
    % retract beavs
    ExtnDesire = 0;
else
    % write c1 and k as a function of Cd
    k  = @(Cd) rho*Cd*Aref/2/m*sin(alpha);
    c1 = @(Cd) -1/sqrt(g)/sqrt(k(Cd))*atan(V0*sqrt(k(Cd))/sqrt(g));
    
    % use fzero() to determine desired Cd
    Cd_guess = 1;
    Cd_desire = fzero(@(Cd) h0 - 1/k(Cd)*log(cos(sqrt(g)*sqrt(k(Cd))*c1(Cd))) - hTarg, Cd_guess);
    
    % calculate desired area from desired Cd
    % use idea behind InterpCd()
%     ExtnDesire = (Cd_desire-Cd_rocket)*Aref/Cd_BEAVS/BladeWdth/BladeCnt;
    ExtnDesire = Aref/BladeWdth/BladeCnt*( (Cd_desire-Cd_rocket)/4.8 )^(2/3);
end

%% Adjust Blade Extension

% if desired to be beyond extension
if ExtnDesire > BEAVSExtnMAX
    % set to max extension
    Extn = BEAVSExtnMAX;
% check if unable to adjust extension to desired in time step
elseif abs(ExtnDesire-ExtnCurrent) > BladeExtnRate*(Time(2)-Time(1))
    % adjust extension by maximum amount in time step
    Extn = ExtnCurrent + sign(ExtnDesire-ExtnCurrent)*BladeExtnRate*(Time(2)-Time(1));
% check if desired is below zero (cannot be negative)
elseif ExtnDesire < 0
    % set to zero
    Extn = 0;
else
    % set to desired extension
    Extn = ExtnDesire;
end

% calculate Cd
ABEAVS = Extn*BladeWdth*BladeCnt;
Cd_BEAVS = InterpCd(ABEAVS,Aref);
Cd = Cd_rocket + Cd_BEAVS*(ABEAVS/Aref);

end

