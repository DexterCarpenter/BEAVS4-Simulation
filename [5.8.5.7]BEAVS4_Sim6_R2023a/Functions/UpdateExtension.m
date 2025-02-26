%% UpdateExtension()
function [Cd,BladeExtn] = UpdateExtension(Time, A_ref, Cd_rocket, BladeExtnCurrent, BladeWdth, BladeCnt, BladeExtnRate, BEAVSExtnMAX, BladeExtnDesire)
%% SUMMARY
%   This function takes the current blade state along with a desired blade
%   extension. When using PID, the desired blade extension is controlled by
%   the control function, u. The function ensures that the blades cannot
%   extend faster than possible dictated by BladeExtnRate. It also ensures
%   the blades don't extend beyond their max. Outputs are the resulting
%   BEAVS Cd and extension.
%
%% INPUTS
%   Time                double, [prev time, next time] in sec
%   A_ref               double, rocket reference area, m2
%   Cd_rocket           double, Cd of the rocket from openrocket
%   BladeExtnCurrent    double, current extension of the blades, m
%   BladeWdth           double, width of blades, m
%   BladeCnt            double, number of blades
%   BladeExtnRate       double, max rate at which blades can extend, m/s
%   BEAVSExtnMAX        double, max distance blades can extend, m
%   BladeExtnDesire     double, desired extension, m
%
%% OUTPUTS
%   Cd                  double, resulting TOTAL Cd created
%   BladeExtn           double, resulting extension, m

%% Adjust Blade Extension

% if desired to be beyond extension max
if BladeExtnDesire > BEAVSExtnMAX
    % set to max extension
    BladeExtn = BEAVSExtnMAX;
% check if desired is below zero (cannot be negative)
elseif BladeExtnDesire < 0
    % set to zero
    BladeExtn = 0;
% check if unable to adjust extension to desired in time step
elseif abs(BladeExtnDesire-BladeExtnCurrent) > BladeExtnRate*(Time(2)-Time(1))
    % adjust extension by maximum amount in time step
    BladeExtn = BladeExtnCurrent + sign(BladeExtnDesire-BladeExtnCurrent)*BladeExtnRate*(Time(2)-Time(1));
else
    % set to desired extension
    BladeExtn = BladeExtnDesire;
end

%% Calculate Cd

A_BEAVS = BladeExtn*BladeWdth*BladeCnt; % calculate total surface area
Cd_BEAVS = InterpCd(A_BEAVS,A_ref); % interpolate the cd according to surface area
Cd = Cd_rocket + Cd_BEAVS*(A_BEAVS/A_ref); % add Cds

end

