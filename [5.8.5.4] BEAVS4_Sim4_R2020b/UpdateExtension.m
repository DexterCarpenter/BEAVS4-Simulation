%% UpdateExtension()
function [Cd,BladeExtn] = UpdateExtension(Time, A_ref, Cd_rocket, BladeExtnCurrent, BladeWdth, BladeCnt, BladeExtnRate, BEAVSExtnMAX, BladeExtnDesire)
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

