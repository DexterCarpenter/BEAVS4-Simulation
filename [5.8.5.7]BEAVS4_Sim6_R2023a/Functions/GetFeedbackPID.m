%% GetFeedbackPID()
function [u, varargout] = GetFeedbackPID(Time, Kp, Ki, Kd, u, V, Vtarg)
%% SUMMARY
%   The purpose of this function is to provide a package that encapsulates
%   the PID functionality at its core. It takes K values, the control
%   function, time, and velocities and outputs an updated control function
%   and optional error
%
%% INPUTS
%   Time        double, [current time, next time] in seconds
%   Kp          double, Proportional PID Constant
%   Ki          double, Integral PID Constant
%   Kd          double, Derivative PID Constant
%   u           double, current control function value
%   V           double, [prev prev vel, prev vel, current vel], m/s
%   VTarg       double, target velocity from lookup table, m/s
%
%% OUTPUTS
%   u           double, new control function value
%   varargout   optional output
%       err     double, current error between velocities

% define delta time
dt = Time(2) - Time(1);

% find error
% [prev prev err, prev err, current err]
err = Vtarg - V;

% update control function
u = u + (Kp+Ki*dt+Kd/dt)*err(3) + abs((-Kp-2*Kd/dt)*err(2)) + (Kd/dt)*err(1);

% resolve varargout
varargout{1} = err(3); % output only current error

end

