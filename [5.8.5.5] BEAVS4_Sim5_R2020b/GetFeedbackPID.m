%% GetFeedbackPID()
function [u, varargout] = GetFeedbackPID(Time, Kp, Ki, Kd, u, V, Vtarg)

% define delta time
dt = Time(2) - Time(1);

% find error
err = Vtarg - V;

% update control function
u = u + (Kp+Ki*dt+Kd/dt)*err(3) + abs((-Kp-2*Kd/dt)*err(2)) + (Kd/dt)*err(1);

% resolve varargout
varargout{1} = err(3);

end

