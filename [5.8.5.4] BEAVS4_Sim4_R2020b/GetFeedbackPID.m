%% GetFeedbackPID()
function u = GetFeedbackPID(Time, Kp, Ki, Kd, u, h)

% 10k ft Target Apogee, convert to meters
h_targ = 10000*0.3048;

% define delta time
dt = Time(2) - Time(1);

% find error
err = h_targ - h;

% update control function
u = u + (Kp+Ki*dt+Kd/dt)*err(3) + (-Kp-2*Kd/dt)*err(2) + (Kd/dt)*err(1);

end

