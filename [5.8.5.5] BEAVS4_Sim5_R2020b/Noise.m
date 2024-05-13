%% Noise()
function [varout] = Noise(varin,degree)
%% INPUTS
    % varin  | 2 x 1 vector containing variable to be adjusted
    %        | varin(1) is previous value, varin(2) is value to be adjusted
    % degree | value from 0-100 corresponding to a precent.
    %        | adjustes the severity of the noise
%% OUTPUTS
    % varout | adjusted value, scalar

stdev = abs(varin(2) - varin(1))*degree/100;
varout = normrnd(varin(2),stdev);

end