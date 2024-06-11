%% OSU ESRA 2024
% This script is intended to output sample data via serial to the BEAVS4
% system
% Dexter Carpenter

% This script calls upon BEAVS Simulation 5 to obtain a sample flight path
% Switch directory to BEAVS 5 Simulation
Folder = cd;
cd(fullfile(Folder, '..', '[5.8.5.5]BEAVS4_Sim5_R2020b'));

% Run Sim5
BEAVS4_Sim5();

% Switch back to SerialOutput
Folder = cd;
cd(fullfile(Folder, '..', '[5.8.5.6]BEAVS4_SerialOutput'));

%% Serial Communication Setup

fprintf('\n\n---------------------------------------------------------\n');
fprintf('BEAVS4_SerialOutput.m\n');
fprintf('Select Communication Port\n');
fprintf('Example: COM1\n');

% move cursor to command window
commandwindow;

% Setup Serial port
port = input('Serial Port: ', 's');
if isempty(port); port = 'COM1'; end % if no input, use COM1
fprintf('%s Selected\n',port);
BaudRate = 115200;
SerialPort = serialport(port, BaudRate);

%% Start Data Transfer

% Wait for user input to begin sending data
fprintf('\n\nPress any key to begin transferring data\n'); pause

snum = 4; % use simulation number 4 (PID)
idxprev = 0; % set prev value
tic; % begin stopwatch
while (toc < Time(end,snum)) % while still within the time vector
    timer = round(toc,3); % round the timer to 3rd decimal place
    idx = find(Time(:,snum)==timer); % get the index for given time
    if (idx==idxprev); continue; end % if same index as last, skip sending
    idxprev = idx; % update new prev. index
    if (~isempty(idx)) % if the index isn't empty
        WriteData(SerialPort,Time(idx,snum),h(idx,snum)); % send data to Serial
    end
    if (V(idx,4) < 0); break; end % if past apogee, stop
end
toc;
fprintf('\nData transfer complete! Yippie!\n');

%% Functions

function WriteData(SerialPort,T,H)
    % send serial data
    write(SerialPort,T,'uint8');
    %data = read(SerialPort,2,'uint32');
    % print data that was sent to cmd window
    fprintf('\n    T:%0.3fs    H:%0.3fft    ',T,H*3.28084);
end
