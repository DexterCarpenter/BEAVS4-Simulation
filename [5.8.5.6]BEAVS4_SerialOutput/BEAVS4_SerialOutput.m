%% OSU ESRA 2024
% This script is intended to output sample data via serial to the BEAVS4
% system
% Dexter Carpenter

% This script calls upon BEAVS Simulation 5 to obtain a sample flight path
% Switch directory to BEAVS 5 Simulation
Folder = cd;
cd(fullfile(Folder, '..', '[5.8.5.5]BEAVS4_Sim5_R2020b'));

% run Sim5
BEAVS4_Sim5();

% Switch back to SerialOutput
Folder = cd;
cd(fullfile(Folder, '..', '[5.8.5.6]BEAVS4_SerialOutput'));



