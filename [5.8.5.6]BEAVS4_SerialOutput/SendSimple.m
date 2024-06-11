
SerialPort = serialport('COM5', 115200);

for i = 1:10
    write(SerialPort,i,'double');
    fprintf('\n%0.0f',i);
    pause;
end