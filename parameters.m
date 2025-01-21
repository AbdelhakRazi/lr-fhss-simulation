function [Payload,Header_N_DR8,Header_duration,F_duration,Header_ToA_DR8,Nodes,Simulation_T,pkct_p_h,Pt,wavelength,Gr,Gt,eta,E,R,H] = parameters
%% Gains and Pt are converted into linear form
Payload = 10;           % Message payload
Header_N_DR8 = 3;       % Header replicas
Header_duration = 0.233; %233 ms long headers
F_duration = 0.05;       %50 ms payload data fragments
Header_ToA_DR8 = Header_N_DR8*Header_duration;
%Nodes = [50 60 70 80 90 100].*1e3;
Nodes = 50e3;
Simulation_T = 3600; % 1 hour duration
pkct_p_h = 4;      % Packets per hour per end-device
Pt = 10^(14/10)/1000;      % Transmit Power of LoRa 14 dBm
Freq_Band = 868e6;         % 868 MHz (frequency band Europe)
SpeedLight  = 3e8;         % Speed of light
wavelength = SpeedLight/Freq_Band;

Gr=(10.^((22.6)/10));      %22.6 dBi: Gateway at Satellite
Gt=(10.^((2.15)/10));      %2.15 dBi: End-device
eta = 2;

%% parameters for satellites
E = 10:1:90;               %Elevation Angles
R = 6378e3;                % Radius of earth
H = 780e3;                 %Orbital height
end