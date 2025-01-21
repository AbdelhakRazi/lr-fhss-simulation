function [PA_S_LoRa, PS_LoRa, PS_LoRa_Capture, Distance] = LoRa(SF, BW, CR, MonteCarlo)
    % Set default parameters that would normally come from parameters()
    Payload = 51;           % Payload size in bytes
    Nodes = 1000;          % Number of nodes
    Simulation_T = 3600;   % Simulation time in seconds
    pkct_p_h = 1;         % Packets per hour per node
    Pt = 0.025;           % Transmit power (25mW)
    wavelength = 0.34;     % Wavelength at 868MHz
    Gr = 1;               % Receiver antenna gain
    Gt = 1;               % Transmitter antenna gain
    eta = 2;              % Path loss exponent
    R = 6371e3;           % Earth radius in meters
    E = 10:1:90;               %Elevation Angles
    H = 780e3;            % Satellite altitude in meters
    
    %% Calculate Time on Air for standard LoRa
    n_preamble = 8;       % Number of preamble symbols
    
    % Calculate number of symbols
    n_payload = 8 + ceil((8*Payload - 4*SF + 28 + 16 - 20*0)/(4*(SF-2*0)))*(CR + 4);
    
    % Calculate symbol duration
    T_sym = (2^SF)/BW;
    
    % Total packet duration
    ToA_LoRa = (n_preamble + n_payload) * T_sym;
    
    %% Distance from user to satellite as function of elevation angle
    [Distance, Elevation_Angles, Ground_distance,FootPrint_R]= Satellite_Geometry(H,E);
    X = [1 8 11 14 17 19 22 25]; %To simulate fewer points
    Distance=Distance(X);
    E_angles = [10 20 30 40 50 60 70 80 90];
    K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
    k = sort(interp1(E_angles,K_factor,Elevation_Angles),'descend');
    % %% Distance calculation
    % % Calculate distances for elevation angles from E to 90 degrees
    % elevation_range = E:1:90;
    % Distance = H ./ sind(elevation_range);
    % Elevation_Angles = elevation_range;
    
    % % Sample specific points like in the original code
    % X = [1 8 11 14 17 19 22 25];
    % Distance = Distance(X);
    % Elevation_Angles = Elevation_Angles(X);
    
    %% Rician K-factor calculation
    E_angles = [10 20 30 40 50 60 70 80 90];
    K_factor = [1.24 3.07 3.24 3.6 3.89 5.63 9.77 17.06 25.11];
    k = sort(interp1(E_angles, K_factor, Elevation_Angles), 'descend');
    
    %% Main simulation loop
    PA_S_LoRa = zeros(1, length(Distance));
    PS_LoRa = zeros(1, length(Distance));
    PS_LoRa_Capture = zeros(1, length(Distance));
    
    for c = 1:length(Distance)
        N = Nodes;
        decoded = 0;
        decoded_Capture = 0;
        
        for m = 1:MonteCarlo
            % Generate random transmission time for target node
            Ts = rand(1,1) * Simulation_T;
            
            % Generate interfering transmissions using Poisson process
            lambda = N * pkct_p_h / 3600; % Convert packets per hour to packets per second
            num_interferers = poissrnd(lambda * 2 * ToA_LoRa);
            
            if num_interferers > 0
                % Generate random times for interferers within vulnerable period
                interferer_times = Ts + (rand(1, num_interferers) * 2 - 1) * ToA_LoRa;
                
                % Count overlapping packets
                overlaps = sum(abs(interferer_times - Ts) < ToA_LoRa);
                
                if overlaps > 0
                    % Calculate received powers for capture effect
                    % Target signal
                    kD = k(c);
                    muD = sqrt(kD/(2*(kD+1)));
                    sigmaD = sqrt(1/(2*(kD+1)));
                    h1D = abs(sigmaD*randn(1,1) + muD + 1j*(sigmaD*randn(1,1) + muD))^2;
                    Pr_desired = Pt*h1D*Gr*Gt*(wavelength/(4*pi*Distance(c)))^eta;
                    
                    % Interfering signals
                    total_interference = 0;
                    for i = 1:overlaps
                        % Random distance for interferer within footprint
                        r_int = sqrt(rand()) * Distance(c);
                        k_int = interp1(E_angles, K_factor, asind(H/sqrt(H^2 + r_int^2)));
                        
                        mu_int = sqrt(k_int/(2*(k_int+1)));
                        sigma_int = sqrt(1/(2*(k_int+1)));
                        h1_int = abs(sigma_int*randn(1,1) + mu_int + 1j*(sigma_int*randn(1,1) + mu_int))^2;
                        
                        Pr_int = Pt*h1_int*Gr*Gt*(wavelength/(4*pi*r_int))^eta;
                        total_interference = total_interference + Pr_int;
                    end
                    
                    % Check if packet survives with capture effect
                    if Pr_desired >= 6 * total_interference % 6 dB threshold for capture
                        decoded_Capture = decoded_Capture + 1;
                    end
                else
                    % No overlaps - packet succeeds
                    decoded = decoded + 1;
                    decoded_Capture = decoded_Capture + 1;
                end
            else
                % No interferers - packet succeeds
                decoded = decoded + 1;
                decoded_Capture = decoded_Capture + 1;
            end
        end
        
        % Calculate probabilities
        PS_LoRa(c) = decoded;
        PS_LoRa_Capture(c) = decoded_Capture;
        
        % Analytical probability (simplified collision model)
        lambda = N * pkct_p_h / 3600;
        PA_S_LoRa(c) = exp(-2 * lambda * ToA_LoRa);
    end
    
    % Normalize simulation results
    PS_LoRa = PS_LoRa / MonteCarlo;
    PS_LoRa_Capture = PS_LoRa_Capture / MonteCarlo;
end

