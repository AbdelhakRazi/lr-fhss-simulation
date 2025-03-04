function [Distance, Elevation_Angles, Ground_distance,FootPrint_R] = Satellite_Geometry (H,E)

     %% Maximum distance from user to satellite 
    R = 6378e3;            %Radius of earth 

    X = cosd(E).*cosd(E);  %cos^2
    V = ((H + R)./R)^2;

    Slant_Range = R.*(sqrt(V-X) - sind(E));  % Slant_Range(E): distance from user to satellite
    %Sorting distance 780 to 2325 km
    Sort_SR=sort( Slant_Range);  
 
    %% Points (uniform step sizes) for base of coverage triangle 
    
    %b = sqrt(c^2 - a^2)  Right angle triangle formula

    % dMAX is the maximum propogation link (user to satellite)

    FootPrint_R = sqrt(max( Slant_Range)^2 - H^2);
    
    Ground_distance = 1:90e3:FootPrint_R; %Step size = 90e3 km

    Distance = zeros(1,length(Ground_distance));
    
    for loop=1:length(Ground_distance)
    
        Distance(1,loop) = sqrt(H^2 + Ground_distance(loop)^2);
    
    end
    
    %% Finding Elevation angle as function of maximum distance
    %https://www.researchgate.net/publication/220099290 ...
    ..._The_Range_and_Horizon_Plane_Simulation_for_Ground_Stations_of_Low_Earth_Orbiting_LEO_Satellites
    
    %Equation 8
    
    Elevation_Angle_Ground_Distance=zeros(1,length(Distance));
    
    for p=1:length(Distance)
    
        
        Elevation_Angle_Ground_Distance(1,p) = (H*(H+2*R) - Distance(p)^2)/(2*Distance(p)*R);
    
    end

    Elevation_Angles=asind(Elevation_Angle_Ground_Distance);

end