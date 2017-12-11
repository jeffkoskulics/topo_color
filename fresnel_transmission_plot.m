%Fresnel transmission

water_index = 1.333;
plastic_index = 1.575;
air_index = 1;

degrees = -48.6:0.01:48.6;
theta_i = degrees * pi / 180;

incident_index = water_index;
refracted_index = air_index;

Rs = ((incident_index*cos(theta_i) - refracted_index*sqrt(1 - ((incident_index/refracted_index)*sin(theta_i)).^2))./...
    ((incident_index*cos(theta_i) + refracted_index*sqrt(1 - ((incident_index/refracted_index)*sin(theta_i)).^2)))).^2;

Rp = ((incident_index*sqrt(1-((incident_index/refracted_index)*sin(theta_i)).^2) - refracted_index*cos(theta_i))./...
    (incident_index*sqrt(1-((incident_index/refracted_index)*sin(theta_i)).^2) + refracted_index*cos(theta_i))).^2;

R = (Rs + Rp)/2;

%Transmittance = 1 - Reflectance
T_water_air = 1 - R;

figure;
plot(degrees,T,'color','blue')
hold on

transmission_water_to_air = T_water_air(round(length(degrees)/2))

%Water-to-plastic transmission

incident_index = water_index;
refracted_index = plastic_index;

Rs = ((incident_index*cos(theta_i) - refracted_index*sqrt(1 - ((incident_index/refracted_index)*sin(theta_i)).^2))./...
    ((incident_index*cos(theta_i) + refracted_index*sqrt(1 - ((incident_index/refracted_index)*sin(theta_i)).^2)))).^2;

Rp = ((incident_index*sqrt(1-((incident_index/refracted_index)*sin(theta_i)).^2) - refracted_index*cos(theta_i))./...
    (incident_index*sqrt(1-((incident_index/refracted_index)*sin(theta_i)).^2) + refracted_index*cos(theta_i))).^2;

R = (Rs + Rp)/2;

T_water_plastic = 1 - R;

%Transmission from plastic to air

%new refracted angle
theta_i = asin(water_index * sin(theta_i) / plastic_index);

incident_index = plastic_index;
refracted_index = air_index;

Rs = ((incident_index*cos(theta_i) - refracted_index*sqrt(1 - ((incident_index/refracted_index)*sin(theta_i)).^2))./...
    ((incident_index*cos(theta_i) + refracted_index*sqrt(1 - ((incident_index/refracted_index)*sin(theta_i)).^2)))).^2;

Rp = ((incident_index*sqrt(1-((incident_index/refracted_index)*sin(theta_i)).^2) - refracted_index*cos(theta_i))./...
    (incident_index*sqrt(1-((incident_index/refracted_index)*sin(theta_i)).^2) + refracted_index*cos(theta_i))).^2;

R = (Rs + Rp)/2;

T_plastic_air = 1 - R;

T_water_plastic_air = T_water_plastic.*T_plastic_air;

plot(degrees,T_water_plastic_air,'color','black')

transmission_water_to_air_to_plastic = T_water_plastic_air(round(length(degrees)/2))