function out = get_perturbing_acceleration_J2_SRP_deterministic(t,x)
global mu;
global R J2 Asc S_B m c;
global day month year UT initial_time
% global t_data nu_data;

a = x(1);
f = x(2);
g = x(3);
h = x(4);
k = x(5);
L = x(6);
CR = 1.5;

% convert to Keplerian elements
e = sqrt(f^2 + g^2);
p = a*(1 - e^2);

inc = 2*atan2(sqrt(h^2 + k^2),ones(size(h)));
omega_plus_Omega = atan2(g,f);  
Omega = atan2(k,h);

omega = omega_plus_Omega - Omega;
theta = L - omega_plus_Omega;

% convert to ECI cartesian coordinates
r = p/(1+e*cos(theta));
x = (cos(Omega)*cos(omega+theta) - sin(Omega)*sin(omega+theta)*cos(inc)) * r;
y = (sin(Omega)*cos(omega+theta) + cos(Omega)*sin(omega+theta)*cos(inc)) * r;
z = (sin(omega+theta)*sin(inc)) * r;

% scaling factor which will multiply all components of the J2 perturbation
scaling_factor = 3*J2*mu*R^2/(2*r^4);

% components of the J2 perturbation in ECI frame
px_ECI = scaling_factor * x/r * (5*z^2/r^2 - 1);
py_ECI = scaling_factor * y/r * (5*z^2/r^2 - 1);
pz_ECI = scaling_factor * z/r * (5*z^2/r^2 - 3);

% convert the components to satellite Local Vertical Local Horizontal frame
% which essentially captures the satellite orbital plane and its normal

C3_omega_plus_theta = [cos(omega+theta) sin(omega+theta) 0;
                                       -sin(omega+theta)  cos(omega+theta) 0;
                                       0      0      1];
                                   
C1_i = [1    0   0;
            0  cos(inc)  sin(inc);
            0  -sin(inc) cos(inc)];
       
C3_Omega = [cos(Omega) sin(Omega) 0;
                     -sin(Omega)  cos(Omega) 0;
                      0      0      1];
                  
ap_J2 = C3_omega_plus_theta * C1_i * C3_Omega * [px_ECI; py_ECI; pz_ECI];

% Need to find Julian Day for the current time, and eventually Solar Radiation Pressure (SRP) at the current time
current_time = UT + (initial_time + t)/3600; % Note that UT is in hour, so converting current time to hour.
[~,epsilon, lambda, nu] = JulianDay_calculation(day,month,year,current_time,x,y,z,r,R); 

ap_SRP = - C3_omega_plus_theta * C1_i * C3_Omega * [cosd(lambda); sind(lambda)*cosd(epsilon); sind(lambda)*sind(epsilon)] * nu*S_B/c*CR*Asc/m; 
    
out =  ap_J2 + ap_SRP;