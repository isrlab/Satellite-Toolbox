function [r_es, epsilon, lambda, nu] = JulianDay_calculation(D, M, Y, current_time, xpos, ypos, zpos, r, R)

% Julian Day Number
% JDN = round((1461*(Y + 4800 + (M-14)/12)/4 + (367 * (M - 2 - 12*((M - 14)/12)))/12 - (3 * ((Y + 4900 + (M - 14)/12)/100))/4 + D - 32075; 
 JDN = round(D - 32075 + 1461*( Y + 4800 + ( M - 14 ) / 12 ) / 4 + 367*( M - 2 - ( M - 14 ) / 12 * 12 ) / 12 - 3*( ( Y + 4900 + ( M - 14 ) / 12 ) / 100 ) / 4);

% Add the Universal Time (UT) for the fractional part
% UT is in hours
JD = JDN + (current_time-12)/24;

% Number of days since J2K
n = JD - 2451545.0;

% angle between geocentric equatorial and ecliptic planes
% also known as earth's obliquity
% epsilon = r2d*meaningful_angle(d2r*(23.439 - 3.56*10^-7*n));
epsilon = mod(23.439 - 3.56*10^-7*n, 360);

% solar longitude 
% L = r2d*meaningful_angle(d2r*(280.459 + 0.98564736*n));
L = mod(280.459 + 0.98564736*n, 360);
% M = r2d*meaningful_angle(d2r*(357.529 + 0.98560023*n));
M = mod(357.529 + 0.98560023*n, 360);
% lambda = r2d*meaningful_angle(d2r*(L + 1.915*sind(M) + 0.0020*sind(2*M)));
lambda = mod(L + 1.915*sind(M) + 0.0020*sind(2*M), 360);
% distance between earth and sun
% r_es = (1.00014 - 0.01671*cosd(M) - 0.000140*cosd(2*M))*1.496e+11;
r_es = (1.00014 - 0.01671*cosd(M) - 0.000140*cosd(2*M))*149597870.691e+3;


% unit vector from earth to sun
earth_sun_unit_vec = [cosd(lambda); sind(lambda)*cosd(epsilon); sind(lambda)*sind(epsilon)];

% unit vector from earth to satellite
earth_satellite_unit_vec = [xpos/r; ypos/r; zpos/r];

% calculate three angles
theta0 = acos(earth_sun_unit_vec' * earth_satellite_unit_vec);
theta1 = acos(R/r);
theta2 = acos(R/r_es);

if theta1 + theta2 <= theta0
    nu = 0;
else
    nu = 1;
end