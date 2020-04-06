function [] = RCS_main_subroutine(userInput)

tle = userInput.tle;   % Two-Line Element filename as a string

count = 0; % initializing figure count, will be incremented every time a new figure is generated
lw = 2.2;   % line width in figures 
fs = 22;     % font size in figures

disp('-----------------------------------------------');
disp('Two-Line Element of the Satellite ... ');
disp('-----------------------------------------------');
fileId=fopen(tle,'r');    % GEO: Geosynchronous Earth Orbit

lineData=fgetl(fileId); % read first line which contains the name of the satellite
disp(lineData);

global Asc m;
Asc = userInput.a * userInput.b;
m = userInput.m;
a_dim = 0.5* userInput.a;
b_dim = 0.5* userInput.b;

% Declare a few other global variables needed to compute J2 perturbation 
% and solar radiation pressure during orbit propagation
global J2 R S_B c 
S_B = 1367;            % Energy flux of radiation photons on Earth's orbit, in unit of W/m^2.
% S_B = 0;             % A way to artificially turn off solar radiation pressure

c = 299792458;        % Speed of light in m/s
J2 = 0.00108263;     % the second zonal harmonic
% J2 = 0;                 % A way to artificially turn off J2 perturbation

R = 6371071.03;       % earth's equatorial radius in meters

lineData=fgetl(fileId); % read the second line which contains epoch information
disp(lineData);

line2 = split(lineData);
epoch = str2double(line2(4));

formatOut = 'yyyy mm dd HH MM SS'; 
start = str2num(datestr(datetime(epoch, 'convertfrom','posixtime'), formatOut)); 

global year month day UT
year = start(1);
month = start(2);
day = start(3);
UT = start(4) + start(5)/60 + start(6)/3600; % universal time in hour


lineData=fgetl(fileId); % read the third line which contains the satellite initial condition information
disp(lineData);
fclose(fileId);

data = str2num(lineData);

% Using degree-to-radian and radian-to-degree conversions where necessary
global r2d d2r;
d2r = pi/180;
r2d = 180/pi;

% initial condition for inclination
i0 = d2r*data(1,3);
% i0 = 0;       % A way to artificially restrict the satellite to earth's equatorial plane

% initial condition for right ascension
Omega0 = d2r*data(1,4);

% initial condition for eccentricity
e0 = data(1,5)*1e-7;  % the value shown in the TLE is eccentricity scaled up by 10^7
% e0 = 0;       % A way to artificially make the orbit circular

% initial condition for argument of perigee
omega0 = d2r*data(1,6);

% initial condition for mean anomaly
M0 = d2r*data(1,7);

% initial condition for mean motion
n0 = data(1,8);         % n0 is in the unit of rev/day
% n0 = 1;                 % A way to artificially make the satellite orbit the earth in exactly 24 hours

% convert initial condition for mean motion to that for semi-major axis
% First convert rev/day to rad/s

% 24 hrs per day, 1 hr = 60 mins, 1 min = 60 secs
seconds_per_day = 24*60*60;
n0_scaled = 2*pi/seconds_per_day * n0;

% Use gravitational constant and mass of earth
G = 6.67408e-11;
M = 5.792e24; % the mass of earth in kg

global mu;
mu = G*M;
a0bar = mu^(1/3) / (n0_scaled)^(2/3);


% Make sure that the initial conditions are correct.

disp('--------------------------------------------------------------------');
disp('Initial conditions in the Keplerian elements: ');
disp('--------------------------------------------------------------------');
disp(['Average semi-major axis = ', num2str(a0bar/1000), ' km']);
disp(['Eccentricity = ', num2str(e0)]);
disp(['Inclination = ', num2str(r2d*i0), ' deg']);
disp(['Right ascension = ', num2str(r2d*Omega0), ' deg']);
disp(['Argument of perigee = ', num2str(r2d*omega0), ' deg']);
disp(['Mean anomaly = ', num2str(r2d*M0), ' deg']);
disp('---------------------------------------------------------------------');


% --------------------------------------------------------------------
%       Convert initial conditions to equinoctial coordinates
%---------------------------------------------------------------------
f0 = e0*cos(Omega0 + omega0);
g0 = e0*sin(Omega0 + omega0);
h0 = tan(i0/2)*cos(Omega0);
k0 = tan(i0/2)*sin(Omega0);

% True anomaly is needed to compute the initial condition of true longitude, L0
theta0 = M0 + (2*e0 - 1/4*e0^3) * sin(M0) + 5/4*e0^2*sin(2*M0) + 13/12*e0^3*sin(3*M0); % formula obtained from Wikipedia
% L0 = Omega0 + omega0 + theta0;
L0 = meaningful_angle(Omega0 + omega0 + theta0);

disp('--------------------------------------------------------------------');
disp('             Initial conditions in the equinoctial elements: ');
disp('--------------------------------------------------------------------');
disp(['Average semi-major axis = ', num2str(a0bar/1000), ' km']);
disp(['f = ', num2str(f0)]);
disp(['g = ', num2str(g0)]);
disp(['h = ', num2str(h0)]);
disp(['k = ', num2str(k0)]);
disp(['True longitude = ', num2str(r2d*L0), ' deg']);
disp('---------------------------------------------------------------------');

% Convert the initial conditions back to Keplerian elements 
% to make sure there is no error in coding so far

e0_check = sqrt(f0^2 + g0^2);
% p0_check = a0bar*(1 - e0_check^2);
i0_check = 2*atan2(sqrt(h0^2 + k0^2),ones(size(h0)));
omega0_plus_Omega0_check = atan2(g0,f0);  
Omega0_check = atan2(k0,h0);
omega0_check = omega0_plus_Omega0_check - Omega0_check;
theta0_check = L0 - omega0_plus_Omega0_check;


disp('---------------------------------------------------------------------------------');
disp('  Recomputed initial conditions in the Keplerian elements for sanity check: ');
disp('---------------------------------------------------------------------------------');
disp('           Showing true anomaly, which should be close enough  ...');
disp('                   ... to the mean anomaly read from TLE: ');
disp('----------------------------------------------------------------------------------');
disp(['Average semi-major axis = ', num2str(a0bar/1000), ' km']);
disp(['Eccentricity = ', num2str(e0_check)]);
disp(['Inclination = ', num2str(r2d*meaningful_angle(i0_check)), ' deg']);
disp(['Right ascension = ', num2str(r2d*meaningful_angle(Omega0_check)), ' deg']);
disp(['Argument of perigee = ', num2str(r2d*meaningful_angle(omega0_check)), ' deg']);
disp(['True anomaly = ', num2str(r2d*meaningful_angle(theta0_check)), ' deg']);
disp('---------------------------------------------------------------------');


disp('-----------------------------------------------------------------------------------');
disp('            Running the simulation for satellite orbital dynamics: ...                              ');
disp('-----------------------------------------------------------------------------------');
disp('The following are included: ');
disp('1. Semi-major axis as read from TLE. No uncertainty.');
disp('2. J2 perturbation included');
disp('3. Solar radiation pressure included with coefficient = 1.5 ...');
disp('    which is halfway between black body and perfect reflector.');
disp('------------------------------------------------------------------------------------');

global initial_time
initial_time = userInput.t0;
number_of_orbits = userInput.n;
time_step = userInput.dt;

final_time = initial_time + (86400/n0 * number_of_orbits + time_step - rem(seconds_per_day/n0 * number_of_orbits, time_step));  % time needed to complete 5 orbits
t = initial_time : time_step : final_time;
x0 = [a0bar f0 g0 h0 k0 L0];
disp('-----------------------------------------------------------------------');
disp('   Propagating in equinoctial coordinates ...');
disp('-----------------------------------------------------------------------');
[t,x] = ode45(@propagation_satellite_J2_SRP_deterministic,t,x0);
disp('Propagation complete. Equinoctial data ready.');
disp('------------------------------------------------------------------------');

time = t./(seconds_per_day/n0);  % this is time in orbit periods

% equinoctial variables
a = x(:,1);
f = x(:,2);
g = x(:,3);
h = x(:,4);
k = x(:,5);
L = x(:,6);

disp('Converting data to Keplerian coordinates ...');
disp('--------------------------------------------------------------------------');
% Converting to Keplerian
ecc = sqrt(f.^2 + g.^2);
p = a.*(1 - ecc.^2);
inc = 2*atan2(sqrt(h.^2 + k.^2),ones(size(h)));
omega_plus_Omega = atan2(g,f);  
Omega = atan2(k,h);
omega = omega_plus_Omega - Omega;
theta = L - omega_plus_Omega;

disp('---------------------------------------------------------------------------------------------------');
disp('Keplerian data ready. Obtaining satellite positions in Earth-Centric Inertial (ECI) frame ...');
disp('---------------------------------------------------------------------------------------------------');

% positions
r = p./(1+ecc.*cos(theta));
x = (cos(Omega).*cos(omega+theta) - sin(Omega).*sin(omega+theta).*cos(inc)) .* r;
y = (sin(Omega).*cos(omega+theta) + cos(Omega).*sin(omega+theta).*cos(inc)) .* r;
z = (sin(omega+theta).*sin(inc)) .* r;

% Satellite position vector in ECI frame
rT_ECI = [x y z];

disp('----------------------------------------------------------------------------------------------');
disp('SATELLITE position vector components in the ECI frame are now available at all times.');
disp('----------------------------------------------------------------------------------------------');

disp('------------------------------------------------------------');
disp('Plotting the shadow function nu as function of time ... ');
disp('------------------------------------------------------------');

    
disp('-----------------------------------------------------------------------');
disp('Computing SUN position vector components in the ECI frame ... ');
disp('-----------------------------------------------------------------------');

rN_ECI = zeros(length(time),3);
rN_polar = zeros(length(time),3);

% Compute and plot shadow function 'nu' for solar radiation pressure
nu = zeros(length(time),1);
    
   for i = 1:length(time)

       current_time = UT + (initial_time + t(i))/3600;
       [rNmag, epsilon, lambda, nu_curr] = JulianDay_calculation(day,month,year,current_time,x(i),y(i),z(i),r(i),R);
       
       % store earth-sun distance, solar longitude and earth's obliquity
       rN_polar(i,1) = rNmag;
       rN_polar(i,2) = lambda;
       rN_polar(i,3) = epsilon;
       
       % epsilon and lambda are in degrees
       rN_ECI(i,1) = rNmag * cosd(lambda);
       rN_ECI(i,2) = rNmag * sind(lambda)*cosd(epsilon);
       rN_ECI(i,3) = rNmag * sind(lambda)*sind(epsilon);
       
       nu(i) = nu_curr;
   end
   
   figure(count+1); clf; count = count + 1;
    plot(time, nu, 'b','linewidth',lw);
    grid on;
    set(gca,'linewidth',lw,'fontsize',fs);
    set(gcf,'color',[1 1 1]);
    ylabel('shadow function \nu');
    xlabel('orbit period');
    xlim([0 min(180,number_of_orbits)]);
    ylim([-0.2 1.2]);     
   
  % plot the sun position to see how much the sun is moving w.r.t. ECI frame
  % sun's movement should be very small
  
  % plot the polar quantities describing sun position
  
  figure(count+1); clf; count = count+1;
  
  subplot(3,1,1);
  plot(time, rN_polar(:,1)/1e3, '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('earth-sun distance');
  ylabel('r_N (km)');
  xlim([0 number_of_orbits]);
  
  subplot(3,1,2);
  plot(time, rN_polar(:,2), '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('solar longitude');
  ylabel('\lambda_N (deg)');
  xlim([0 number_of_orbits]);
  
  subplot(3,1,3);
  plot(time, rN_polar(:,3), '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('earth obliquity');
  ylabel('\epsilon (deg)');
  xlabel('orbit period');
  xlim([0 number_of_orbits]);
  
  %plot cartesian quantities describing sun position
  figure(count+1); clf; count = count+1;
  
  subplot(3,1,1);
  plot(time, rN_ECI(:,1)/1e3, '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('x-position of sun in ECI frame');
  ylabel('x_N (km)');
  xlim([0 number_of_orbits]);
  
  subplot(3,1,2);
  plot(time, rN_ECI(:,2)/1e3, '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('y-position of sun in ECI frame');
  ylabel('y_N (km)');
  xlim([0 number_of_orbits]);
  
  subplot(3,1,3);
  plot(time, rN_ECI(:,3)/1e3, '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('z-position of sun in ECI frame');
  ylabel('z_N (km)');
  xlabel('orbit period');
  xlim([0 number_of_orbits]);
  
  % Now plot the deviations with time
  figure(count+1); clf; count = count+1;
  
  subplot(3,1,1);
  plot(time, (rN_polar(:,1)-rN_polar(1,1)*ones(size(rN_polar,1),1))/1e3, '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('variation in earth-sun distance from initial condition');
  ylabel('\delta r_N (km)');
  xlim([0 number_of_orbits]);
  
  subplot(3,1,2);
  plot(time, rN_polar(:,2)-rN_polar(1,2)*ones(size(rN_polar,1),1), '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('variation in solar longitude from initial condition');
  ylabel('\delta \lambda_N (deg)');
  xlim([0 number_of_orbits]);
  
  subplot(3,1,3);
  plot(time, rN_polar(:,3)-rN_polar(1,3)*ones(size(rN_polar,1),1), '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('variation in earth obliquity from initial condition');
  ylabel('\delta \epsilon (deg)');
  xlabel('orbit period');
  xlim([0 number_of_orbits]);
  
  %plot cartesian quantities describing sun position
  figure(count+1); clf; count = count+1;
  
  subplot(3,1,1);
  plot(time, (rN_ECI(:,1)-rN_ECI(1,1)*ones(size(rN_ECI,1),1))/1e3, '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('variation in sun x-position from initial condition');
  ylabel('\delta x_N (km)');
  xlim([0 number_of_orbits]);
  
  subplot(3,1,2);
  plot(time, (rN_ECI(:,2)-rN_ECI(1,2)*ones(size(rN_ECI,1),1))/1e3, '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('variation in sun y-position from initial condition');
  ylabel('\delta y_N (km)');
  xlim([0 number_of_orbits]);
  
  subplot(3,1,3);
  plot(time, (rN_ECI(:,3)-rN_ECI(1,3)*ones(size(rN_ECI,1),1))/1e3, '-b', 'linewidth', lw); grid on;
  set(gca,'linewidth', lw, 'fontsize', fs);
  set(gcf, 'color', [1 1 1]);
  title('variation in sun z-position from initial condition');
  ylabel('\delta z_N (km)');
  xlabel('orbit period');
  xlim([0 number_of_orbits]);
 

disp('----------------------------------------------------------------------------------------------');
disp('SUN position vector components in the ECI frame are now available at all times.');
disp('----------------------------------------------------------------------------------------------');

disp('--------------------------------------------------------------------------------------');
disp('Computing RADAR STATION position vector components in the ECI frame ... ');
disp('--------------------------------------------------------------------------------------');

rD_EFE_latlong = zeros(length(time),2);
rD_ECI_latlong = zeros(length(time),2);
rD_ECI = zeros(length(time),3);
rD_EFE = zeros(length(time),3);

rT_EFE_latlong = zeros(length(time),2);
rT_ECI_latlong = zeros(length(time),2);
rT_EFE = zeros(length(time),3);

lambdaD = d2r*userInput.lambdaD;
phiD = d2r*userInput.phiD;

% find initial lambda_E, which is how much the earth-fixed equatorial frame
% is off from the ECI frame at t = 0.

omegaE = 2*pi/seconds_per_day;       % spin angular speed of earth in rad/s
lambdaE_0 = 0 + omegaE * ((UT+initial_time/3600 - 12)*3600);    % essentially t = 0 corresponds to (UT+initial_time/3600) in hours, 
                                                                         % and lambdaE = 0 at 1200 hours
                                                                         
for i = 1:length(time)
    lambdaE = lambdaE_0 + omegaE*(initial_time + t(i));   % 't' has time in seconds, 'time' has time in orbit periods
    
    A = [cos(lambdaE) -sin(lambdaE)  0;
           sin(lambdaE)   cos(lambdaE)  0;
           0                              0           1];
       

    B = [cos(phiD)*cos(lambdaD);
            cos(phiD)*sin(lambdaD);
            sin(phiD)];
        
    rD_ECI(i,:) = (A * B)'*R;
    
    rT_EFE(i,:) = (A'*rT_ECI(i,:)')';
    
    rD_EFE(i,:) = B'*R;
    
    rD_ECI_latlong(i,:) = [atan2d(rD_ECI(i,2), rD_ECI(i,1))     asind(rD_ECI(i,3) / norm(rD_ECI(i,:)))];
    rD_EFE_latlong(i,:) = [atan2d(rD_EFE(i,2), rD_EFE(i,1))     asind(rD_EFE(i,3) / norm(rD_EFE(i,:)))];
    
    rT_ECI_latlong(i,:) = [atan2d(rT_ECI(i,2), rT_ECI(i,1))     asind(rT_ECI(i,3) / norm(rT_ECI(i,:)))];
    rT_EFE_latlong(i,:) = [atan2d(rT_EFE(i,2), rT_EFE(i,1))     asind(rT_EFE(i,3) / norm(rT_EFE(i,:)))];
    
    % Need to transpose because each ROW is saving the values at each time-instant.
end

disp('--------------------------------------------------------------------------------------------------------');
disp('RADAR STATION position vector components in the ECI frame are now available at all times ... ');
disp('--------------------------------------------------------------------------------------------------------');


figure(count+1); clf; count = count+1;
subplot(2,1,1);
plot(time, rT_ECI_latlong(:,1), '-b', time, rT_EFE_latlong(:,1), '-.b', ...
       time, rD_ECI_latlong(:,1), '-r', time, rD_EFE_latlong(:,1), '-.r', ...
       'linewidth', lw);
grid on;
title('Longitude (deg)');
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
xlim([0 min(5,number_of_orbits)]);  % 5 orbits are sufficient to observe the periodic behavior
ylim([-200 200]);

subplot(2,1,2);
plot(time, rT_ECI_latlong(:,2), '-b', time, rT_EFE_latlong(:,2), '-.b', ...
       time, rD_ECI_latlong(:,2), '-r', time, rD_EFE_latlong(:,2), '-.r', ...
       'linewidth', lw);
title('Latitude (deg)');
grid on;
legend('satellite, w.r.t. inertial earth frame', ... 
           'satellite, w.r.t. rotating earth frame', ...
           'radar, w.r.t. inertial earth frame', ...
           'radar, w.r.t. rotating earth frame');
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
xlabel('orbit period');
xlim([0 min(5,number_of_orbits)]);
ylim([-100 100]);

% Satellite longitude w.r.t. rotating earth varies periodically throughout
% the year, albeit of small magnitude. 
figure(count+1); clf; count = count+1;
subplot(2,1,1);
plot(time, rT_EFE_latlong(:,1), '-.b', ...
       time, rD_EFE_latlong(:,1), '-.r', ...
       'linewidth', lw);
grid on;
title('Longitude (deg)');
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
 xlim([0 number_of_orbits]);
 ylim([-200 200]);


subplot(2,1,2);
plot(time, rT_EFE_latlong(:,2), '-.b', ...
       time, rD_EFE_latlong(:,2), '-.r', ...
       'linewidth', lw);
title('Latitude (deg)');
grid on;
legend(...
           'satellite, w.r.t. rotating earth frame', ...
           'radar, w.r.t. rotating earth frame');
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
xlabel('orbit period');
xlim([0 number_of_orbits]);
ylim([-100 100]);


disp('-------------------------------------------------------------------------------------------');
disp('Computing desired radar azimuth and elevation to track the satellite at all times ... ');
disp('-------------------------------------------------------------------------------------------');

rDT_ECI = rT_ECI - rD_ECI;
rDT_EFE = rT_EFE - rD_EFE;

rDT_LD = zeros(length(time),3);
azimuth = zeros(length(time),1);
elevation = zeros(length(time),1);
elevation_acceptable = zeros(length(time),1);

for i = 1:length(time)
    
    lambdaE = lambdaE_0 + omegaE*(initial_time + t(i));
    % convert rDT_ECI to local frame at D
    A = [-sin(lambdaD)                      cos(lambdaD)                      0;
            -sin(phiD)*cos(lambdaD)    -sin(phiD)*sin(lambdaD)    cos(phiD);
            cos(phiD)*cos(lambdaD)     cos(phiD)*sin(lambdaD)    sin(phiD)];
        
    B = [cos(lambdaE)    sin(lambdaE)     0;
            -sin(lambdaE)   cos(lambdaE)    0;
            0                        0                    1];
        
     
    rDT_LD(i,:) = (A*B*rDT_ECI(i,:)')';
    
    xDT = rDT_LD(i,1);  yDT = rDT_LD(i,2);  zDT = rDT_LD(i,3);
    
    azimuth(i) = atan2d(xDT,yDT);
    elevation(i) = asind(zDT/norm(rDT_ECI(i,:)));
    
    if elevation(i) >= 0
        elevation_acceptable(i) = elevation(i);
    else 
        elevation_acceptable(i) = NaN;
    end
    
end

figure(count+1); clf; count = count+1;
[xE,yE,zE] = sphere;
% Scale to desire radius.
radius = R/1e3;
xE = xE * radius;
yE = yE * radius;
zE = zE * radius;
% Translate sphere to new location.
offset = 0;
% Plot as surface.
surf(xE+offset,yE+offset,zE+offset);
set(gca,'linewidth',lw,'fontsize',fs);
set(gcf,'color',[1 1 1]);
xlabel('x (km)');
ylabel('y (km)');
zlabel('z (km)');
axis equal;
hold on;
plot3(x/1e3, y/1e3, z/1e3, '-b', ...
        rD_ECI(:,1)/1e3, rD_ECI(:,2)/1e3, rD_ECI(:,3)/1e3, '-r','linewidth', lw); grid on;
legend('earth', 'satellite trajectory', 'rotating sensor position');

% Plot relative x, y, z positions of satellite w.r.t. radar
figure(count+1); clf; count = count + 1;
subplot(3,1,1);
plot(time, rDT_ECI(:,1)/1e3, '-.r', ...
       time, rDT_EFE(:,1)/1e3, '-.b', ...
       time, rDT_LD(:,1)/1e3, '-.k', ...
       'linewidth', lw); grid on;
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
title('relative x-position of satellite w.r.t. radar');
ylabel('x_{DT} (km)');
legend('inertial earth frame', 'rotating earth frame', 'radar local frame');
xlim([0 min(5,number_of_orbits)]);

subplot(3,1,2);
plot(time, rDT_ECI(:,2)/1e3, '-.r', ...
       time, rDT_EFE(:,2)/1e3, '-.b', ...
       time, rDT_LD(:,2)/1e3, '-.k', ...
       'linewidth', lw); grid on; 
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
title('relative y-position of satellite w.r.t. radar');
ylabel('y_{DT} (km)');
xlim([0 min(5,number_of_orbits)]);

subplot(3,1,3);
plot(time, rDT_ECI(:,3)/1e3, '-.r', ...
       time, rDT_EFE(:,3)/1e3, '-.b', ...
       time, rDT_LD(:,3)/1e3, '-.k', ...
       'linewidth', lw); grid on; 
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
title('relative z-position of satellite w.r.t. radar');
ylabel('z_{DT} (km)');
xlim([0 min(5,number_of_orbits)]);
xlabel('orbit period');


% Plot only the slow variations
% Plot relative x, y, z positions of satellite w.r.t. radar
figure(count+1); clf; count = count + 1;
subplot(3,1,1);
plot(time, rDT_EFE(:,1)/1e3, '-.b', ...
       time, rDT_LD(:,1)/1e3, '-.k', ...
       'linewidth', lw); grid on;
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
title('relative x-position of satellite w.r.t. radar');
ylabel('x_{DT} (km)');
legend('rotating earth frame', 'radar local frame');
xlim([0 number_of_orbits]);

subplot(3,1,2);
plot(time, rDT_EFE(:,2)/1e3, '-.b', ...
       time, rDT_LD(:,2)/1e3, '-.k', ...
       'linewidth', lw); grid on; 
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
title('relative y-position of satellite w.r.t. radar');
ylabel('y_{DT} (km)');
xlim([0 number_of_orbits]);

subplot(3,1,3);
plot(time, rDT_EFE(:,3)/1e3, '-.b', ...
       time, rDT_LD(:,3)/1e3, '-.k', ...
       'linewidth', lw); grid on; 
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
title('relative z-position of satellite w.r.t. radar');
ylabel('z_{DT} (km)');
xlim([0 number_of_orbits]);
xlabel('orbit period');

% Now plot azimuth and elevation
figure(count+1); clf; count = count + 1;
subplot(2,1,1);
plot(time, azimuth, '-.r', ... 
      time, -180*ones(length(azimuth),1), '--k', ...
      time, 180*ones(length(azimuth),1), '-.k','linewidth', lw); grid on;
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
title('azimuth angle');
ylabel('\lambda_{DT} (deg)');
xlim([0 number_of_orbits]);
ylim([-200 200]);

subplot(2,1,2);
plot(time, elevation, '-.b', time, elevation_acceptable, '-.r', ...
      time, 0*ones(length(elevation),1), '--k', ...
      time, 90*ones(length(elevation),1), '-.k', ...
      'linewidth', lw); grid on; 
legend('computed', 'acceptable', 'lower bound', 'upper bound'); 
set(gca, 'linewidth', lw, 'fontsize', fs);
set(gcf, 'color', [1 1 1]);
title('elevation angle');
ylabel('\phi_{DT} (deg)');
xlabel('orbit period');
xlim([0 number_of_orbits]);
ylim([-100 100]);


if all(isnan(elevation_acceptable))
    disp('ERROR: The satellite is ALWAYS outside the reach of the radar.');
    disp('Please change the latitude and longitude of the radar station.');
    disp('---------------------END OF CODE -------------------------------------');
    
else
    disp('-------------------------------------------------------------------------------------------');
    disp('Computing vertical and horizontal backscattered RCS at all times ... ');
    disp('-------------------------------------------------------------------------------------------');    

    rTN_ECI = rN_ECI - rT_ECI;

    alpha = zeros(length(time),1);

    lambda_R = userInput.lambdaR;

    rcsdb_v = zeros(length(time),1);
    rcsdb_h = zeros(length(time),1);

    for i = 1:length(time)
        if elevation(i) < 0  % discard computation of rcs, because radar will not be able to track the satellite
            alpha(i) = NaN;
            rcsdb_v(i) = NaN;
            rcsdb_h(i) = NaN;
        else
            % compute rcs in decibels per square meter
            alpha(i) = acosd((rDT_ECI(i,:)*rTN_ECI(i,:)')/(norm(rDT_ECI(i,:))*norm(rTN_ECI(i,:))));
            if alpha(i) > 90
                alpha(i) = 180 - alpha(i);  % symmetric about 90 deg.
            end
            [rcsdb_v(i), rcsdb_h(i)] = rcs_rect_plate_modified(a_dim, b_dim, lambda_R, alpha(i));
        end
    end

    figure(count+1); clf; count = count + 1;
    subplot(3,1,1);
    plot(time, alpha, '-b', 'linewidth', lw); grid on;
    set(gca, 'linewidth', lw, 'fontsize', fs);
    set(gcf, 'color', [1 1 1]);
    title('aspect angle (deg)');
    ylabel('\alpha');
    xlim([0 number_of_orbits]);
    ylim([-10 100]);

    subplot(3,1,2);
    plot(time, rcsdb_v, '-b', 'linewidth', lw); grid on;
    set(gca, 'linewidth', lw, 'fontsize', fs);
    set(gcf, 'color', [1 1 1]);
    ylabel('RCS_V');
    title('vertical polarization (dBsm)');
    xlim([0 number_of_orbits]);

    subplot(3,1,3);
    plot(time, rcsdb_h, '-b', 'linewidth', lw); grid on;
    set(gca, 'linewidth', lw, 'fontsize', fs);
    set(gcf, 'color', [1 1 1]);
    ylabel('RCS_H');
    xlabel('orbit period');
    title('horizontal polarization (dBsm)');
    xlim([0 number_of_orbits]);


    % Zoom in to see RCS in a few orbits

    % find an orbit during which the radar will definitely track the satellite
    [~,idx] = max(elevation_acceptable);
    orbit_index = floor(time(idx));

    figure(count+1); clf; 
    subplot(3,1,1);
    plot(time, alpha, '-b', 'linewidth', lw); grid on;
    set(gca, 'linewidth', lw, 'fontsize', fs);
    set(gcf, 'color', [1 1 1]);
    title('aspect angle (deg)');
    ylabel('\alpha');
    xlim([orbit_index orbit_index+2]);
    ylim([-10 100]);

    subplot(3,1,2);
    plot(time, rcsdb_v, '-b', 'linewidth', lw); grid on;
    set(gca, 'linewidth', lw, 'fontsize', fs);
    set(gcf, 'color', [1 1 1]);
    ylabel('RCS_V');
    title('vertical polarization (dBsm)');
    xlim([orbit_index orbit_index+2]);

    subplot(3,1,3);
    plot(time, rcsdb_h, '-b', 'linewidth', lw); grid on;
    set(gca, 'linewidth', lw, 'fontsize', fs);
    set(gcf, 'color', [1 1 1]);
    ylabel('RCS_H');
    xlabel('orbit period');
    title('horizontal polarization (dBsm)');
    xlim([orbit_index orbit_index+2]);
        
end