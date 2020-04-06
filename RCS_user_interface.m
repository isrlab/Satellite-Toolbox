clc; clear; close all;

idx = input('Press 1 to run the simulation with default inputs, 2 to enter new inputs: ');

switch idx
    case 1
        % default inputs
        userInput.tle = 'TLE_Intelsat36.txt';
        userInput.m = 3253;         % satellite mass 3253 kg
        userInput.a = 5.2;              % one dimension 5.2 m
        userInput.b = 3.4;              % another dimension 3.4 m
        userInput.t0 = 0;               % starting at t = 0
        userInput.n = 750;              % 750 orbits, spanning over 2 years
        userInput.dt = 60;              % data are available every 1 minute
        
        % Goldstone Deep Space Communications Complex
        userInput.lambdaD = -116.7846;     % radar station longitude 
        userInput.phiD = 35.2824;            % radar station latitude 


        userInput.lambdaR = 50e-2;  % radar wavelength 50 cm
        
        
    case 2
        disp('=======================================================');
        disp('Please provide the following data towards computation of Radar Cross-Section (RCS).');
        disp('=======================================================');
        disp(' ');                                                                                                                                
        disp('It is assumed that the satellite TLE data is saved in a text file.');
        disp(' ');
        userInput.tle = input('Name of the satellite TLE file with extension ".txt": ', 's');
        disp('-----------------------------------------------------------------------------------------------');
        disp(' ');   
        userInput.m = input('Mass of the satellite in kg: ');
        disp('-----------------------------------------------------------------------------------------------');
        disp(' ');   
        disp('It is assumed that the satellite is a rectangular flat plate with dimensions 2a and 2b.');
        disp('Please provide the plate dimensions as measured, without multiplication or division by 2.'); 
        disp(' ');
        userInput.a = input('2a in meters: ');
        userInput.b = input('2b in meters: ');
        disp('-----------------------------------------------------------------------------------------------');
        disp(' ');  
        disp('Please provide the initial time in seconds.');
        disp('Suggested initial time = 0. This corresponds to the epoch time in the TLE.');
        disp(' ');
        userInput.t0 = input('Initial time in seconds: ');
        disp('---------------------------------------------------------------------');
        disp(' ');  
        disp('Please provide the number of satellite orbits.');
        disp('For a geostationary satellite, the suggested number of orbits = 750. This spans a little over two years.');
        disp(' ');
        userInput.n = input('Number of orbits: ');
        disp('---------------------------------------------------------------------');
        disp(' '); 
        disp('Time-step is the time interval at which the computed data are available.');
        disp('Please provide the time-step in seconds. A good choice is 60 seconds, i.e. 1 minute.');
        disp(' ');
        userInput.dt = input('Time-step in seconds: ');
        disp('---------------------------------------------------------------------');
        disp(' '); 
        disp('Please provide the longitude and latitude of the radar station in degrees.');
        disp('--------------------------------------------------------------------');
        disp(' '); 
        disp('Longitude positive means eastern hemisphere, negative means western hemisphere.');
        disp(' '); 
        userInput.lambdaD = input('Longitude of the radar station in degrees (between -180 and 180): ');
        disp('----------------------------------------------------------------------------');
        disp(' '); 
        disp('Latitude positive means northern hemisphere, negative means southern hemisphere.');
        disp(' '); 
        userInput.phiD = input('Latitude of the radar station in degrees (between -90 and 90): ');
        disp('----------------------------------------------------------------------------------------------');
        disp(' ');
        userinput.lambdaR = input('Radar wavelength in meters: ');
        disp('========================================================');
        
        
end

disp(' ');
disp('All set to run the main subroutine.');
disp('=======================================================');

RCS_main_subroutine(userInput);

disp(' ');
disp('=======================================================');
disp('Figures are ready.');
disp('The last two figures plot the RCS data. The other figures show the intermediate analysis.');
disp('=======================================================');
