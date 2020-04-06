function [rcsdb_v, rcsdb_h] = rcs_rect_plate_modified(a, b, lambda, theta_deg)

%----------------------------------------------------------------------------
% Credit: 
% This function is modified from "rcs_rect_plate" in Chapter 2, Listing
% 2.8, pages 40 - 41 of Chapman and Hall (CRC, 2000).
%----------------------------------------------------------------------------

% This function computes the backscattered RCS for a rectangular flat plate.
% The RCS is computed for vertical and horizontal polarization based on
% Eq.s(2.50)through (2.60). For aspect angle close to 90 deg, 
% horizontal polarization is zero, and vertical polarization is computed
% using equation (2.61).

eps = 0.000001;
% lambda = 3.0e+8 / freq;
ka = 2 * pi * a / lambda;

% aspect angle converted to radians
theta = (pi/180) * theta_deg;

% Use (2.50) - (2.60) for theta in [0, 85 deg]

if theta_deg <= 85
    
        % compute the individual components of vertical polarization
        sigma1v = cos(ka *sin(theta)) - 1i * sin(ka *sin(theta)) / sin(theta);
        sigma2v = exp(1i * ka - (pi /4)) / (sqrt(2 * pi) *(ka)^1.5);
        sigma3v = (1 + sin(theta)) * exp(-1i * ka * sin(theta)) / ...
        (1 - sin(theta))^2;
        sigma4v = (1 - sin(theta)) * exp(1i * ka * sin(theta))/ ...
        (1 + sin(theta))^2;
        sigma5v = 1. - (exp(1i * 2 * ka - (pi / 2)) / (8 * pi * (ka)^3));

        % compute the individual components of vertical polarization 
        sigma1h = cos(ka *sin(theta)) + 1i * sin(ka *sin(theta)) / sin(theta);
        sigma2h = 4 * exp(1i * ka * (pi / 4)) / (sqrt(2 * pi * ka));
        sigma3h = exp(-1i * ka * sin(theta)) / (1 - sin(theta));
        sigma4h = exp(1i * ka * sin(theta)) / (1 + sin(theta));
        sigma5h = 1 - (exp(1j * 2 * ka + (pi / 4)) / 2 * pi * ka);

        % compute vertical polarization RCS
        rcs_v = (b^2 / pi) * (abs(sigma1v - sigma2v *((1/ cos(theta)) ...
        + .25 * sigma2v* (sigma3v + sigma4v)) * (sigma5v)^-1))^2 + eps;

        % compute horizontal polarization RCS
        rcs_h = (b^2 / pi) * (abs(sigma1h - sigma2h *((1/ cos(theta)) ...
        - .25* sigma2h* (sigma3h + sigma4h)) * (sigma5h)^-1))^2 + eps;    


else
%     % Use (2.61) for vertical polarization, horizontal polarization = 0
%     rcs_v = a*b^2/lambda * ((1 + (pi/2)/(2*a/lambda)^2)) + (1 - (pi/2)/(2*a/lambda)^2)*cos(2*ka - 3*pi/5);
%     rcs_h = 0;

        rcs_v = NaN;
        rcs_h = NaN;
end

% convert to dBsm
rcsdb_v = 10*log10(rcs_v);
rcsdb_h = 10*log10(rcs_h);

