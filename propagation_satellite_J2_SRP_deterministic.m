function out = propagation_satellite_J2_SRP_deterministic(t,in)

global mu;
a = in(1);
f = in(2);
g = in(3);
h = in(4);
k = in(5);
L = in(6);

 ap = get_perturbing_acceleration_J2_SRP_deterministic(t,in);
 S = ap(1);
 T = ap(2);
 N = ap(3);

s = sqrt(1-f^2-g^2);
W = 1 + f*cos(L) + g*sin(L);
A = f + cos(L*(1+W));
B = g + sin(L*(1+W));
X = 1 + h^2 + k^2;
e = sqrt(f^2 + g^2);
p = a*(1-e^2);

dadt = 2*sqrt(a/mu)*a/s*((f*sin(L) - g*cos(L))*S + W*T);
dfdt = sqrt(p/mu)*1/W* (W*sin(L)*S + A*T - g*(h*sin(L) - k*cos(L))*N);
dgdt = sqrt(p/mu)*1/W* (-W*cos(L)*S + B*T + f*(h*sin(L) + k*cos(L))*N);
dhdt = 1/2*sqrt(p/mu) * X/W * cos(L) * N;
dkdt = 1/2*sqrt(p/mu) * X/W * sin(L) * N;
dLdt = sqrt(mu/p^3)*W^2 + sqrt(p/mu) * 1/W *(h*sin(L) - k*cos(L))*N;

out = [dadt; dfdt; dgdt; dhdt; dkdt; dLdt];

