function [x1,y1,z1] = RK4GLO(te, ti, x0, y0, z0, Vx0, Vy0, Vz0, AxLS, AyLS, AzLS)
% RK4GLO calculates position at updated epoch ti of GLONASS satellite given initial
% positions and velocities of satellite at a reference epoch te using 4th order Runge-Kutta
% or RK4 numerical integration method of the differential equations
% modelling the satellite's motion.
% x0, y0, z0 - positions of satellite at reference epoch in PZ90 reference
% system
% Vx0, Vy0, Vz0 - velocities of SV at reference epoch in PZ90 reference
% system
% AxLS, AyLS, AzLS -  lunar  +  solar accelerations at reference epoch in PZ90
% refrence system
% ti - current epoch in seconds
% te - reference epoch in seconds

% Notes:
%   Luni-solar accelerations considered constant during integration period
h = ti-te;
% per GLONASS ICD,  integration period h must be less than 15 min
if abs(h)>900
    warning('Time difference between navigation message reference epoch and observation epoch exceeds recommended 15min');
end

%% Recalculate ephemeris parameters

% m, n, o - RK4 parameters for x, y, z respectively
% J, K, L - RK4 parameters for Vx, Vy, Vz respectively

% initial values
m1 = Vx0; 
n1 = Vy0; 
o1 = Vz0;
J1 = getAx(x0, y0, z0, AxLS);
K1 = getAy(x0, y0, z0, AyLS);
L1 = getAz(x0, y0, z0, AzLS);

m2 = Vx0 + 0.5*h*J1;
n2 = Vy0 + 0.5*h*K1;
o2 = Vz0 + 0.5*h*L1;

xstep = x0 + 0.5*h*m1;
ystep = y0 + 0.5*h*n1;
zstep = z0 + 0.5*h*o1;
J2 = getAx(xstep, ystep, zstep, AxLS);
K2 = getAy(xstep, ystep, zstep, AyLS);
L2 = getAz(xstep, ystep, zstep, AzLS);

m3 = Vx0 + 0.5*h*J2;
n3 = Vy0 + 0.5*h*K2;
o3 = Vz0 + 0.5*h*L2;

xstep = x0 + 0.5*h*m2;
ystep = y0 + 0.5*h*n2;
zstep = z0 + 0.5*h*o2;
J3 = getAx(xstep, ystep, zstep, AxLS);
K3 = getAy(xstep, ystep, zstep, AyLS);
L3 = getAz(xstep, ystep, zstep, AzLS);

m4 = Vx0 + h*J3;
n4 = Vy0 + h*K3;
o4 = Vz0 + h*L3;

% 4th order velocity parameters not needed for position
% xstep = x0 + h*m3;
% ystep = y0 + h*n3;
% zstep = z0 + h*o3;
% J4 = getAx(xstep, ystep, zstep, AxLS);
% K4 = getAy(xstep, ystep, zstep, AyLS);
% L4 = getAz(xstep, ystep, zstep, AzLS);

x1 = x0 + (h/6)*(m1 + 2*m2 + 2*m3 + m4);
y1 = y0 + (h/6)*(n1 + 2*n2 + 2*n3 + n4);
z1 = z0 + (h/6)*(o1 + 2*o2 + 2*o3 + o4);

    function Ax = getAx(x, y, z, AxLS)
        
        r = sqrt(x^2 + y^2 + z^2);
        rho = Const.AE/r;
        xbar = x/r;
        zbar = z/r;
        mubar = Const.MU/(r^2);
        Ax = -mubar*xbar + 1.5*Const.C20*mubar*xbar*(rho^2)*(1-5*zbar^2) + AxLS;
    end

    function Ay = getAy(x, y, z, AyLS)
        
        r = sqrt(x^2 + y^2 + z^2);
        rho = Const.AE/r;
        ybar = y/r;
        zbar = z/r;
        mubar = Const.MU/(r^2);
        Ay = -mubar*ybar + 1.5*Const.C20*mubar*ybar*(rho^2)*(1-5*zbar^2) + AyLS;
    end

    function Az = getAz(x, y, z, AzLS)
        
        r = sqrt(x^2 + y^2 + z^2);
        rho = Const.AE/r;
        zbar = z/r;
        mubar = Const.MU/(r^2);
        Az = -mubar*zbar + 1.5*Const.C20*mubar*zbar*(rho^2)*(3-5*zbar^2) + AzLS;
    end
end
