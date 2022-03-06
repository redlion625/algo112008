function x = simpleRK4(h, x, als)
% RK4GLO calculates position at updated epoch ti of GLONASS satellite given initial
% positions and velocities of satellite at a reference epoch te using 4th order Runge-Kutta
% or RK4 numerical integration method of the differential equations
% modelling the satellite's motion.
% h = tsv - tb
% x = [x0;y0,z0,vx0,vy0,vz0]
% als = [axls;ayls;azls]
% MATLAB implemtation of RTKLIB c source

% Notes:
%   Luni-solar accelerations considered constant during integration period
%h = ti - te;
% per GLONASS ICD,  integration period h must be less than 15 min
if abs(h)>900
    warning("h = "+num2str(h)+" > 15min");
end

%% Recalculate ephemeris parameters

x=reshape(x,numel(x),1);
k1=desat(x,als); w=x+k1*h/2;
k2=desat(w,als); w=x+k2*h/2;
k3=desat(w,als); w=x+k3*h;
k4=desat(w,als);
x=x+(k1+2*k2+2*k3+k4)*h/6;

    function xdot = desat(x,als)
        xdot=zeros(6,1);
        r2=dot(x(1:3),x(1:3)); r3=sqrt(r2)*r2; omge2=Const.OMEGAE_GLO^2;
        a=1.5*Const.J20_GLO*Const.MU_GLO*(Const.AE_GLO^2)/r2/r3;
        b=5*x(3)*x(3)/r2;
        c=-Const.MU_GLO/r3-a*(1.0-b);
        
        xdot(1)=x(4); xdot(2)=x(5); xdot(3)=x(6);
        xdot(4)=(c+omge2)*x(1)+2*Const.OMEGAE_GLO*x(5)+als(1);
        xdot(5)=(c+omge2)*x(2)-2*Const.OMEGAE_GLO*x(4)+als(2);
        xdot(6)=(c-2*a)*x(3)+als(3);
    end
end