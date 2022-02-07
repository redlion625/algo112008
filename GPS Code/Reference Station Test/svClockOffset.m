function dt = svClockOffset(t1,t0,a0,a1,a2)
% svClockOffset(t1,t0,a0,a1,a2)
% Given:
%     a0,a1,a2    2nd order SV clock offset polynomial terms from ephemeris
%     t1          uncorrected transmission time
%     t0          time of ephemeris epoch
%
% Returned: 
%     dt          satellite clock offset correction
%
% Notes:
%   For GLONASS, -Taun = a0, Gamman = a1 and a2 = 0
dt = a0 + a1*(t1-t0) + a2*(t1-t0)^2;
