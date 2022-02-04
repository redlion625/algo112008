function gmst = gmstGLO(y,m,d,hh,mm,ss)
% gmstGLO(y,m,d,hh,mm,ss)
% Given:
%     y, m, d, hh, mm, ss
%
% Returned:
%     Greenwich mean sidereal time in radians
% 
% Reference:
%     GLONASS CDMA ICD (2016)
jd0 = greg2jd(y,m,d,hh,mm,ss);

tdelta = (jd0-Const.DJ00)/Const.DJC;
gmst = eraGLO(jd0) + 0.0000000703270726          + 0.0223603658710194*tdelta + ...
                   + 0.0000067465784654*tdelta^2 - 0.0000000000021332*tdelta^3 - ...
                   - 0.0000000001452308*tdelta^4 - 0.0000000000001784*tdelta^5;

