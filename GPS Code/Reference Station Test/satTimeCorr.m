function tsv=satTimeCorr(p,tr,tref,a,rootA,Ek)
% https://gssc.esa.int/navipedia/index.php/Clock_Modelling
% times in seconds of the week of a consistent time scale
% tref - time of reference epoch
% tr - time of receiver
% p = pseudorange
% a - polynomial coefficients, max 2nd order

% e - orbit eccentricity
% rootA - sqrt of orbit semi-major axis
% Ek - orbit eccentric anomaly

% Unlike in GPS, relativistic corrections to GLONASS orbits eccentricity are 
% transmitted within the navigation message into the satellite clock corrections (τn, γn).
% Thence, (2) is not needed with such broadcast message.

% In GLONASS, as the message is updated every 1/2h, only a first order
% polynomial is considered, being (a0=−τn) the clock offset and (a1=γn)
% the relative frequency offset, see table 2. Note that in the RINEX files 
% the −τn value is given instead of τn.
c=299792458;
tsv0=tr-p/c;

dtref=t-tref;
dtpoly=a(1)+a(2)*dtref+a(3)*dtref^2; % polynomial correction due to satellite and receiver clock differences
dtrel=-4.442807633e-10*e*rootA*sin(Ek); % relativistic correction due to relativity

tsv=tsv0-dtpoly-dtrel;



