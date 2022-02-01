function [dpsi,deps]=nut06a(date1,date2)
t = ((date1 - Const.DJ00) + date2) / Const.DJC;

% Factor correcting for secular variation of J2. 
fj2 = -2.7774e-6 * t;

% Obtain IAU 2000A nutation. 
[dp,de]=nut00a(date1, date2);

% Apply P03 adjustments (Wallace & Capitaine, 2006, Eqs.5). 
dpsi = dp + dp * (0.4697e-6 + fj2);
deps = de + de * fj2;
