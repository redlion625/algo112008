function rbpn = pnm06a(date1,date2)


% Fukushima-Williams angles for frame bias and precession.
[gamb, phib, psib, epsa]=pfw06(date1, date2);

% Nutation components.
[dp, de]=nut06a(date1, date2);

% Equinox based nutation x precession x bias matrix.
rbpn=fw2m(gamb, phib, psib + dp, epsa + de);
