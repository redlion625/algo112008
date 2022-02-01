function mjd = greg2mjd(y,m,d,hh,mm,ss)
mjd=jd2mjd(greg2jd(y,m,d,hh,mm,ss));