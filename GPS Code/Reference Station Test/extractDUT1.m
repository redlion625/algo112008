function [DUT1,DUT1err] = extractDUT1(utc)
% extractDUT1
% Given UTC in Modified Julian Days
% 
% Returned (UT1-UTC)

load('aeroiersdata20220131.mat','-mat','mjd','ut1utc','ut1utcerr')
ind=find(mjd==utc);
DUT1=ut1utc(ind);
DUT1err=ut1utcerr(ind);
