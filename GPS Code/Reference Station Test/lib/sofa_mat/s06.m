function s = s06(date1,date2,x,y)
% % Time since J2000.0, in Julian centuries */
%    double t;
% 
% % Miscellaneous */
%    int i, j;
%    double a, w0, w1, w2, w3, w4, w5;
% 
% % Fundamental arguments */
%    double fa[8];
% 
% % Returned value */
%    double s;

% --------------------- */
% The series for s+XY/2 */
% --------------------- */
% 
%    typedef struct {
%       int nfa[8];      % coefficients of l,l',F,D,Om,LVe,LE,pA */
%       double s, c;     % sine and cosine coefficients */
%    } TERM;

% Polynomial coefficients */
sp = [ ...

   % 1-6 */
          94.00e-6, ...
        3808.65e-6, ...
        -122.68e-6, ...
      -72574.11e-6, ...
          27.98e-6, ...
          15.62e-6 ...
   ];

% Terms of order t^0 */
   %static const TERM s0[] = [
   s0=[ ...

   % 1-10 */
      0,  0,  0,  0,  1,  0,  0,  0, -2640.73e-6,   0.39e-6; ...
      0,  0,  0,  0,  2,  0,  0,  0,   -63.53e-6,   0.02e-6; ...
      0,  0,  2, -2,  3,  0,  0,  0,   -11.75e-6,  -0.01e-6; ...
      0,  0,  2, -2,  1,  0,  0,  0,   -11.21e-6,  -0.01e-6; ...
      0,  0,  2, -2,  2,  0,  0,  0,     4.57e-6,   0.00e-6; ...
      0,  0,  2,  0,  3,  0,  0,  0,    -2.02e-6,   0.00e-6; ...
      0,  0,  2,  0,  1,  0,  0,  0,    -1.98e-6,   0.00e-6; ...
      0,  0,  0,  0,  3,  0,  0,  0,     1.72e-6,   0.00e-6; ...
      0,  1,  0,  0,  1,  0,  0,  0,     1.41e-6,   0.01e-6; ...
      0,  1,  0,  0, -1,  0,  0,  0,     1.26e-6,   0.01e-6; ...

   % 11-20 */
      1,  0,  0,  0, -1,  0,  0,  0,     0.63e-6,   0.00e-6; ...
      1,  0,  0,  0,  1,  0,  0,  0,     0.63e-6,   0.00e-6; ...
      0,  1,  2, -2,  3,  0,  0,  0,    -0.46e-6,   0.00e-6; ...
      0,  1,  2, -2,  1,  0,  0,  0,    -0.45e-6,   0.00e-6; ...
      0,  0,  4, -4,  4,  0,  0,  0,    -0.36e-6,   0.00e-6; ...
      0,  0,  1, -1,  1, -8, 12,  0,     0.24e-6,   0.12e-6; ...
      0,  0,  2,  0,  0,  0,  0,  0,    -0.32e-6,   0.00e-6; ...
      0,  0,  2,  0,  2,  0,  0,  0,    -0.28e-6,   0.00e-6; ...
      1,  0,  2,  0,  3,  0,  0,  0,    -0.27e-6,   0.00e-6; ...
      1,  0,  2,  0,  1,  0,  0,  0,    -0.26e-6,   0.00e-6; ...

   % 21-30 */
      0,  0,  2, -2,  0,  0,  0,  0,     0.21e-6,   0.00e-6; ...
      0,  1, -2,  2, -3,  0,  0,  0,    -0.19e-6,   0.00e-6; ...
      0,  1, -2,  2, -1,  0,  0,  0,    -0.18e-6,   0.00e-6; ...
      0,  0,  0,  0,  0,  8,-13, -1,     0.10e-6,  -0.05e-6; ...
      0,  0,  0,  2,  0,  0,  0,  0,    -0.15e-6,   0.00e-6; ...
      2,  0, -2,  0, -1,  0,  0,  0,     0.14e-6,   0.00e-6; ...
      0,  1,  2, -2,  2,  0,  0,  0,     0.14e-6,   0.00e-6; ...
      1,  0,  0, -2,  1,  0,  0,  0,    -0.14e-6,   0.00e-6; ...
      1,  0,  0, -2, -1,  0,  0,  0,    -0.14e-6,   0.00e-6; ...
      0,  0,  4, -2,  4,  0,  0,  0,    -0.13e-6,   0.00e-6; ...

   % 31-33 */
      0,  0,  2, -2,  4,  0,  0,  0,     0.11e-6,   0.00e-6; ...
      1,  0, -2,  0, -3,  0,  0,  0,    -0.11e-6,   0.00e-6; ...
      1,  0, -2,  0, -1,  0,  0,  0,    -0.11e-6,   0.00e-6  ...
   ];

% Terms of order t^1 */
   %static const TERM s1[] = {
s1 = [ ...
   % 1 - 3 */
      0,  0,  0,  0,  2,  0,  0,  0,    -0.07e-6,   3.57e-6; ...
      0,  0,  0,  0,  1,  0,  0,  0,     1.73e-6,  -0.03e-6; ...
      0,  0,  2, -2,  3,  0,  0,  0,     0.00e-6,   0.48e-6  ...
   ];

% Terms of order t^2 */
   %static const TERM s2[] = {
   s2 = [ ...

   % 1-10 */
      0,  0,  0,  0,  1,  0,  0,  0,   743.52e-6,  -0.17e-6; ...
      0,  0,  2, -2,  2,  0,  0,  0,    56.91e-6,   0.06e-6; ...
      0,  0,  2,  0,  2,  0,  0,  0,     9.84e-6,  -0.01e-6; ...
      0,  0,  0,  0,  2,  0,  0,  0,    -8.85e-6,   0.01e-6; ...
      0,  1,  0,  0,  0,  0,  0,  0,    -6.38e-6,  -0.05e-6; ...
      1,  0,  0,  0,  0,  0,  0,  0,    -3.07e-6,   0.00e-6; ...
      0,  1,  2, -2,  2,  0,  0,  0,     2.23e-6,   0.00e-6; ...
      0,  0,  2,  0,  1,  0,  0,  0,     1.67e-6,   0.00e-6; ...
      1,  0,  2,  0,  2,  0,  0,  0,     1.30e-6,   0.00e-6; ...
      0,  1, -2,  2, -2,  0,  0,  0,     0.93e-6,   0.00e-6; ...

   % 11-20 */
      1,  0,  0, -2,  0,  0,  0,  0,     0.68e-6,   0.00e-6; ...
      0,  0,  2, -2,  1,  0,  0,  0,    -0.55e-6,   0.00e-6; ...
      1,  0, -2,  0, -2,  0,  0,  0,     0.53e-6,   0.00e-6; ...
      0,  0,  0,  2,  0,  0,  0,  0,    -0.27e-6,   0.00e-6; ...
      1,  0,  0,  0,  1,  0,  0,  0,    -0.27e-6,   0.00e-6; ...
      1,  0, -2, -2, -2,  0,  0,  0,    -0.26e-6,   0.00e-6; ...
      1,  0,  0,  0, -1,  0,  0,  0,    -0.25e-6,   0.00e-6; ...
      1,  0,  2,  0,  1,  0,  0,  0,     0.22e-6,   0.00e-6; ...
      2,  0,  0, -2,  0,  0,  0,  0,    -0.21e-6,   0.00e-6; ...
      2,  0, -2,  0, -1,  0,  0,  0,     0.20e-6,   0.00e-6; ...

   % 21-25 */
      0,  0,  2,  2,  2,  0,  0,  0,     0.17e-6,   0.00e-6; ...
      2,  0,  2,  0,  2,  0,  0,  0,     0.13e-6,   0.00e-6; ...
      2,  0,  0,  0,  0,  0,  0,  0,    -0.13e-6,   0.00e-6; ...
      1,  0,  2, -2,  2,  0,  0,  0,    -0.12e-6,   0.00e-6; ...
      0,  0,  2,  0,  0,  0,  0,  0,    -0.11e-6,   0.00e-6 ...
   ];

% Terms of order t^3 */
   %static const TERM s3[] = {
s3 = [ ...
   % 1-4 */
      0,  0,  0,  0,  1,  0,  0,  0,     0.30e-6, -23.42e-6; ...
      0,  0,  2, -2,  2,  0,  0,  0,    -0.03e-6,  -1.46e-6; ...
      0,  0,  2,  0,  2,  0,  0,  0,    -0.01e-6,  -0.25e-6; ...
      0,  0,  0,  0,  2,  0,  0,  0,     0.00e-6,   0.23e-6  ...
   ];

% Terms of order t^4 */
   %static const TERM s4[] = {
   s4 = [ ...

   % 1-1 */
      0,  0,  0,  0,  1,  0,  0,  0,    -0.26e-6,  -0.01e-6 ...
   ];

% Number of terms in the series */
%    static const int NS0 = (int) (sizeof s0 / sizeof (TERM));
%    static const int NS1 = (int) (sizeof s1 / sizeof (TERM));
%    static const int NS2 = (int) (sizeof s2 / sizeof (TERM));
%    static const int NS3 = (int) (sizeof s3 / sizeof (TERM));
%    static const int NS4 = (int) (sizeof s4 / sizeof (TERM));

NS0 = size(s0,1);
NS1 = size(s1,1);
NS2 = size(s2,1);
NS3 = size(s3,1);
NS4 = size(s4,1);

% ------------------------------------------------------------------ */

% Interval between fundamental epoch J2000.0 and current date (JC). */
   t = ((date1 - Const.DJ00) + date2) / Const.DJC;

% Fundamental Arguments (from IERS Conventions 2003) */

fa=zeros(8,1);

% Mean anomaly of the Moon. */
   fa(1) = fal03(t);

% Mean anomaly of the Sun. */
   fa(2) = falp03(t);

% Mean longitude of the Moon minus that of the ascending node. */
   fa(3) = faf03(t);

% Mean elongation of the Moon from the Sun. */
   fa(4) = fad03(t);

% Mean longitude of the ascending node of the Moon. */
   fa(5) = faom03(t);

% Mean longitude of Venus. */
   fa(6) = fave03(t);

% Mean longitude of Earth. */
   fa(7) = fae03(t);

% General precession in longitude. */
   fa(8) = fapa03(t);

% Evaluate s. */

   w0 = sp(1);
   w1 = sp(2);
   w2 = sp(3);
   w3 = sp(4);
   w4 = sp(5);
   w5 = sp(6);

%    for (i = NS0-1; i >= 0; i--) {
%    a = 0.0;
%    for (j = 0; j < 8; j++) {
%       a += (double)s0[i].nfa[j] * fa[j];
%    }
%    w0 += s0[i].s * sin(a) + s0[i].c * cos(a);
%    }

for i = NS0:-1:1
    a=0;
    for j = 1:8
        a = a + s0(i,j) * fa(j);
    end
    w0 = w0 + s0(i,9)*sin(a) + s0(i,10)*cos(a);
end
%     
%    for (i = NS1-1; i >= 0; i--) {
%       a = 0.0;
%       for (j = 0; j < 8; j++) {
%          a += (double)s1[i].nfa[j] * fa[j];
%       }
%       w1 += s1[i].s * sin(a) + s1[i].c * cos(a);
%    }

for i = NS1:-1:1
    a=0;
    for j = 1:8
        a = a + s1(i,j)*fa(j);
    end
    w1 = w1 + s1(i,9)*sin(a) + s1(i,10)*cos(a);
end

%    for (i = NS2-1; i >= 0; i--) {
%       a = 0.0;
%       for (j = 0; j < 8; j++) {
%          a += (double)s2[i].nfa[j] * fa[j];
%       }
%       w2 += s2[i].s * sin(a) + s2[i].c * cos(a);
%    }

for i = NS2:-1:1
    a=0;
    for j = 1:8
        a = a + s2(i,j)*fa(j);
    end
    w2 = w2 + s2(i,9)*sin(a) + s2(i,10)*cos(a);
end

%    for (i = NS3-1; i >= 0; i--) {
%       a = 0.0;
%       for (j = 0; j < 8; j++) {
%          a += (double)s3[i].nfa[j] * fa[j];
%       }
%       w3 += s3[i].s * sin(a) + s3[i].c * cos(a);
%    }

   
for i = NS3:-1:1
    a=0;
    for j = 1:8
        a = a + s3(i,j)*fa(j);
    end
    w3 = w3 + s3(i,9)*sin(a) + s3(i,10)*cos(a);
end
   
%    for (i = NS4-1; i >= 0; i--) {
%       a = 0.0;
%       for (j = 0; j < 8; j++) {
%          a += (double)s4[i].nfa[j] * fa[j];
%       }
%       w4 += s4[i].s * sin(a) + s4[i].c * cos(a);
%    }

for i = NS4:-1:1
    a=0;
    for j = 1:8
        a = a + s4(i,j)*fa(j);
    end
    w4 = w4 + s4(i,9)*sin(a) + s4(i,10)*cos(a);
end

   s = (w0 + ...
       (w1 + ...
       (w2 + ...
       (w3 + ...
       (w4 + ...
        w5 * t) * t) * t) * t) * t) * Const.DAS2R - x*y/2.0;

%   return s;
