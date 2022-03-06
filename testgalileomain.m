%% ESSE 3670 driver file
clear all; 
close all;
clc;
format long g


%Reading the observation and navigation data into the program
Nav_data = Nav_readergal("D:\Third Year\ESSE 3670\Project 3\Data Downloaded\datasets to try with\algo\algo112008\brdm1190.21p"); %Use brdm1180.21p    %brdm1190.21p %brdm0840.21p  %brdm0840.21p  %brdm2260.20p  %brdm1180.21p  %brdm1560.20p  %brdm1180.21p"); %brdm1360.20p
Obs_data = testgalobsreader("D:\Third Year\ESSE 3670\Project 3\Data Downloaded\datasets to try with\algo\algo112008\Proper0429SamsungS20Ultra_GnssLog.21o"); %CorrectSamsungS20Ultra_GnssLog0428.21o  %Use SamsungS20Ultra_GnssLog0429.21o   % SamsungS20Ultra_GnssLog0429.21o  %SamsungS20Ultra_GnssLog0325.21o   %Pixel4XL_GnssLog0604.20o  %Mi8_GnssLog.20o  % %SamsungS20Ultra_GnssLog.21o  %Pixel4_GnssLog.21o     %Pixel4_GnssL.20o  Pixel4_GnssLog.21o"); %Pixel4_GnssLog.20o" %Pixel4_GnssLog.21o %39ea118x.21o


%% Match observation and navgitation data
Observations = Obs_data.data;
obssize = size(Observations{1,2}); %Determining the number of observations
navsize = size(Nav_data.data); %Determining the number of navigation

totalnum = 0; 
%Iterating through the observations vector
for i = 1:length(Observations)
    obssize = size(Observations{i,2}); %Entering the current submatrix for the current observation epoch
    totalnum = totalnum + 1; 
    for j = 1:obssize(1) %Iterating through the number of satellites at each epoch 
        timediff = [];
        I = [];
        satnum = Observations{i,2}(j,2); %The current satellite number at the current epoch 
        navsize = size(Nav_data.data); %The total number of satellites at specific GPS time
        for k = 1:navsize(1) %Iterating through the navigation data satellite at a particular epoch observations
            if Nav_data.data(k,1) == satnum %If the satellite number in the navigation data is the same as that in the current satellite in the observation data
                I = [I,k]; %Saving the index at which the match was made
                timediff = [timediff,(Observations{i,1} - Nav_data.data(k,2))]; %Taking the difference between the navigation time block and the observation time for the same satellite
            end
        end
        [M,P] = min(abs(timediff)); %Determining the minimum difference in time between the navigation and observation time epochs for a particular satellite for matching
        testm(j,1) = M; 
        index = I(P); %Saving the minimum matched index for current reference
        Observations{i,2}(j,obssize(2)+1) = index; %Concatenating the best satellite epoch and their corresponding satellite information to the observation epoch it was matched to in a submatrix
    end
end

%% Loop through each satellite
f1 = 1575.42e6; %L1
f2 = 1176.45e6; %L5
GPSA = Nav_data.GPSA; 
GPSB = Nav_data.GPSB;
gamma = ((f1/f2)^2);
gravitationalparameter = 3.986004418e14; %3.986005e14; %WGS84 value of the Earth's gravitational constant for GPS
earthrotation = 7.2921151467e-5; %WGS84 value of the Earth's rotation rate
c = 2.99792458e8; %Speed of light
for i = 1:length(Observations) %Iterating through the number of observations
    obssize = size(Observations{i,2}); %Determining the size of the submatrix attached to the epoch of observations
    clear elevationangle azimuth dTsv
     remove = [];
     if(i == Obs_data.sum - 1)
         break;
     else
         
    for j = 1:obssize(1) %Iterating through the submatrix attached to the observation epoch
        P2 = Observations{i,2}(j,7); %P2 code
        P1 = Observations{i,2}(j,3); %P1 code
        if(P1 ~= 0)
            if(P2 ~= 0)
                Pseudorange(i,j) = ((f1^2)/((f1^2) - (f2^2)))*(P1) - ((f2^2)/((f1^2) - (f2^2)))*(P2); % ((f1^2)*P1 - (f2^2)*P2)/(f1^2-f2^2);
                Pseudoranget(i,j) = (P2 - gamma*P1)/(1 - gamma);
                Pseudorangeagain(i,j) = ((f1^2)*P1 - (f2^2)*P2)/(f1^2-f2^2);
                diffpseudo(i,j) = P1 - Pseudorange(i,j); 
                testpseudo(i,j) = Pseudorange(i,j) - Pseudoranget(i,j);
                uo = 1;
                indL1 = Observations{i,2}(j,6);
                indL2 = Observations{i,2}(j,10);
                if(indL1 > indL2)
                    indexL = Observations{i,2}(j,6); %indexL(i,j)
                    bL = 293;
                else
                    indexL = Observations{i,2}(j,10); %indexL(i,j)
                    bL = 29.3;
                end
            else
                Pseudorange(i,j)= P1; 
                indexL = Observations{i,2}(j,6); %indexL(i,j)
                bL = 293;
            end
        end
        if(P2 ~= 0)
            if(P1 == 0)
                Pseudorange(i,j)= P2; 
                indexL = Observations{i,2}(j,10); %indexL(i,j)
                bL = 29.3;
            end
        end        
        
        %Pseudorange(i,j) = Observations{i,2}(j,3); %((f1^2)*P1 - (f2^2)*P2)/(f1^2-f2^2); %Calculating the uncorrected pseudorange
        index = Observations{i,2}(j,11); 
        sqrtA = Nav_data.data(index, 13);  %Finding the semi-major axis
        A = (sqrtA)^2; %Check %Finding the semi-major axis
        N0 = sqrt(gravitationalparameter/A^3); %Check %Computing the mean motion
        Toe = Nav_data.data(index, 14); %Reading in the ephemeris data reference time of week
        Toc = Nav_data.data(index, 24)*60*60*24*7 + Nav_data.data(index, 14); %Finding the clock data reference time of week
        T = Pseudorange(i,j)/c; %Calculating the travel time from the uncorrected pseudorange and the speed of light
        %%Not sure if T is applicable and Trec along with time
        time = cell2mat(Observations(i,1)); 
        Trec = time - T; 
        Tk = Trec - Toc; %Calculating the time of ephemeris reference time
        %Checking that the time of ephemeris is between -302400 and 302400
        if Tk > 302400
            Tk = Tk - 604800;
        end
        if Tk < -302400
            Tk = Tk + 604800;
        end
        N = N0 + Nav_data.data(index, 8); %The corrected mean motion
        Mk = Nav_data.data(index, 9) + N*Tk; %Calculating the mean anomaly
        ecc = Nav_data.data(index, 11); %Reading in the eccentricity
        count = 1;
        %Determining the eccentric anomaly value from Kepler's equation for
        %Eccentric Anomaly. The Ek will first be approximated by Mk. The Ek
        %will then be found through a series of iterations
        Ek0 = Mk;
        diff = 1;
        Ek = Ek0 + ecc*sin(Ek0);
        diff = abs(Ek-Ek0);
        while count <25 && diff > 1e-12
            Ek = Mk + ecc*sin(Ek0);
            diff = abs(Ek-Ek0);
            Ek0 = Ek;
            count = count +1;
        end 
%         Ekexp = (ecc + (cos(Vk(i,j))))/(1 + ecc*(cos(Vk(i,j))));
%         Ektest = acos(Ekexp);         
%         diffEktest = Ek - Ektest;
        Toc = Nav_data.data(index, 24)*60*60*24*7 + Nav_data.data(index, 14); %Determining the clock data reference time of week
        dTr = (-4.442807633e-10)*ecc*sqrtA*sin(Ek); %Ek %Relativistic clock correction term
        dTsv(j) = Nav_data.data(index, 3) + Nav_data.data(index, 4)*(time-Toc) + Nav_data.data(index, 5)*((time-Toc)^2) + dTr; %Determining the satellite clock offset
        Vk(i,j) = atan2((sqrt(1-ecc^2)*sin(Ek)/(1-ecc*cos(Ek))),((cos(Ek)-ecc)/(1-ecc*cos(Ek)))); %Ek %Determining the true anomaly
        ArgLat(1) = Vk(i,j) + Nav_data.data(index, 20); %Finding the argument of latitude
        %The second harmonic perturbations
        ArgLat(2) = Nav_data.data(index, 12)*sin(2*ArgLat(1)) + Nav_data.data(index, 10)*cos(2*ArgLat(1)); %The argument of latitude correction
        ArgLat(3) = Nav_data.data(index, 7)*sin(2*ArgLat(1)) + Nav_data.data(index, 19)*cos(2*ArgLat(1)); %The radial correction
        ArgLat(4) = Nav_data.data(index, 17)*sin(2*ArgLat(1)) + Nav_data.data(index, 15)*cos(2*ArgLat(1)); %The inclination correction
        Corr(1) = ArgLat(1) + ArgLat(2); %The corrected argument of latitude
        Corr(2) = A*(1-ecc*cos(Ek))+ArgLat(3); %Ek %The corrected radius
        Corr(3) = Nav_data.data(index, 18) + Nav_data.data(index, 22)*Tk + ArgLat(4); %The corrected inclination angle
        %Determining the positions in the orbital plane
        Xkprime = Corr(2)*cos(Corr(1)); 
        Ykprime = Corr(2)*sin(Corr(1));
        
        OMEGAdot = Nav_data.data(index, 21); %Reading in the rate of right ascension
        OMEGAk = Nav_data.data(index, 16) + (OMEGAdot - earthrotation)*Tk - earthrotation*Toe; %Determining the corrected longitude of ascending node
        xk = Xkprime*cos(OMEGAk)-Ykprime*cos(Corr(3))*sin(OMEGAk); %The satellite's X component in ECEF coordinates
        yk = Xkprime*sin(OMEGAk)+ Ykprime*cos(Corr(3))*cos(OMEGAk); %The satellite's Y component in ECEF coordinates
        Zk = Ykprime*sin(Corr(3)); %The satellite's Z component in ECEF coordinates
        
        %Making corrections for the Earth's rotation
        theta = earthrotation * T;
         Xk = xk*cos(theta) + yk*sin(theta); %The X component in ECEF coordinates
         Yk = -xk*sin(theta) + yk*cos(theta); %The Y component in ECEF coordinates
         
        VectorXYZ = [Xk, Yk, Zk]; %The ECEF coordinates
        Approx = [Obs_data.pos_xyz(1,1) ; Obs_data.pos_xyz(1,2) ; Obs_data.pos_xyz(1,3)]; %The approximate position of the receiver
        Approxcoord = Approx'; 
        LLHtrans = EC2LLH(VectorXYZ); %Calculating the latitude, longitude, and height of the ECEF coordinates
        LLHtrans(1:2) = LLHtrans(1:2)*pi/180; %Converting the result to radians
        longit = LLHtrans(1); 
        latit = LLHtrans(2); 
        LLH = LLHtrans';
        ENU = EC2ENU(VectorXYZ, Approxcoord, LLH); %Determining the easting, northing, and up from the ECEF coordinates

        azimuth(j) = atan2(ENU(1,1),ENU(1,2)); %Determining the azimuth 
        slantd = sqrt(ENU(1)^2 + ENU(2)^2 + ENU(3)^2); %Finding the slant distance
        zenith = acos(ENU(3)/slantd); %Finding the zenith angle
        elevationangle(j) = pi/2 - zenith; %Determining the elevation angle from the zenith angle
        if(elevationangle(j) < (10/180*pi)) %If the elevation angle is less tha 10 degrees, store the index
           remove = [remove,j]; 
        end
        
        Observations{i,2}(j,18) = elevationangle(j); %Storing the elevation angle
        Observations{i,2}(j,19) = azimuth(j); %Storing the zenith angle

        %%ionospheric effect calculations
        
        if(P2 == 0)
            orbitalw = 0.0137/((Observations{i,2}(j,18)) + 0.11) - 0.022; %Calculating the earth's central angle between the user position and the earth projection of ionospheric intersection point
            %This is calculated using the elevation angle 
        
            thetai = latit + orbitalw*(cos(azimuth(j))); %Calculating the geodetic latitude of the earth projection fo the ionospheric intersection point 
            %This is calculated using the user geodetic latitude, the earth's
            %central angle, and the azimuth angle between the user and
            %satellite
        
            %If the geodetic latitude of the earth projection falls between
            %the constraints, it has to be corrected
            if(thetai > 0.416)
                thetai = 0.416; 
            elseif(thetai < -0.416)
                thetai = -0.416; 
            else
                thetai = thetai; 
            end
        
            longitudei = longit + (orbitalw*(sin(azimuth(j))))/(cos(thetai)); %Calculating the geodetic longitude of the earth projection of the ionospheric intersection point
            %This is calculated based on the user geodetic longitude in
            %radians, the azimuth, and the geodetic lattiude of the earth
            %projection of the ionospheric intersection point
            geomaglat = thetai + 0.064*(cos(longitudei - 1.617)); %Geomagnetic latitude of the earth projection of the ionospheric intersection point
            %This is calculated with the geodetic latitude of the earth projection and the geodetic longitude of the earth projection of the ionospheric intersection
        
            add = 0;
            %Calulating PER based on the ionospheric beta coefficients and the
            %geomagnetic magnitude
            for t = 1:4
                PER = (GPSB(t))*(geomaglat^(t-1));
                PER = PER + add; 
                add = PER;
            end
            if(PER < 72000)
                PER = 72000; 
            end
        
            timet = (4.32e4)*(longitudei) + Observations{j,1};
            x = (2*pi*(timet - 50400))/PER; 
        
            sumadd = 0;
            %Calulating AMP based on the ionospheric alpha coefficients and the
            %geomagnetic magnitude
            for t = 1:4
                AMP = (GPSA(t))*(geomaglat^(t-1));
                AMP = AMP + sumadd; 
                sumadd = AMP;
            end
            if(AMP < 0)
                AMP = 0; 
            end
        
            F = 1.0 + 16.0*((0.53 - (elevationangle(j))*(1/pi))^3); %Calculating the obliquity factor based on the elevation angle
        
            if(abs(x) >= 1.57)
                Tiono(i,j) = (3e8)*(F*(5.0e-9)); %Calculating the ionospheric correction model 
            else
                Texp = AMP*(1 - (x^2)/2 + (x^4)/24);
                Tiono(i,j) = (3e8)*(F*((5.0e-9) + Texp)); %Calculating the ionospheric correction model
            end
        
            corrPseudorange(i,j) = Pseudorange(i,j) + Tiono(i,j);
        else
            corrPseudorange(i,j) = Pseudorange(i,j); 
        end

         %Constants to determine the tropospheric error
        Pwv = 16; 
        pd = 1013 - Pwv; 
        Ttrop = 300; 
        %Calculating and storing the tropospheric error
        dtrop = (0.002277/(cos(zenith)))*(pd + (1255/Ttrop + 0.05)*Pwv - 1.16*((tan(zenith))));
        Observations{i,2}(j,20) = dtrop;
        
        %The corrected pseudorange
        correctedpseudorange(i,j) = corrPseudorange(i,j) - dtrop; % + Tiono(i,j);
        
        %Storing the corrected pseudorange, satellite clock offest, and
        %ECEF Xk, Yk, and Zk values in the overall Observations matrix for
        %each satellite at each epoch
        Observations{i,2}(j,12) = correctedpseudorange(i,j);
        Observations{i,2}(j,13) = dTsv(j);
        Observations{i,2}(j,14) = Xk;
        Observations{i,2}(j,15) = Yk;
        Observations{i,2}(j,16) = Zk;
        
        factorL = (-0.5*(indexL)/10);
        testpower = power(10,factorL);
        scalingfactor(i,j) = (bL)*testpower; %(-0.5*(factorL/10)));       
    end
    
    %Removing all satellites that are below 10 degrees
    for k = length(remove):-1:1
        Observations{i,2}(remove(k),:) = [];
    end
    
     check = isnan(Observations{i,2});
     for j = obssize:-1:1
         if check(j) == 1
             Observations{i,2} = [];
         end
     end
    
    obssize = size(Observations{i,2});
       
    %%Signal to noise ratio weighting
    
    %% Least squares portion
  
    %Clearing the vectors involved in the adjustment process
    L = []; 
    Xs = []; 
    Ys = []; 
    Zs = [];
    
    L = Observations{i,2}(:,12); %The current psuedorange measurement for the current satellite and epoch
    %The initial approximate of the user satellite position for the observation file
    xapprox = 0;  % Obs_data.pos_xyz(1,1); 
    yapprox = 0;   %Obs_data.pos_xyz(1,2); 
    zapprox = 0;   %Obs_data.pos_xyz(1,3);     
    
    Xs = Observations{i,2}(:,14); %The current x Earth-fixed coordinates of SV antenna phase center
    Ys = Observations{i,2}(:,15); %The current y Earth-fixed coordinates of SV antenna phase center
    Zs = Observations{i,2}(:,16); %The current z Earth-fixed coordinates of SV antenna phase center
        
    %Setting the initial approximations of the receiver positions and the
    %receiver clock offset
    xcurrent = xapprox;   ycurrent = yapprox;   zcurrent = zapprox;   dtrcurrent = 0;
    count = 0;   
    tab = 0;
    
    %Performing a nonlinear least squares adjustment to determine the least
    %squares estimates of the receiver x, y, z position and the receiver
    %clock offset
    while(count < 40)
        clear PR A mis
        x0 = xcurrent;   y0 = ycurrent;   z0 = zcurrent;  dtr0 = dtrcurrent; %Setting the current estimates of the unknown parameters to the current least squares estimates of the parameters
    
    for j = 1:obssize(1) %A for-loop to create the first design matrix
        %Determining the geometric range
        %x0 = xcurrent, y0 = ycurrent, z0 = zcurrent
        geometricexpression = (Xs(j) - xcurrent)^2 + (Ys(j)- ycurrent)^2 + (Zs(j) - zcurrent)^2;
        geometricinitial = sqrt(geometricexpression); %+c(dts-dtr)

        %The deriviatives for the X, Y, and Z receiver components and the
        %receiver clock offset
        derivx = -((Xs(j) - x0)/geometricinitial); 
        derivy = -((Ys(j) - y0)/geometricinitial); 
        derivz = -((Zs(j) - z0)/geometricinitial); 
        derivdtr = c;
 
        %Populating the first design matrix
        A(j,1) = derivx;
        A(j,2) = derivy;
        A(j,3) = derivz;
        A(j,4) = derivdtr;
        
        %Determining the current uncorrected pseudorange for the
        %approximatations of the receiver position and clock offset
        PR(j) = sqrt(geometricexpression) + c * (dtr0 - Observations{i,2}(j,13));  %- dTsv(j)
    end
        
        mis = L - PR'; %Determining the misclosure between the current estimate of the pseudorange and the corrected pseudorange
        %P = eye(obssize(1));
        
        %Determining the updating value (deltax) to the current
        %approximations of the receiver position and clock offset
        Inverseterm = (A')*A; 
        Inverse = inv(Inverseterm); 
        deltax = Inverse*(A')*mis; 
        xo = [xcurrent; ycurrent; zcurrent; dtr0];  %[x0; y0; z0; dtr0]  %The current approximations of the receiver position and clock offset
        xhat = xo + deltax; %The current least squares estimates of the receiver position and clock offset
        
        %Setting the current approximations of the receiver position and clock offset to the current least squares estimates 
        xcurrent = xhat(1); 
        ycurrent = xhat(2); 
        zcurrent = xhat(3); 
        dtrcurrent = xhat(4); 
        disp(j); 
        thresh = [0.01; 0.01; 0.01; 1e-4]; %Determining the threshold at which the least squares loop with break
        count = count + 1; %Keeping track of how many iterations the loop makes
        if(all(abs(deltax) < thresh)) %If the delta values are less than the prescribed threshold values, the loop will break
            fprintf("Loop is broken");
            resultcheck(i) = count;
            v = A*deltax - mis; %Calculating the residuals
            Final = xhat; %Storing the least squares estimates of the X, Y, and Z components of the receiver's position and the receiver clock offset
            Observations{i,2}(:,17) = v; %Storing the residuals 
            aposteriori = ((v')*v)/(j - 4); %Determining the a-posteriori variance factor
%             Final(i,2) = aposteriori*Inverse;
            break; 
        end
%         if(all(abs(deltax) < thresh)) %If the delta values are less than the prescribed threshold values, the loop will break
%             fprintf("Loop is broken");
%             resultcheck(i) = count; 
%             Final = xhat; %Storing the least squares estimates of the X, Y, and Z components of the receiver's position and the receiver clock offset
%             v = A*deltax - mis; %Calculating the residuals
%             Observations{i,2}(:,17) = v; %Storing the residuals 
%             aposteriori = ((v')*v)/(j - 4); %Determining the a-posteriori variance factor
% %             Final(i,2) = aposteriori*Inverse;
%             break; 
%         end
     end
    
    DOP = inv(A'*A); %Determining the DOP values 
    OverallDOP = diag(DOP); 
    DOP1 = OverallDOP(1); %Finding the accuracy of the Xu component
    DOP2 = OverallDOP(2); %Finding the accuracy of the Yu component
    DOP3 = OverallDOP(3); %Finding the accuracy of the Zu component
    DOP4 = OverallDOP(4); %Finding the accuracy of the time (t) component
    GDOP(i) = sqrt(DOP1 + DOP2 + DOP3 + DOP4); %Determining the GDOP
    PDOP(i) = sqrt(DOP1 + DOP2 + DOP3); %Finding the PDOP
    TDOP(i) = sqrt(DOP4); %Determining the TDOP
    HDOP(i) = sqrt(DOP1 + DOP2); %Finding the HDOP 
    VDOP(i) = sqrt(DOP3); %Determining the VDOP
    
    %Converting the final X, Y, and Z components of the satellite position
    %into Easting, Northing, and Up
    TrueXYZ = [Final(1:3)];
    
    TrueTestENU = EC2ENU(TrueXYZ); 
    
    TrueLLHtrans = EC2LLH(TrueXYZ);
    XYZTable = EC2LLH(TrueXYZ);
    Currentepoch = Observations{i,1};
    CombinedTable(i,:) = [XYZTable(1); XYZTable(2); XYZTable(3); Currentepoch];
    XYZComparisonM(i,:) = [TrueTestENU(1); TrueTestENU(2); TrueTestENU(3); Currentepoch];  % [TrueXYZ(1); TrueXYZ(2); TrueXYZ(3); Currentepoch]; 
    
    
    TrueLLHtrans(1:2) = TrueLLHtrans(1:2)*pi/180;
    TrueLLH = TrueLLHtrans';
    TrueENU(i, :) = EC2ENU(TrueXYZ', 0, TrueLLH);     %EC2ENU(TrueXYZ', Approxcoord, TrueLLH);
    
     end
end
    
VT = table('Size',[size(CombinedTable)], 'VariableTypes', {'double', 'double', 'double', 'double'}, 'VariableNames',{'Lat','Long','Height','Epoch'});
VT{:,:}=CombinedTable;
%writetable(VT,'residual_stats.csv');

referenceresults = readmatrix('0429ground_truth.csv'); 

referencecoordinates = referencecoordinates(referenceresults);

measuredcoordinates = table2array(VT);

for d = 1:length(measuredcoordinates)
    latdiff(d,1) = abs(referencecoordinates(d,1) - measuredcoordinates(d,1)); 
    longdiff(d,1) = abs(referencecoordinates(d,2) - measuredcoordinates(d,2)); 
end


convertreference = zeros(length(measuredcoordinates),4);
for d = 1:length(measuredcoordinates)
    convertreference(d,1:3) = lla2ecef([referencecoordinates(d,1),referencecoordinates(d,2),referencecoordinates(d,3)]);
    convertreference(d,4) = referencecoordinates(d,4); 
end



for d = 1:length(measuredcoordinates)
    errorinpositionX(d) = convertreference(d,1) - XYZComparisonM(d,1); %measuredcoordinates(d,1); 
    errorinpositionY(d) = convertreference(d,2) - XYZComparisonM(d,2); %measuredcoordinates(d,2); 
    errorinpositionZ(d) = convertreference(d,3) - XYZComparisonM(d,3); %measuredcoordinates(d,3); 
    errorinpositionT(d) = convertreference(d,4) - XYZComparisonM(d,4); %measuredcoordinates(d,4); 
end

errorinpositionX = errorinpositionX.'; 
errorinpositionY = errorinpositionY.'; 
errorinpositionZ = errorinpositionZ.'; 
errorinpositionT = errorinpositionT.'; 


f7 = figure; 
plot(errorinpositionX);
hold on
title('Error In X Position');
ylabel('Error (Metres)'); 
xlabel('Epoch Number');
hold off


f8 = figure; 
plot(errorinpositionY);
hold on
title('Error In Y Position');
ylabel('Error (Metres)'); 
xlabel('Epoch Number');
hold off


f8 = figure; 
plot(errorinpositionZ);
hold on
title('Error In Z Position');
ylabel('Error (Metres)'); 
xlabel('Epoch Number');
hold off

%M = groundtruth; %readmatrix('ground_truth.csv'); 



%Holds all of azimuth, elevation, and tropospheric data and sorts them by satellite rather
%than by time for plotting purposes
lengthplot = length(errorinpositionX);
for i = 1:lengthplot    %2880
	obs = size(Observations{i,2});
	for j = 1:obs(1)
		sat = Observations{i,2}(j,2);
		Observations{sat,3}(i) = Observations{i,2}(j,18);
		Observations{sat,4}(i) = Observations{i,2}(j,19);
        Observations{sat,5}(i) = Observations{i,2}(j,20);
        Observations{sat,7}(i) = Observations{i,2}(j,17);
    end
    temp = size(Observations{i,2});
    Observations{i,6} = temp(1);
end

% i = 0;
% for i = 1:Obs_data.sum  %1748     %2880
% 	obs = size(Observations{i,2});
% 	for j = 1:obs(1)
% 		sat = Observations{i,2}(j,2);
% 		Observations{sat,3}(i) = Observations{i,2}(j,18);
% 		Observations{sat,4}(i) = Observations{i,2}(j,19);
%         Observations{sat,5}(i) = Observations{i,2}(j,20);
%         Observations{sat,7}(i) = Observations{i,2}(j,17);
%     end
%     temp = size(Observations{i,2});
%     Observations{i,6} = temp(1);
% end

%The code here is purely for plotting purposes. Any alterations are just to
%represent the results more clearly and accurately.
%Steps through each of the satellites at each epoch
for i = 1:32  %32   %1748
    size = length(Observations{i,3});
    for j = size:-1:1 %Remove the zero values from the azimuth and elevations
        if Observations{i,3}(j) == 0
            Observations{i,3}(j) = [];
            Observations{i,4}(j) = [];
        end
    end
    diff = 2880-size;
    Observations{i,5} = [Observations{i,5},zeros(1,diff)]; %Filling out each satellite observation
    for j = 1:2880  %Inserting NaN where there would be zeros in the tropospheric error for plotting purposes
        if Observations{i,5}(j) == 0
            Observations{i,5}(j) = NaN;
        end
    end
    size = length(Observations{i,7});
    diff = 2880-size;
    Observations{i,7} = [Observations{i,7},zeros(1,diff)]; %Filling out each residual observation
    for j = 1:2880  %Inserting NaN where there would be zeros
        if Observations{i,7}(j) == 0
            Observations{i,7}(j) = NaN;
        end
    end
end
%% Plots

%Plotting the number of satellites over time
satquan = cell2mat(Observations(:,6));
f1 = figure; 
plot(satquan);
title('The Number Of Satellites Versus Epoch');
ylabel('Number Of Satellites'); 
xlabel('Epoch Number');


%Plotting the azimuth angle vs elevation angle
f2 = figure; 
scatter(Observations{1,4},Observations{1,3},20,'.')
hold on
title('The Azimuth Angle Versus The Elevation Angle');
ylabel('Elevation Angle (Radians)'); 
xlabel('Azimuth Angle (Radians)');
for i = 2:32  %1748   %32
    scatter(Observations{i,4},Observations{i,3},20,'filled')
end
hold off

%Plotting the tropospheric error versus epoch
f3 = figure; 
plot(Observations{1,5});
hold on
title('Tropospheric Error Per Epoch');
ylabel('Error (Metres)'); 
xlabel('Epoch Number');
for i = 2:32  %1748  %32
    plot(Observations{i,5});
end
hold off

%Plotting the error as the difference between the least squares estimated
%and initial approximation of the reference receiver position
f4 = figure;
subplot(2,2,1)
plot(TrueENU(:,1),'.'); grid on;
title({
    ['User E Position Error With Respect']
    ['To The Reference Receiver Position']
    });
ylabel('Error (Metres)');
xlabel('Epoch Number');
subplot(2,2,2)
plot(TrueENU(:,2),'.'); grid on;
title({
    ['User N Position Error With Respect']
    ['To The Reference Receiver Position']
    });
ylabel('Error (Metres)');
xlabel('Epoch Number');
subplot(2,2,3)
plot(TrueENU(:,3),'.'); grid on;
title({
    ['User U Position Error With Respect']
    ['To The Reference Receiver Position']
    });
ylabel('Error (Metres)');
xlabel('Epoch Number');


%Plotting The DOP and DOP components
f5 = figure;
plot(GDOP);
hold on
title('Dilution Of Precision (DOP)');
ylabel('DOP (Metres)'); 
xlabel('Epoch Number)');
plot(HDOP);
plot(VDOP);
plot(TDOP);
plot(PDOP);
legend('GDOP', 'HDOP', 'VDOP', 'TDOP', 'PDOP');
hold off

%Plot residuals
f6 = figure; 
plot(Observations{1,7});
hold on
title('Residuals From The Least Squares Computations');
ylabel('Residuals (Metres)'); 
xlabel('Epoch Number');
countn = 0;  
 for i = 2:32   %32
         plot(Observations{i,7});
 end
% title('Residuals From The Least Squares Computations');
% ylabel('Residuals (Metres)'); 
% xlabel('Epoch Number');
hold off

exportgraphics(f1,'testsatnum.png')
exportgraphics(f2,'testlocalcoord.png')
exportgraphics(f3,'testtroperr.png')
exportgraphics(f4,'testenuerror.png')
exportgraphics(f5,'testdop.png')
exportgraphics(f6,'testv.png')
    
    