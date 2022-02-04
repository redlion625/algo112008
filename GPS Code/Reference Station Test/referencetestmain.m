%% ENG4000 GNSS Processing Main

% Column names of per epoch glonass array
% 1-5
% PRN,  Code,   Phase,  Doppler,    Signal Strength,    
%6-8
% Xsat,Ysat,Zsat

clearvars;
close all;
clc;
format long g

tic

plotRequested = 1;

addpath(genpath('lib'));

%Reading the observation and navigation data into the program

obsfilename = "algo1180.21o";

Nav = Nav_reader("brdm1180.21p");
Obs = obs_read_rinex211(obsfilename);  %Pixel4_GnssLog.21o %39ea118x.21o

%dateOfObs = dateOfFirstObs(obsfilename);

% GAST at 0h of the observation date
% gast0h = utc2gast(dateOfObs(1),dateOfObs(2),dateOfObs(3),0,0,0);

% GMST at 0h of observation date
obsdatevec = Obs.date_obs;
obs = Obs.epochdata;
if isempty(Obs.date_obs)
else
gmst0 = gmstGLO(obsdatevec(1),obsdatevec(2),obsdatevec(3),0,0,0);
end
% Julian day number of observation
jd0 = greg2jd(obsdatevec(1),obsdatevec(2),obsdatevec(3),0,0,0);

%% Match observation and navigation data
gps_observations = Obs.epochdata(:,1:2);
gloobs = Obs.epochdata(:,[1 3]);
%obssize = size(Observations{1,2}); % Determining the number of observations
n_gpsrecords = size(Nav.data,1); % The total number of GPS ephemeris data records


%Iterating through the observations vector
for i = 1:length(gps_observations)
    obssize = size(gps_observations{i,2}); %Entering the current submatrix for the current observation epoch
    
    for j = 1:obssize(1) %Iterating through the number of satellites at each epoch
        
        satnum = gps_observations{i,2}(j,2);    %The current satellite number at the current epoch
        
        index = ephMatch(satnum,gps_observations{i,1},Nav.data(:,1),Nav.data(:,2));
        
        %         timediff = [];
        %         I = [];
        %         %Iterating through the navigation data satellite at a particular
        %         %epoch observations
        %         for k = 1:n_gpsrecords
        %
        %             %If the satellite number in the navigation data is the same as
        %             %that in the current satellite in the observation data
        %             if Nav.data(k,1) == satnum
        %
        %                 %Saving the index at which the match was made
        %                 I = [I,k];
        %
        %                 %Taking the difference between the navigation time block and the observation time for the same satellite
        %                 timediff = [timediff,(Observations{i,1} - Nav.data(k,2))];
        %
        %             end
        %         end
        %         minInd = ephMatch(satnum,Observations{i,1},Nav.data(:,1),Nav.data(:,2));
        %         [M,P] = min(abs(timediff)); %Determining the minimum difference in time between the navigation and observation time epochs for a particular satellite for matching
        %         testm(j,1) = M;
        %
        %         index = I(P); %Saving the minimum matched index for current reference
        %         minInd == index;
        
        % super secret test case
        %         if minInd ~= index
        %             error('aha');
        %         end
        gps_observations{i,2}(j,obssize(2)+1) = index; %Concatenating the best satellite epoch and their corresponding satellite information to the observation epoch it was matched to in a submatrix
    end
end

%% Loop through each satellite
f1 = 1575.42e6; %L1
f2 = 1227.6e6; %L2
gravitationalparameter = Const.MU; %WGS84 value of the Earth's gravitational constant for GPS
earthrotation = Const.OMEGAE; %WGS84 value of the Earth's rotation rate
c = Const.CMPS; %Speed of light
xapprox=0;yapprox=0;zapprox=0;
num_epoch = length(gps_observations);

% Reserve index:
%   i   for iterating through each epoch of observation
%   j   for iterating through every observation in the epoch of a
%       particular constellation

for i = 1:num_epoch %Iterating through every epoch
    obssize = size(gps_observations{i,2}); %Determining the size of the submatrix attached to the epoch of observations
    clear elevationangle azimuth dTsv
    remove = [];
    if(i == 1749)
        break;
    else
        
        for j = 1:obssize(1) %Iterating through the GPS submatrix attached to the observation epoch
            %P2 = Observations{i,2}(j,5); %P1 code
            %P1 = Observations{i,2}(j,6); %P2 code
            Pseudorange(i,j) = gps_observations{i,2}(j,3); %((f1^2)*P1 - (f2^2)*P2)/(f1^2-f2^2); %Calculating the uncorrected pseudorange
            index = gps_observations{i,2}(j,11);
            sqrtA = Nav.data(index, 13); %Finding the semi-major axis
            A = (sqrtA)^2; %Finding the semi-major axis
            N0 = sqrt(gravitationalparameter/A^3); %Computing the mean motion
            Toe = Nav.data(index, 14); %Reading in the ephemeris data reference time of week
            Toc = Nav.data(index, 24)*60*60*24*7 + Nav.data(index, 14); %Finding the clock data reference time of week
            T = Pseudorange(i,j)/c; %Calculating the travel time from the uncorrected pseudorange and the speed of light
            %%Not sure if T is applicable and Trec along with time
            time = cell2mat(gps_observations(i,1));
            Trec = time - T;
            Tk = Trec - Toc; %Calculating the time of ephemeris reference time
            %Checking that the time of ephemeris is between -302400 and 302400
            if Tk > 302400
                Tk = Tk - 604800;
            end
            if Tk < -302400
                Tk = Tk + 604800;
            end
            N = N0 + Nav.data(index, 8); %The corrected mean motion
            Mk = Nav.data(index, 9) + N*Tk; %Calculating the mean anomaly
            ecc = Nav.data(index, 11); %Reading in the eccentricity
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
            Toc = Nav.data(index, 24)*60*60*24*7 + Nav.data(index, 14); %Determining the clock data reference time of week
            dTr = (-4.442807633e-10)*ecc*sqrtA*sin(Ek); %Relativistic clock correction term
            dTsv(j) = Nav.data(index, 3) + Nav.data(index, 4)*(time-Toc) + Nav.data(index, 5)*((time-Toc)^2) + dTr; %Determining the satellite clock offset
            Vk(i,j) = atan2((sqrt(1-ecc^2)*sin(Ek)/(1-ecc*cos(Ek))),((cos(Ek)-ecc)/(1-ecc*cos(Ek)))); %Determining the true anomaly
            ArgLat(1) = Vk(i,j) + Nav.data(index, 20); %Finding the argument of latitude
            %The second harmonic perturbations
            ArgLat(2) = Nav.data(index, 12)*sin(2*ArgLat(1)) + Nav.data(index, 10)*cos(2*ArgLat(1)); %The argument of latitude correction
            ArgLat(3) = Nav.data(index, 7)*sin(2*ArgLat(1)) + Nav.data(index, 19)*cos(2*ArgLat(1)); %The radial correction
            ArgLat(4) = Nav.data(index, 17)*sin(2*ArgLat(1)) + Nav.data(index, 15)*cos(2*ArgLat(1)); %The inclination correction
            Corr(1) = ArgLat(1) + ArgLat(2); %The corrected argument of latitude
            Corr(2) = A*(1-ecc*cos(Ek))+ArgLat(3); %The corrected radius
            Corr(3) = Nav.data(index, 18) + Nav.data(index, 22)*Tk + ArgLat(4); %The corrected inclination angle
            %Determining the positions in the orbital plane
            Xkprime = Corr(2)*cos(Corr(1));
            Ykprime = Corr(2)*sin(Corr(1));
            
            OMEGAdot = Nav.data(index, 21); %Reading in the rate of right ascension
            OMEGAk = Nav.data(index, 16) + (OMEGAdot - earthrotation)*Tk - earthrotation*Toe; %Determining the corrected longitude of ascending node
            xk = Xkprime*cos(OMEGAk)-Ykprime*cos(Corr(3))*sin(OMEGAk); %The satellite's X component in ECEF coordinates
            yk = Xkprime*sin(OMEGAk)+ Ykprime*cos(Corr(3))*cos(OMEGAk); %The satellite's Y component in ECEF coordinates
            Zk = Ykprime*sin(Corr(3)); %The satellite's Z component in ECEF coordinates
            
            %Making corrections for the Earth's rotation
            theta = earthrotation * T;
            Xk = xk*cos(theta) + yk*sin(theta); %The X component in ECEF coordinates
            Yk = -xk*sin(theta) + yk*cos(theta); %The Y component in ECEF coordinates
            
            VectorXYZ = [Xk, Yk, Zk]; %The ECEF coordinates
            Approx = [Obs.pos_xyz(1,1) ; Obs.pos_xyz(1,2) ; Obs.pos_xyz(1,3)]; %The approximate position of the receiver
            Approxcoord = Approx';
            LLHtrans = EC2LLH(VectorXYZ); %Calculating the latitude, longitude, and height of the ECEF coordinates
            LLHtrans(1:2) = LLHtrans(1:2)*pi/180; %Converting the result to radians
            LLH = LLHtrans';
            ENU = EC2ENU(VectorXYZ, Approxcoord, LLH); %Determining the easting, northing, and up from the ECEF coordinates
            
            azimuth(j) = atan2(ENU(1,1),ENU(1,2)); %Determining the azimuth
            slantd = sqrt(ENU(1)^2 + ENU(2)^2 + ENU(3)^2); %Finding the slant distance
            zenith = acos(ENU(3)/slantd); %Finding the zenith angle
            elevationangle(j) = pi/2 - zenith; %Determining the elevation angle from the zenith angle
            if(elevationangle(j) < (10/180*pi)) %If the elevation angle is less tha 10 degrees, store the index
                remove = [remove,j];
            end
            
            gps_observations{i,2}(j,18) = elevationangle(j); %Storing the elevation angle
            gps_observations{i,2}(j,19) = azimuth(j); %Storing the zenith angle
            
            
            %Constants to determine the tropospheric error
            Pwv = 16;
            pd = 1013 - Pwv;
            Ttrop = 300;
            %Calculating and storing the tropospheric error
            dtrop = (0.002277/(cos(zenith)))*(pd + (1255/Ttrop + 0.05)*Pwv - 1.16*((tan(zenith))));
            gps_observations{i,2}(j,20) = dtrop;
            
            %The corrected pseudorange
            correctedpseudorange(i,j) = Pseudorange(i,j)- dtrop;
            
            %Storing the corrected pseudorange, satellite clock offest, and
            %ECEF Xk, Yk, and Zk values in the overall Observations matrix for
            %each satellite at each epoch
            gps_observations{i,2}(j,12) = correctedpseudorange(i,j);
            gps_observations{i,2}(j,13) = dTsv(j);
            gps_observations{i,2}(j,14) = Xk;
            gps_observations{i,2}(j,15) = Yk;
            gps_observations{i,2}(j,16) = Zk;
            
        end
        
        %Removing all satellites that are below 10 degrees
        for k = length(remove):-1:1
            gps_observations{i,2}(remove(k),:) = [];
        end
        
        %     Removing any rows with NaN in them from the observation file (rows
        %     that are missing GNSS code)
        %    check = isnan(Observations{i,2}(:,12));
        %     for j = obssize:-1:1
        %         if check(j) == 1
        %             Observations{i,2}(j,:) = [];
        %         end
        %     end
        
        check = isnan(gps_observations{i,2});
        for j = obssize:-1:1
            if check(j) == 1
                gps_observations{i,2} = [];
            end
        end
        
        obssize = size(gps_observations{i,2});
        
        num_gloobs = size(gloobs{i,2},1);
        
        % Processing GLONASS observations at epoch i
        
        trec = gloobs{i,1}; % receiver time in seconds since GPS epoch
        
        for j = 1:num_gloobs
                        
            % find navigation parameters closest in time with receiver time
            iNav = ephMatch(gloobs{i,2}(j,1),trec,...
                Nav.glonav.satnum,...
                Nav.glonav.gpst);
            % (Nav.glonav.gpst(indexnav) - trec)/60
            PR = gloobs{i,2}(j,2);
            
            % time of transmission
            ti = emissionTime(PR,trec,...
                Nav.glonav.gpst(iNav),...
                Nav.glonav.mTauN(iNav),...
                Nav.glonav.GammaN(iNav),...
                0);
            
            % emission time in seconds of the day
            ti_sod = ti - (jd0 - Const.DJGPS)*Const.DAYSEC;
            
            % Greenwich sidereal time at receiver time
            S = gmst0 + Const.OMEGAE*ti_sod;
            
            % transform ephemeris parameters from ECEF PZ90.11 system to inertial
            [xinert,yinert,zinert] = rotatez(Nav.glonav.xpos(iNav),...
                Nav.glonav.ypos(iNav),...
                Nav.glonav.zpos(iNav),...
                S);
            
            [vxinert,vyinert,vzinert] = rotatez(Nav.glonav.xvel(iNav),...
                Nav.glonav.yvel(iNav),...
                Nav.glonav.zvel(iNav),...
                S);
            
            vxinert = vxinert + Const.OMEGAE*xinert;
            vyinert = vyinert + Const.OMEGAE*yinert;
            
            [AxLS,AyLS,AzLS] = rotatez(Nav.glonav.xacc(iNav),...
                Nav.glonav.yacc(iNav),...
                Nav.glonav.zacc(iNav),...
                S);
            
            % Numerical integration of satellite motion via RK4
            [xi,yi,zi] = RK4GLO(Nav.glonav.gpst(iNav),ti,...
                xinert,yinert,zinert,...
                vxinert,vyinert,vzinert,...
                AxLS,AyLS,AzLS);
            
            % rotate back to PZ90.11
            [xi,yi,zi] = rotatez(xi,yi,zi,-S);
            
            % Transformation between PZ90.11 to ITRF2008
            % using cartesian transformation parameters from
            % https://eng.mil.ru/files/PZ-90.11_final-v8.pdf
            itrfcoord = cart2cart([xi; yi; zi], ...
                                  [-0.003; -0.001; 0],...
                                  [0.019; -0.042; 0.002]*Const.DMAS2R,...
                                  0);
                              
            gloobs{i,2}(j,6:8)=itrfcoord';
                              
        end
        
        
        
        %%%%%
        %% Least squares portion
        num_gpsobs = size(gps_observations{i,2}(:,12),1);
        L = [gps_observations{i,2}(:,12);gloobs{i,2}(:,2)]; %The current psuedorange measurement for the current satellite and epoch
        Xs = [gps_observations{i,2}(:,14);gloobs{i,2}(:,6)]; %The current x Earth-fixed coordinates of SV antenna phase center
        Ys = [gps_observations{i,2}(:,15);gloobs{i,2}(:,7)]; %The current y Earth-fixed coordinates of SV antenna phase center
        Zs = [gps_observations{i,2}(:,16);gloobs{i,2}(:,8)]; %The current z Earth-fixed coordinates of SV antenna phase center
        dTsv = [gps_observations{i,2}(:,13);zeros(num_gloobs,1)];
        
        [TrueXYZ,v,DOP]=pointPos(L,Xs,Ys,Zs,dTsv,xapprox,yapprox,zapprox);
        xapprox=TrueXYZ(1);yapprox=TrueXYZ(2);zapprox=TrueXYZ(3); % set previous epoch's solution as approximation for next epoch, does not change results it seems
        gps_observations{i,2}(:,17) = v(1:num_gpsobs);
        
        %%%%%
        DOP1 = DOP(1); %Finding the accuracy of the Xu component
        DOP2 = DOP(2); %Finding the accuracy of the Yu component
        DOP3 = DOP(3); %Finding the accuracy of the Zu component
        DOP4 = DOP(4); %Finding the accuracy of the time (t) component
        GDOP(i) = sqrt(DOP1 + DOP2 + DOP3 + DOP4); %Determining the GDOP
        PDOP(i) = sqrt(DOP1 + DOP2 + DOP3); %Finding the PDOP
        TDOP(i) = sqrt(DOP4); %Determining the TDOP
        HDOP(i) = sqrt(DOP1 + DOP2); %Finding the HDOP
        VDOP(i) = sqrt(DOP3); %Determining the VDOP
        
        %Converting the final X, Y, and Z components of the satellite position
        %into Easting, Northing, and Up
        %TrueXYZ = [Final(1:3)];
        TrueLLHtrans = EC2LLH(TrueXYZ);
        XYZTable = EC2LLH(TrueXYZ);
        Currentepoch = gps_observations{i,1};
        CombinedTable(i,:) = [XYZTable(1); XYZTable(2); XYZTable(3); Currentepoch];
        
        TrueLLHtrans(1:2) = TrueLLHtrans(1:2)*pi/180;
        TrueLLH = TrueLLHtrans';
        TrueENU(i, :) = EC2ENU(TrueXYZ', Approxcoord, TrueLLH);
        
    end
end

VT = table('Size',[size(CombinedTable)], 'VariableTypes', {'double', 'double', 'double', 'double'}, 'VariableNames',{'Lat','Long','Height','Epoch'});
VT{:,:}=CombinedTable;
writetable(VT,'LLH_T.csv');

%Holds all of azimuth, elevation, and tropospheric data and sorts them by satellite rather
%than by time for plotting purposes
for i = 1:1748     %2880
    obs = size(gps_observations{i,2});
    for j = 1:obs(1)
        sat = gps_observations{i,2}(j,2);
        gps_observations{sat,3}(i) = gps_observations{i,2}(j,18);
        gps_observations{sat,4}(i) = gps_observations{i,2}(j,19);
        gps_observations{sat,5}(i) = gps_observations{i,2}(j,20);
        gps_observations{sat,7}(i) = gps_observations{i,2}(j,17);
    end
    temp = size(gps_observations{i,2});
    gps_observations{i,6} = temp(1);
end

%The code here is purely for plotting purposes. Any alterations are just to
%represent the results more clearly and accurately.
%Steps through each of the satellites at each epoch
for i = 1:32  %32   %1748
    numobs = length(gps_observations{i,3});
    for j = numobs:-1:1 %Remove the zero values from the azimuth and elevations
        if gps_observations{i,3}(j) == 0
            gps_observations{i,3}(j) = [];
            gps_observations{i,4}(j) = [];
        end
    end
    diff = 2880-numobs;
    gps_observations{i,5} = [gps_observations{i,5},zeros(1,diff)]; %Filling out each satellite observation
    for j = 1:2880  %Inserting NaN where there would be zeros in the tropospheric error for plotting purposes
        if gps_observations{i,5}(j) == 0
            gps_observations{i,5}(j) = NaN;
        end
    end
    numobs = length(gps_observations{i,7});
    diff = 2880-numobs;
    gps_observations{i,7} = [gps_observations{i,7},zeros(1,diff)]; %Filling out each residual observation
    for j = 1:2880  %Inserting NaN where there would be zeros
        if gps_observations{i,7}(j) == 0
            gps_observations{i,7}(j) = NaN;
        end
    end
end

toc


if plotRequested
    %% Plots
    
    %Plotting the number of satellites over time
    satquan = cell2mat(gps_observations(:,6));
    f1 = figure;
    plot(satquan);
    title('The Number Of Satellites Versus Epoch');
    ylabel('Number Of Satellites');
    xlabel('Epoch Number');
    
    
    %Plotting the azimuth angle vs elevation angle
    f2 = figure;
    scatter(gps_observations{1,4},gps_observations{1,3},20,'.')
    hold on
    title('The Azimuth Angle Versus The Elevation Angle');
    ylabel('Elevation Angle (Radians)');
    xlabel('Azimuth Angle (Radians)');
    for i = 2:32  %1748   %32
        scatter(gps_observations{i,4},gps_observations{i,3},20,'filled')
    end
    hold off
    
    %Plotting the tropospheric error versus epoch
    f3 = figure;
    plot(gps_observations{1,5});
    hold on
    title('Tropospheric Error Per Epoch');
    ylabel('Error (Metres)');
    xlabel('Epoch Number');
    for i = 2:32  %1748  %32
        plot(gps_observations{i,5});
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
    plot(gps_observations{1,7});
    exportgraphics(gca,'v.png')
    hold on
    title('Residuals From The Least Squares Computations');
    ylabel('Residuals (Metres)');
    xlabel('Epoch Number');
    countn = 0;
    for i = 2:32   %32
        plot(gps_observations{i,7});
    end
    % title('Residuals From The Least Squares Computations');
    % ylabel('Residuals (Metres)');
    % xlabel('Epoch Number');
    hold off
    exportgraphics(f1,'refsatnum.png')
    exportgraphics(f2,'reflocalcoord.png')
    exportgraphics(f3,'reftroperr.png')
    exportgraphics(f4,'refenuerror.png')
    exportgraphics(f5,'refdop.png')
    exportgraphics(f6,'refv.png')
end





