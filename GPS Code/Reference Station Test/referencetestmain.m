%% ENG4000 GNSS Processing Main

% Column names of per epoch glonass array
% 1-5
% PRN,  Code,   Phase,  Doppler,    Signal Strength,
%6-8
% Xsat,Ysat,Zsat

wid = 'MATLAB:nearlySingularMatrix';
warning('off',wid);

clearvars;
close all;
clc;
format long g

tic

plotRequested = 1;

addpath(genpath('lib'));

%Reading the observation and navigation data into the program

obsfilename = "algo1180.21o";

Nav = glonavread("brdc1180.21g");
Obs = obs_read_rinex211(obsfilename);  %Pixel4_GnssLog.21o %39ea118x.21o

%dateOfObs = dateOfFirstObs(obsfilename);

% GAST at 0h of the observation date
% gast0h = utc2gast(dateOfObs(1),dateOfObs(2),dateOfObs(3),0,0,0);

% GMST at 0h of observation date
obsdatevec = Obs.date_obs;
obs = Obs.epochdata;
% if isempty(Obs.date_obs)
%     error('Calendar date of observations not found');
% else
%     gmst0 = gmstGLO(obsdatevec(1),obsdatevec(2),obsdatevec(3),0,0,0);
% end
% Julian day number of observation
jd0 = greg2jd(obsdatevec(1),obsdatevec(2),obsdatevec(3),0,0,0);

%% Match observation and navigation data
obs = Obs.epochdata;
gloobs = Obs.epochdata(:,[1 2]);
%obssize = size(Observations{1,2}); % Determining the number of observations
n_gpsrecords = size(Nav.data,1); % The total number of GPS ephemeris data records

ne = length(obs);
GDOP = zeros(ne,1);PDOP = zeros(ne,1);TDOP = zeros(ne,1);
HDOP = zeros(ne,1);VDOP = zeros(ne,1);
DX = zeros(ne,1);DY = zeros(ne,1);DZ = zeros(ne,1);

%Iterating through the observations vector
for i = 1:length(obs)
    obssize = size(obs{i,2}); %Entering the current submatrix for the current observation epoch
    
    for j = 1:obssize(1) %Iterating through the number of satellites at each epoch
        break
        satnum = obs{i,2}(j,2);    %The current satellite number at the current epoch
        
        index = ephMatch(satnum,obs{i,1},Nav.data(:,1),Nav.data(:,2));
        
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
        

        %         if minInd ~= index
        %             error('aha');
        %         end
        obs{i,2}(j,obssize(2)+1) = index; %Concatenating the best satellite epoch and their corresponding satellite information to the observation epoch it was matched to in a submatrix
    end
end

%% Loop through each satellite
f1 = 1575.42e6; %L1
f2 = 1227.6e6; %L2
gravitationalparameter = Const.MU; %WGS84 value of the Earth's gravitational constant for GPS
earthrotation = Const.OMEGAE; %WGS84 value of the Earth's rotation rate
c = Const.CMPS; %Speed of light
xapprox=0;yapprox=0;zapprox=0;
num_epoch = length(obs);
% Reserve index:
%   i   for iterating through each epoch of observation
%   j   for iterating through every observation in the epoch of a
%       particular constellation

for i = 1:num_epoch %Iterating through every epoch
    
    obssize = size(obs{i,2},1); %Determining the size of the submatrix attached to the epoch of observations
    clear elevationangle azimuth dTsv L
    remove = [];
    %if(i == 1)
    %   continue
    %else
    numglo = 0;
    XYZsat = zeros(1,3);
    for j = 1:obssize(1) %Iterating through the GPS submatrix attached to the observation epoch
        
        if obs{i,2}(j,1)==2
            numglo = numglo + 1;
            obssize = size(obs{i,2});
            
            %num_gloobs = size(gloobs{i,2},1);
            
            % Processing GLONASS observations at epoch i
            
            trec = obs{i,4}-Nav.LS; % receiver time in sod UTC
            
            %satnumobs = obs{i,2}(j,2);
            
            %for j = 1:num_gloobs
            
            
            % find navigation parameters closest in time with receiver time
            iNav = ephMatch(obs{i,2}(j,2),trec,...
                Nav.data.satnum,...
                Nav.data.sod);
            % (Nav.data.gpst(indexnav) - trec)/60
            PR = obs{i,2}(j,3)-Const.CMPS*(Nav.LS-Nav.mTauC);
            % time of transmission
            
            te = Nav.data.sod(iNav);
            
            ti=trec-PR/c;
            
            dtsv=ti-te;
            for s=1:2
                dtsv=dtsv+Nav.data.mTauN(iNav)-Nav.data.GammaN(iNav)*dtsv;
            end
            dtsv=-Nav.data.mTauN(iNav)+Nav.data.GammaN(iNav)*dtsv;
            
            ti=ti-dtsv;
            
            
            %ti=ti-Nav.data.mTauN(iNav)-Nav.data.GammaN(iNav)*(ti-te);%-Nav.mTauC;
            
            %h=ti-te;
            x=simpleRK4(ti-te,...
                [Nav.data.xpos(iNav);Nav.data.ypos(iNav);Nav.data.zpos(iNav);...
                Nav.data.xvel(iNav);Nav.data.yvel(iNav);Nav.data.zvel(iNav)],...
                [Nav.data.xacc(iNav);Nav.data.yacc(iNav);Nav.data.zacc(iNav)]);
            
            %                 ti = emissionTime(PR,trec,te,...
            %                     Nav.data.mTauN(iNav)+Nav.mTauC,...
            %                     Nav.data.GammaN(iNav),...
            %                     0);
            
            % emission time in seconds of the day
            %ti_sod = ti - (jd0 - Const.DJGPS)*Const.DAYSEC;
            %te_sod = te - (jd0 - Const.DJGPS)*Const.DAYSEC;
            %                 dt = ti_sod-te_sod;
            %                 % Greenwich sidereal time at receiver time
            %
            %                 Se = gmst0 + Const.OMEGAE*te;%Se=-Se;
            %
            %                 x0 = Nav.data.xpos(iNav);
            %                 y0 = Nav.data.ypos(iNav);
            %                 z0 = Nav.data.zpos(iNav);
            %                 vx0 = Nav.data.xvel(iNav);
            %                 vy0 = Nav.data.yvel(iNav);
            %                 vz0 = Nav.data.zvel(iNav);
            %                 Ax0 = Nav.data.xacc(iNav);
            %                 Ay0 = Nav.data.yacc(iNav);
            %                 Az0 = Nav.data.zacc(iNav);
            %
            %                 % transform ephemeris parameters from ECEF PZ90.11 system to inertial
            %                 [x0,y0,z0] = rotatez(x0,y0,z0,Se);
            %                 [vx0,vy0,vz0] = rotatez(vx0,vy0,vz0,Se);
            %                 vx0 = vx0-Const.OMEGAE*y0;
            %                 vy0 = vy0+Const.OMEGAE*x0;
            %                 [Ax0,Ay0,Az0] = rotatez(Ax0,Ay0,Az0,Se);
            %
            %                 % Numerical integration of satellite motion via RK4
            %                 [x,y,z] = RK4GLO(te_sod,ti_sod,x0,y0,z0,vx0,vy0,vz0,0,0,0);
            %                 dx = Ax0*dt^2/2;dy = Ay0*dt^2/2; dz = Az0*dt^2/2;
            %                 x=x+dx;y=y+dy;z=z+dz;
            %                 % rotate back to PZ90.11
            %                 Si = gmst0 + Const.OMEGAE*ti_sod;%Si=-Si;
            %                 [x,y,z] = rotatez(x,y,z,-Si);
            
            % Transformation between PZ90.11 to ITRF2008
            % using cartesian transformation parameters from
            % https://eng.mil.ru/files/PZ-90.11_final-v8.pdf
%             xyz = cart2cart([x; y; z], ...
%                 [-0.003; -0.001; 0],...
%                 [0;0;0],...
%                 ...%[0.019; -0.042; 0.002]*Const.DMAS2R,...
%                 0);
            %xyz = [x;y;z];
            
            obs{i,2}(j,11:13)=x(1:3)';
            XYZsat(numglo,:) = x(1:3)';
            L(numglo) = PR;
        end
        
        %end
    end
    %%%%%
    %% Least squares portion
    %num_gpsobs = size(obs{i,2}(:,12),1);
    %L = [obs{i,2}(:,12);gloobs{i,2}(:,2)]; %The current psuedorange measurement for the current satellite and epoch
    %Xs = [obs{i,2}(:,14);gloobs{i,2}(:,6)]; %The current x Earth-fixed coordinates of SV antenna phase center
    %Ys = [obs{i,2}(:,15);gloobs{i,2}(:,7)]; %The current y Earth-fixed coordinates of SV antenna phase center
    %Zs = [obs{i,2}(:,16);gloobs{i,2}(:,8)]; %The current z Earth-fixed coordinates of SV antenna phase center
    %dTsv = [obs{i,2}(:,13);zeros(num_gloobs,1)];
    Xs = XYZsat(:,1);Ys = XYZsat(:,2);Zs = XYZsat(:,3);
    nL = numel(L);
    dTsv = zeros(nL,1);
    
    %L=reshape(L,[nL 1]);
    L=L';
    
    [TrueXYZ,v,DOP]=pointPos(L,Xs,Ys,Zs,dTsv,xapprox,yapprox,zapprox);
    xapprox=TrueXYZ(1);yapprox=TrueXYZ(2);zapprox=TrueXYZ(3); % set previous epoch's solution as approximation for next epoch, does not change results it seems
    %obs{i,2}(:,17) = v(1:nL);
    
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
    %     xbar = 918129.40;%.3
    %     ybar = -4346071.2;%.281
    %     zbar = 4561977.880;
    xbar = 918129.3;
    ybar = -4346071.281;
    zbar = 4561977.880;
    DX(i) = TrueXYZ(1)-xbar;
    DY(i) = TrueXYZ(2)-ybar;
    DZ(i) = TrueXYZ(3)-zbar;
    
    %Converting the final X, Y, and Z components of the satellite position
    %into Easting, Northing, and Up
    %TrueXYZ = [Final(1:3)];
    TrueLLHtrans = EC2LLH(TrueXYZ);
    XYZTable = EC2LLH(TrueXYZ);
    Currentepoch = obs{i,1};
    CombinedTable(i,:) = [XYZTable(1); XYZTable(2); XYZTable(3); Currentepoch];
    
    TrueLLHtrans(1:2) = TrueLLHtrans(1:2)*pi/180;
    TrueLLH = TrueLLHtrans';
    %TrueENU(i, :) = EC2ENU(TrueXYZ', Approxcoord, TrueLLH);
    
end

VT = table('Size',[size(CombinedTable)], 'VariableTypes', {'double', 'double', 'double', 'double'}, 'VariableNames',{'Lat','Long','Height','Epoch'});
VT{:,:}=CombinedTable;
writetable(VT,'LLH_T.csv');

plot(DX,DY)

if 1
    
else
    %Holds all of azimuth, elevation, and tropospheric data and sorts them by satellite rather
    %than by time for plotting purposes
    for i = 1:1748     %2880
        obs = size(obs{i,2});
        for j = 1:obs(1)
            sat = obs{i,2}(j,2);
            obs{sat,3}(i) = obs{i,2}(j,18);
            obs{sat,4}(i) = obs{i,2}(j,19);
            obs{sat,5}(i) = obs{i,2}(j,20);
            obs{sat,7}(i) = obs{i,2}(j,17);
        end
        temp = size(obs{i,2});
        obs{i,6} = temp(1);
    end
    
    %The code here is purely for plotting purposes. Any alterations are just to
    %represent the results more clearly and accurately.
    %Steps through each of the satellites at each epoch
    for i = 1:32  %32   %1748
        numobs = length(obs{i,3});
        for j = numobs:-1:1 %Remove the zero values from the azimuth and elevations
            if obs{i,3}(j) == 0
                obs{i,3}(j) = [];
                obs{i,4}(j) = [];
            end
        end
        diff = 2880-numobs;
        obs{i,5} = [obs{i,5},zeros(1,diff)]; %Filling out each satellite observation
        for j = 1:2880  %Inserting NaN where there would be zeros in the tropospheric error for plotting purposes
            if obs{i,5}(j) == 0
                obs{i,5}(j) = NaN;
            end
        end
        numobs = length(obs{i,7});
        diff = 2880-numobs;
        obs{i,7} = [obs{i,7},zeros(1,diff)]; %Filling out each residual observation
        for j = 1:2880  %Inserting NaN where there would be zeros
            if obs{i,7}(j) == 0
                obs{i,7}(j) = NaN;
            end
        end
    end
    
    toc
    
    
    if plotRequested
        %% Plots
        
        %Plotting the number of satellites over time
        satquan = cell2mat(obs(:,6));
        f1 = figure;
        plot(satquan);
        title('The Number Of Satellites Versus Epoch');
        ylabel('Number Of Satellites');
        xlabel('Epoch Number');
        
        
        %Plotting the azimuth angle vs elevation angle
        f2 = figure;
        scatter(obs{1,4},obs{1,3},20,'.')
        hold on
        title('The Azimuth Angle Versus The Elevation Angle');
        ylabel('Elevation Angle (Radians)');
        xlabel('Azimuth Angle (Radians)');
        for i = 2:32  %1748   %32
            scatter(obs{i,4},obs{i,3},20,'filled')
        end
        hold off
        
        %Plotting the tropospheric error versus epoch
        f3 = figure;
        plot(obs{1,5});
        hold on
        title('Tropospheric Error Per Epoch');
        ylabel('Error (Metres)');
        xlabel('Epoch Number');
        for i = 2:32  %1748  %32
            plot(obs{i,5});
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
        plot(obs{1,7});
        exportgraphics(gca,'v.png')
        hold on
        title('Residuals From The Least Squares Computations');
        ylabel('Residuals (Metres)');
        xlabel('Epoch Number');
        countn = 0;
        for i = 2:32   %32
            plot(obs{i,7});
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
    
end



