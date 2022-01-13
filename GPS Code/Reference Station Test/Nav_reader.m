function [Nav,navGLO] = Nav_reader(Navemerisfile);
% Nav_reader reads RINEX v 3.04, mixed navigation file
% Arrangement of Navigation data in the Nav struct
%Nav.data(:,1) = 'PRN',
%Nav.data(:,2) = 'GPST',
%Nav.data(:,3) = 'SV Clock Bias',
%Nav.data(:,4) = 'SV Clock Drift',
%Nav.data(:,5) = 'SV Clock Drift Rate',
%Nav.data(:,6) = 'IODE',
%Nav.data(:,7) = 'Crs',
%Nav.data(:,8) = 'Delta n',
%Nav.data(:,9) = 'Mo',
%Nav.data(:,10) = 'Cuc',
%Nav.data(:,11) = 'Eccentricity',
%Nav.data(:,12) = 'Cus',
%Nav.data(:,13) = 'Sqrt(a)',
%Nav.data(:,14) = 'TOE',
%Nav.data(:,15) = 'Cic',
%Nav.data(:,16) = 'OMEGA',
%Nav.data(:,17) = 'CIS',
%Nav.data(:,18) = 'Io',
%Nav.data(:,19) = 'Crc',
%Nav.data(:,20) = 'Omega',
%Nav.data(:,21) = 'OMEGA DOT',
%Nav.data(:,22) = 'IDOT',
%Nav.data(:,23) = 'L2 Channel Codes',
%Nav.data(:,24) = 'GPS Week',
%Nav.data(:,25) = 'L2 P Data Flag',
%Nav.data(:,26) = 'SV Accuracy',
%Nav.data(:,27) = 'SV Health',
%Nav.data(:,28) = 'TGD',
%Nav.data(:,29) = 'IODC',
%Nav.data(:,30) = 'Transmission Time',
%Nav.data(:,31) = 'Fit Interval'
disp('------------------Begin reading navigation file---------------------');
RI = 1;         % receiver indicator defined here for now, later pull out from header.
C1=[]; L1=[]; D1=[]; S1=[]; P2=[]; L2=[]; D2=[];  S2=[];
fidobs = fopen(Navemerisfile);    % open observation file
% obsdata_interval=10;    % late pull out from header
head_lines = 0;         % initilize
Nav=struct();
navGLO=struct();

%% Read header
countTimeCorr=0;
tSysCorr=struct();
countIonCorr=0;
ionCorr=struct();
line=fgets(fidobs);
while ~contains(line,'END OF HEADER')
    head_lines = head_lines+1;
    
    %end_found = findstr(line,'END OF HEADER'); % search to find end of header
    if contains(line,'TIME SYSTEM CORR')
        temp=fixedWidth(line,[4 18 16 7 5]);
        countTimeCorr=countTimeCorr+1;
        % where no time correction is mentioned, constellation is assumed
        % to conform to UTC
        tSysCorr(countTimeCorr).corrType=strip(temp(1));
        tSysCorr(countTimeCorr).a0=str2num(temp(2));
        tSysCorr(countTimeCorr).a1=str2num(temp(3));
        tSysCorr(countTimeCorr).Tref=str2num(temp(4));
        tSysCorr(countTimeCorr).refWeekNum=str2num(temp(5));

    elseif contains(line,'IONOSPHERIC CORR')
        temp=fixedWidth(line,[4 13 12 12 12]);
        countIonCorr=countIonCorr+1;
        ionCorr(countIonCorr).corrType=strip(temp(1));
        ionCorr(countIonCorr).params=str2double(temp(2:5));% see RINEX 304 for meanings of params
        Nav.ionCorr=ionCorr;
    elseif ~isempty(strfind(line,'LEAP SECONDS'))
        temp = strsplit(line);
        Nav.leap_sec = str2double(temp(2));
    elseif contains(line,'BRD_SUM G   -')
        Nav.numepochsGPS=str2num(line(50:58));
    elseif contains(line,'BRD_SUM R   -')
        Nav.numepochsGLO=str2num(line(50:58));
    elseif contains(line,'BRD_SUM E   -')
        Nav.numepochsGAL=str2num(line(50:58));
    elseif contains(line,'BRD_SUM C   -')
        Nav.numepochsBDS=str2num(line(50:58));
        %     elseif ~isempty(strfind(line,'ION ALPHA'))
        %         temp = strsplit(line);
        %         for j=1:4
        %             ion_alpha(j) = str2num([temp{j+1}]);
        %         end
        %         Nav.ion_alpha = ion_alpha;
        %     elseif ~isempty(strfind(line,'ION BETA'))
        %         temp = strsplit(line);
        %         for j=1:4
        %             ion_beta(j) = str2num([temp{j+1}]);
        %         end
        %         Nav.ion_beta = ion_beta;
        %     elseif ~isempty(strfind(line,'DELTA-UTC'))
        %         temp = strsplit(line);
        %         for j=1:3
        %             delta_utc(j) = str2num([temp{j+1}]);
        %         end
        %         Nav.delta_utc = delta_utc;
        %     end
    end
    line = fgets(fidobs);
end
Nav.tSysCorr=tSysCorr;
%% Store Navigation parameters for each constellation of every epoch
line=fgets(fidobs); % get the first line of the body
i_epoch = 0;        % index for measurement data epoch
iGLO=0;
i_obs = 0;          % index for every observation (SV and epoch)
while line~=-1      % do until end of file is reach where fgets returns -1
    i_epoch = i_epoch +1;
    if(line(1) == 'G')
        
        data(i_epoch,1) = str2num(line(2:3)); %str2num(line(1:3));
        %elseif(line(1) == '
        %data(i_epoch,1) = (line(1:3)); %PRN
        %data(i_epoch,1) = str2num(line(1:3)); %PRN
        year = str2num(line(5:8)) - 2000; %Epoch Year
        month = str2num(line(10:11)); %Epoch Month
        day = str2num(line(13:14)); %Epoch Day
        hour = str2num(line(16:17)); %Epoch Hour
        minute = str2num(line(19:20)); %Epoch Minute
        second = str2num(line(22:23)); %Epoch Second
        
        % convert the Epoch Year,Month, Day, Hour, Minute and Second to GPST
        data(i_epoch,2) = toGPST(year+2000,month,day,hour,minute,second);
        %data(i_epoch,2) = toGPST(year+2000,month,day,hour,minute,second);
        data(i_epoch,3) = str2num(line(24:42)); %SV Clock Bias
        data(i_epoch,4) = str2num(line(43:61)); %SV Clock Drift
        data(i_epoch,5) = str2num(line(63:end)); %SV Clock Drift Rate
        % remaining lines of data
        for i = 6:4:27
            line = fgetl(fidobs);
            data(i_epoch,i) = str2num(line(5:23));
            data(i_epoch,i+1) = str2num(line(24:42));
            data(i_epoch,i+2) = str2num(line(43:61));
            data(i_epoch,i+3) = str2num(line(62:end));
        end
        line = fgetl(fidobs);
        data(i_epoch,30) = str2num(line(5:23)); %Transmission Time
        if length(strfind(line,'.')) ~=1 %check if Fit Interval provided
            data(i_epoch,31) = str2num(line(25:42)); %Fit Interval
        end
        line = fgets(fidobs);

    elseif line(1)=='R'
        iGLO=iGLO+1;
        % give nav msg struct per epoch i header info

        
        % read line 1 of msg at epoch i
        temp=fixedWidth(line,[1 2 5 3 3 3 3 3 19 19 19]);
        navGLO(iGLO).satnum=str2num(temp(2));
        navGLO(iGLO).year=str2num(temp(3));
        navGLO(iGLO).month=str2num(temp(4));
        navGLO(iGLO).day=str2num(temp(5));
        %navGLO(iGLO).glonass2utc=glonass2utc;

        
        navGLO(iGLO).hr=str2num(temp(6));
        navGLO(iGLO).min=str2num(temp(7));
        navGLO(iGLO).sec=str2num(temp(8));
        navGLO(iGLO).clkbias=str2num(temp(9));
        navGLO(iGLO).freqbias=str2num(temp(10));
        navGLO(iGLO).msgframetime=str2num(temp(11));
        
        % read line 2
        line=fgetl(fidobs);
        temp=fixedWidth(line,[23 19 19 19]);
        navGLO(iGLO).xpos=str2num(temp(1))*1000; % km->m
        navGLO(iGLO).xvel=str2num(temp(2))*1000;
        navGLO(iGLO).xacc=str2num(temp(3))*1000;
        navGLO(iGLO).sathealth=str2num(temp(4));
        
        % read line 3
        line=fgetl(fidobs);
        temp=fixedWidth(line,[23 19 19 19]);
        navGLO(iGLO).ypos=str2num(temp(1))*1000;
        navGLO(iGLO).yvel=str2num(temp(2))*1000;
        navGLO(iGLO).yacc=str2num(temp(3))*1000;
        navGLO(iGLO).freqnum=str2num(temp(4));
        
        % read line 4
        line=fgetl(fidobs);
        temp=fixedWidth(line,[23 19 19 19]);
        navGLO(iGLO).zpos=str2num(temp(1))*1000;
        navGLO(iGLO).zvel=str2num(temp(2))*1000;
        navGLO(iGLO).zacc=str2num(temp(3))*1000;
        navGLO(iGLO).ageop=str2num(temp(4));
        line=fgets(fidobs);
    end
    if(line(1) ~= 'G'&&line(1)~='R')
        break
        
    end
end
%     if(line(1) == 'R')
%          data2(i_epoch,1) = 2; %str2num(line(1:3));
%      %elseif(line(1) == '
%      %data(i_epoch,1) = (line(1:3)); %PRN
%      %data(i_epoch,1) = str2num(line(1:3)); %PRN
%      year = str2num(line(5:8)) - 2000; %Epoch Year
%      month = str2num(line(10:11)); %Epoch Month
%      day = str2num(line(13:14)); %Epoch Day
%      hour = str2num(line(16:17)); %Epoch Hour
%      minute = str2num(line(19:20)); %Epoch Minute
%      second = str2num(line(22:23)); %Epoch Second
%
%      % convert the Epoch Year,Month, Day, Hour, Minute and Second to GPST
%      data2(i_epoch,2) = toGPST(year+2000,month,day,hour,minute,second);
%      %data(i_epoch,2) = toGPST(year+2000,month,day,hour,minute,second);
%      data2(i_epoch,3) = str2num(line(24:42)); %SV Clock Bias
%      data2(i_epoch,4) = str2num(line(43:61)); %SV Clock Drift
%      data2(i_epoch,5) = str2num(line(63:end)); %SV Clock Drift Rate
%      % remaining lines of data
%      for i = 6:4:27
%         line = fgetl(fidobs);
%         data2(i_epoch,i) = str2num(line(5:23));
%         data2(i_epoch,i+1) = str2num(line(24:42));
%         data2(i_epoch,i+2) = str2num(line(43:61));
%         data2(i_epoch,i+3) = str2num(line(62:end));
%       end
%       line = fgetl(fidobs);
%       data2(i_epoch,30) = str2num(line(5:23)); %Transmission Time
%       if length(strfind(line,'.')) ~=1 %check if Fit Interval provided
%           data2(i_epoch,31) = str2num(line(25:42)); %Fit Interval
%       end
%       line = fgets(fidobs);
%     end %if statement
%end %while line1~=-1
Nav.glonav=navGLO;
Nav.data = data;
disp('----------------Completed reading nav file-------------------');
fclose('all');
%%%%%%%%% end Nav_reader.m %%%%%%%%%

%Navigationtime = Nav.data.data(:,2);