function Nav = Nav_reader(Navemerisfile);
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
fidobs = fopen('brdm1190.21p'); %brdm1560.20p  %brdm2260.20p   %brdm1180.21p  %brdm1560.20p %brdm1360.20p');    % open observation file
% obsdata_interval=10;    % late pull out from header
head_lines = 0;         % initilize
Nav= struct();
while 1  % Skip header and look for end of header text
    head_lines = head_lines+1;
    line = fgets(fidobs);
    end_found = findstr(line,'END OF HEADER'); % search to find end of header
    if ~isempty(end_found)
        break;
        
     elseif ~isempty(strfind(line,'IONOSPHERIC CORR'))
            temp = strsplit(line);
            GAL = 'GAL'; 
            GPSA = 'GPSA'; 
            GPSB = 'GPSB'; 
            checkpara = temp(1); 
            if(strcmpi(checkpara,GAL) == 1)
                Nav.GAL = str2double(temp(2:5));
                test = GAL;
            elseif(strcmpi(checkpara,GPSA) == 1)
                Nav.GPSA = str2double(temp(2:5));
                test1 = GPSA;
            elseif(strcmpi(checkpara,GPSB) == 1)
                Nav.GPSB = str2double(temp(2:5));
                test2 = GPSB;
            end
            L = temp;
        
    else if ~isempty(strfind(line,'LEAP SECONDS'))
            temp = strsplit(line);
            Nav.leap_sec = str2double(temp(2));
        elseif ~isempty(strfind(line,'ION ALPHA'))
            temp = strsplit(line);
            for j=1:4
                ion_alpha(j) = str2num([temp{j+1}]);
            end
            Nav.ion_alpha = ion_alpha;
        elseif ~isempty(strfind(line,'ION BETA'))
            temp = strsplit(line);
            for j=1:4
                ion_beta(j) = str2num([temp{j+1}]);
            end
            Nav.ion_beta = ion_beta;
        elseif ~isempty(strfind(line,'DELTA-UTC'))
            temp = strsplit(line);
            for j=1:3
                delta_utc(j) = str2num([temp{j+1}]);
            end
            Nav.delta_utc = delta_utc;
        end
    end
end

line=fgets(fidobs); % get the first line of the body
i_epoch = 0;        % index for measurement data epoch
i_obs = 0;          % index for every observation (SV and epoch)
while line~=-1      % do until end of file is reach where fgets returns -1
    i_epoch = i_epoch +1;
    if(line(1) == 'E') %if(line(1) == 'G')
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
    %if statement
        
    elseif(line(1) ~= 'E') % elseif(line(1) ~= 'G')
%         line = fgetl(fidobs);
%         %line = fgets(fidobs);
%         continue
%         %line = fgetl(fidobs);
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
test = line;
Nav.data = data;
disp('----------------Completed reading nav file-------------------');
fclose('all');
%%%%%%%%% end Nav_reader.m %%%%%%%%%

%Navigationtime = Nav.data.data(:,2);