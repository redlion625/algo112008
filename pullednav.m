disp('------------------Begin reading navigation file---------------------');
RI = 1;         % receiver indicator defined here for now, later pull out from header.
C1=[]; L1=[]; D1=[]; S1=[]; P2=[]; L2=[]; D2=[];  S2=[];
%fidobs = fopen('brdm1180.21p');    % open observation file
fidobs = fopen('brdm1190.21p'); %brdm0840.21p  %brdm1560.20p  %brdm2260.20p   %brdm1180.21p  %brdm1560.20p %brdm1360.20p');    % open observation file
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
        
    elseif ~isempty(strfind(line,'LEAP SECONDS'))
            temp = strsplit(line);
            Nav.leap_sec = str2double(temp(2));
    end
end
line=fgets(fidobs);
i_epoch = 0;        % index for measurement data epoch
i_obs = 0;          % index for every observation (SV and epoch)
while line~=-1      % do until end of file is reach where fgets returns -1
    i_epoch = i_epoch +1;
    %if(line(1) == 'G')
    if(line(1) == 'G') %if(line(1) == 'G')
        data(i_epoch,1) = str2num(line(2:3)); %str2num(line(1:3));
        %elseif(line(1) == '
        %data(i_epoch,1) = (line(1:3)); %PRN

        %end
        line = fgets(fidobs);
    %if statement
    
    %elseif(line(1) ~= 'G')
        
    elseif(line(1) ~= 'G') % elseif(line(1) ~= 'G')
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