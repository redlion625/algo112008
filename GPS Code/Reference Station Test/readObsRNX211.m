function obs = readObsRNX211(filename)
%% reference observations follow different RINEX v 2.11
% follows data format from testgpsobs
%clear all
%filename='algo1180.21o';
fclose('all');

fid=fopen(filename);

%% read header
line=fgetl(fid);
%obs=struct([]);
while ~contains(line,'END OF HEADER')
    if contains(line,'APPROX POSITION XYZ')
        %temp=fixedWidth(line,[14 14 14]);
        obs.pos_xyz=str2double(fixedWidth(line,[14 14 14]));
        %obs.approxX=str2num(temp(1));
        %obs.approxY=str2num(temp(2));
        %obs.approxZ=str2num(temp(3));
    elseif contains(line,'LEAP SECONDS')
        obs.leapsec=str2num(fixedWidth(line,6));
    elseif contains(line,'TYPES OF OBS')
        temp=fixedWidth(line,[6 6 6 6 6 6 6 6 6]);
        obs.obsPerEpoch=str2num(temp(1));
        obs.obsname=strtrim(temp(2:end));
    elseif contains(line,'INTERVAL')
        obs.interval=str2num(line(1:10));
    elseif contains(line,'DELTA H/E/N')
        %temp=fixedWidth(line,[14 14 14]);
        obs.ant_delta_HEN = str2double(fixedWidth(line,[14 14 14]));
        %deltaH=str2num(temp(1));
        %deltaE=str2num(temp(2));
        %deltaN=str2num(temp(3));
    end
    line=fgetl(fid);
end
%% read observation records
numEpoch=86400/obs.interval;
data=cell(numEpoch,2);
epoch=1;
while ~feof(fid)
    line=fgetl(fid);
    if ~isspace(line(3))
        temp=fixedWidth(line,[3 3 3 3 3 11 3 3 36]);
        year=str2num(temp(1));
        month=str2num(temp(2));
        day=str2num(temp(3));
        hour=str2num(temp(4));
        minute=str2num(temp(5));
        second=str2num(temp(6));
        flag=str2num(temp(7));
        satnum=str2num(temp(8));
        
        data{epoch,1} = toGPST(year+2000,month,day,hour,minute,second);
        epoch_data=zeros(satnum,obs.obsPerEpoch+2);
        sat=strtrim(temp(9));
        constellation=zeros(satnum,1);
        prn=zeros(satnum,1);
        if satnum>12
            sat=fixedWidth(sat,3*ones(1,12));
            sat=char(sat');
            constellation(1:12)=sat(:,1);
            prn(1:12)=str2double(string(sat(:,2:3)));
            line=fgetl(fid);
            sat=fixedWidth(line,[32 3*ones(1,satnum-12)]);
            sat=char(sat(2:end)');
            constellation(13:end)=sat(:,1);
            prn(13:end,:)=str2double(string(sat(:,2:3)));
        else
            sat=fixedWidth(sat,3*ones(1,satnum));
            sat=char(sat');
            constellation(1:satnum)=sat(:,1);
            prn(1:satnum,:)=str2double(string(sat(:,2:3)));
        end
        epoch_data(:,1)=constellation;
        epoch_data(:,2)=prn;
        for i = 1:satnum
            line=fgetl(fid);
            epoch_data(i,3:7)=str2double(fixedWidth(line,[14 16 17 16 15]));
            line=fgetl(fid);
            epoch_data(i,8:10)=str2double(fixedWidth(line,[14 16 17]));
        end
        data{epoch,2}=epoch_data;
        epoch=epoch+1;
    else
        error("Unexpected data format, check input file for correct RINEX v2.11 format")
        break
    end
%     if epoch==3
%         break
%     end
end
obs.epochdata=data;
fclose(fid);