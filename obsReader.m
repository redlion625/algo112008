%function obs = obsReader(filename)
% obsReader returns struct for every observation line in the observation
% file given by filename in the local directory 
fid=fopen("Pixel4XLModded_GnssLog.20o");  %Pixel4_GnssLog.21o
%fid=fopen(filename);
% reads header lines
while ~feof(fid)
    line=fgetl(fid);
    headerWidths=[60 20];
    epochHead=fixedWidth(line,headerWidths);
    if contains(epochHead(2),'END OF HEADER')
        break
    end
end
i=1;
count=1;
while ~feof(fid)
    % store epoch into every struct instance
    line=fgetl(fid);
    if line(1)=='>'
        epochWidths=[2 4 3 3 3 3 11 3 3];
        epochHead=fixedWidth(line,epochWidths);
%         year=epochHead(2);
%         month=epochHead(3);
%         day=epochHead(4);
%         hour=epochHead(5);
%         min=epochHead(6);
%         sec=epochHead(7);
%         epochFlag=epochHead(8);
%         numSat=epochHead(9);

        % populate struct instance for every observation record
        while ~feof(fid)
            obs(i).year=str2double(epochHead(2));
            obs(i).month=str2double(epochHead(3));
            obs(i).day=str2double(epochHead(4));
            obs(i).hour=str2double(epochHead(5));
            obs(i).min=str2double(epochHead(6));
            obs(i).sec=str2double(epochHead(7));
            obs(i).epochFlag=str2double(epochHead(8));
            obs(i).numSat=str2double(epochHead(9));
            
            line=fgetl(fid);
            obsWidths=[1 2 16 16 16 16 16 16 16 16];
            obsLine=fixedWidth(line,obsWidths);
            obs(i).constellation=obsLine(1);
            % strip() removes leading whitespace
            % replace() replaces the erroneous whitespace between numbers
            % that occurs randomly and assumes a value of 0
            obs(i).PRN=str2double(replace(strip(obsLine(2),'left')," ","0"));
            obs(i).C1C=str2double(replace(strip(obsLine(3),'left')," ","0"));
            obs(i).L1C=str2double(replace(strip(obsLine(4),'left')," ","0"));
            obs(i).D1C=str2double(replace(strip(obsLine(5),'left')," ","0"));
            obs(i).S1C=str2double(replace(strip(obsLine(6),'left')," ","0"));
            obs(i).C5X=str2double(replace(strip(obsLine(7),'left')," ","0"));
            obs(i).L5X=str2double(replace(strip(obsLine(8),'left')," ","0"));
            obs(i).D5X=str2double(replace(strip(obsLine(9),'left')," ","0"));
            obs(i).S5X=str2double(replace(strip(obsLine(10),'left')," ","0"));
            i=i+1;
            count=count+1;
%             if count>3
%                 break
%             end
        end
    end
end