function [t] = robotconv_Effort(fid) 

%ROBOTCONV  Convert binary robot C data to ascii
%  This function reads in data from the specified 
%  fileid and converts it to ascii and saves it to
%  a structure.  
%     --called by : robotdataread.m
%     --inputs    : fileid (fid)
%     --outputs   : data structure, (t)
%     --calls     : n/a
%  Created by Alaa Ahmed
%  Last modified  27-Sep-2021

%read in the header
binheader = fread(fid,300,'char');   
header = char(binheader');

%find number of headerlines
s=regexp(header,'logheadsize ');
n0begin=s+12;
number=1; j=1;

while number
    if header(n0begin+j)~='0'
        number=0;
    else
        j=j+1;
    end
end

nbegin=n0begin+j;
number=1; i=1;
while number
    if header(nbegin+i)=='s'
        number=0;
    else
        i=i+1;
    end
end

nend=nbegin+i-2;
charnumber=header(nbegin:nend);
hlen=str2num(charnumber);

%find number of columns
s2=regexp(header,'logcolumns ');
ncol=str2num(header(s2+11:s2+13));

ncolactual = 21;

colheaderstr = {'frame',...
                'adc2','adc3','adc4','adc5','adc6','adc7',...
                'time_ms','statenumber','avstatenumber','robotstatenumber',...
                'grasp','ft_xworld','ft_yworld','ft_zworld',...
                'mom_xworld','mom_yworld','mom_zworld',...
                'cursorx','cursory','homex'};

%read in data
status = fseek(fid,hlen,-1);
data=fread(fid,inf,'double','ieee-le');

if ncol >0
    nrows=length(data)/ncol;
    if nrows > 0 
        datam = NaN*zeros(nrows,ncol);
        for i=1:ncol
            datam(:,i)=data([i:ncol:length(data)]);
        end
    else
        datam = NaN*zeros(1,ncol);
        fprintf('Warning: File is empty (no rows).\n');
    end

    a=datam;
    framecal=a(:,1)-a(1,1);
    t.frame=framecal+1;
    t.time=t.frame/200;
    if ncolactual > 1
        for i=2:ncolactual
            eval(['t.' colheaderstr{i} '=a(:,i);'])
        end
    end    
else
    fprintf('Warning: File is empty (no columns).\n');
    t=NaN;
end

 
