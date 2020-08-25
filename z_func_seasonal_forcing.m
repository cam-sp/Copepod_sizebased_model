function [temp,Diff,zmld,Ls,dzdt]=func_seasonal_forcing(param, years)

%load file for the mixed layer depth
fash=load('mld_fasham_mat.mat');
%adapt it to the number of days for the run
zmld1=fash.alk3;
zmld1=zmld1(2:end);
zmld1(end)=zmld1(1);
zmld=[zmld1; repmat(zmld1(1:end),years-1,1)]; %mixed layer depth
dzdt=zmld(2:end)-zmld(1:end-1); %change in mixed layer depth
dzdt=[dzdt(1); dzdt]; %correct for vector size

%load temperature file
temp1=load('tempfile_lat55.txt')+3;
temp1=temp1(2:end);
temp=[temp1(1); repmat(temp1(1:end),years,1)];
temp=temp-273;

%mixing rate (input of N)
Diff=(param.Diff_min+max(0,dzdt))./zmld;

%params light function
clouds=1-(6/8);%cloud cover
kyear=0;
day_type=1;
tsol=1:1:365;
Lss=daily_insolation(kyear,55,tsol,day_type); %surface light
Ls=Lss.*clouds.*0.4.*4.6; %correct for cloud cover, convert to PAR (PAR/I=0.4) and convert to microeinsteins (1Wm^-2=4.6microeinteins)
%correct format
Ls=Ls';
Ls=[Ls(1); repmat(Ls(1:end),years,1)];
Ls=Ls';


end