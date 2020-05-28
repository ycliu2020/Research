clc; clear;
% lab1 real data
% 结论: detrend_yan 中月份对应关系检查通过, 趋势的算法未检验
%----------------------------------read data---------------------------------------------------
% constant
sigma = 5.67e-8; % Stefan-Boltzmann constant:Wm - 2K - 4

%modify path first
inputpath1 =  '/data1/liuyincheng/Observe-rawdata/ERAi/Rad_ori/';
inputpath2 =  '/data1/liuyincheng/Observe-rawdata/ERAi/Heat_Flux/Heat_Flux_Accumlate/';
inputpath3 =  '/data1/liuyincheng/Observe-rawdata/ERAi/Ts/';

%% read raw data. 228 mean raw date 19yrs time line
dR_downlw_sfc_cld0 = zeros(360, 181, 228);
dts_sfc0 = zeros(360, 181, 228);
dR_netsw_sfc_cld0 = zeros(360, 181, 228);
time0 = zeros(228, 1);

for ii = 0:18
    % inputpath1 = strcat('E:\xieyan\ERAi\Rad_ori\rad2004.nc');
    filename1 = strcat(inputpath1,  'rad', num2str((ii + 2000),  '%04i'),  '.nc');
    filename2 = strcat(inputpath2,  'hf', num2str((ii + 2000),  '%04i'),  '.nc');
    filename3 = strcat(inputpath3,  't', num2str((ii + 2000),  '%04i'),  '.nc');
    temp1 = ncread(filename2,  'strd'); % Surface thermal radiation downwards
    temp2 = ncread(filename3,  'skt'); % ts
    temp3 = ncread(filename1,  'ssr'); % surface_net_downward_shortwave_flux
    temp4 = ncread(filename1,  'time'); % hours since 1900 - 01 - 01 00:00:00.0
    % sum up the cumulative radiation of time step 0-12 and 12-24
    dR_downlw_sfc_cld0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp1(:, :, 1:2:23);
    dts_sfc0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp2(:, :, :);
    dR_netsw_sfc_cld0(:, :, ii * 12 + 1:(ii + 1) * 12) = temp3(:, :, 1:2:23);
    time0(ii * 12 + 1:(ii + 1) * 12, 1) = temp4(1:2:23);
end

% transfor mat time
temp.name = strcat(inputpath1,  'rad', num2str((4 + 2000),  '%04i'),  '.nc');
timeUnits = ncreadatt(temp.name,  'time',  'units');
% timeCalendar = ncreadatt(temp.name, 'time', 'calendar');
timeCalendar =  'gregorian';
time0 = cdfdate2num(timeUnits, timeCalendar, time0);

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
p_1 = 1;
[readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin,

% find start and end month
timeSpec = tLin.start{p_1};
formatOut =  'yyyy-mm';
timeStr = string(datestr(time0, formatOut));
startT = find(timeStr == timeSpec);
endT = startT + tLin.inter{p_1} - 1;

% read specific time series
dRdownlw_sfcAll = dR_downlw_sfc_cld0(:, :, startT:endT);
dts_sfc = dts_sfc0(:, :, startT:endT);
dRnetsw_sfcAll = dR_netsw_sfc_cld0(:, :, startT:endT);
time = time0(startT:endT, 1);
ntime = length(time);
% transform unite
dRdownlw_sfcAll = dRdownlw_sfcAll ./ (3600 * 24);
dRnetsw_sfcAll = dRnetsw_sfcAll ./ (3600 * 24);
%% Regrid to 2.5*2.5
lon = 0:359;
lat = 90:-1:-90;
lonPlot = 0:2.5:357.5; nlonf = length(lonPlot);
latPlot = 88.75:-2.5:-88.75; nlatf = length(latPlot);
time2 = 1:ntime;

[Xlon, Ylat, Ttime] = meshgrid(lat, lon, time2);
[Xlonf, Ylatf, Ttimef] = meshgrid(latPlot, lonPlot, time2);

dRdownlw_sfcAll = interp3(Xlon, Ylat, Ttime, dRdownlw_sfcAll, Xlonf, Ylatf, Ttimef);
dts_sfc = interp3(Xlon, Ylat, Ttime, dts_sfc, Xlonf, Ylatf, Ttimef);
dRnetsw_sfcAll = interp3(Xlon, Ylat, Ttime, dRnetsw_sfcAll, Xlonf, Ylatf, Ttimef);

% Deseasonalized
startmonth = tLin.startmonth{p_1};
[dRdownlw_sfcAll, ClimdRdownlw_sfcAll] = monthlyAnomaly3D(144, 72, time, dRdownlw_sfcAll, startmonth);
[dts_sfc, Climdts_sfc] = monthlyAnomaly3D(144, 72, time, dts_sfc, startmonth);
[dRnetsw_sfcAll, ClimdRnetsw_sfcAll] = monthlyAnomaly3D(144, 72, time, dRnetsw_sfcAll, startmonth);
month = 1:12;
% cal the deltaE(4*epis*T^3^delatT)
epis = 1; % emissivity
pnum = startmonth - 2;

for ii = 1:length(dts_sfc)
    dDeltaE(:, :, ii) = 4 * epis * sigma * squeeze(Climdts_sfc(:, :, mod(ii + pnum, 12) + 1)).^3 .* dts_sfc(:, :, ii);
end

%------------------------------cal trend test----------------------------------------------------------------
varb = squeeze(squeeze(dDeltaE(20, 20, :)));
% next step: use real time to compare which is right
% old method
ntime = length(time);
nyy = ntime / 12;
varbShape = reshape(varb, 12, nyy);
varbShapeTest = reshape(varb, nyy, 12);
timeShape = reshape(time, 12, nyy);
varbShape = varbShape'; % nyy*12
varbShapecopy = zeros(18, 12);
timeShape = timeShape'; % nyy*12
% trendMon = zeros(12,1);
trendMon = zeros(12, 2);
trendMon_startJan = zeros(12, 2);
pvalue = zeros(12, 1);
cons = zeros(12, 1);
x0 = ones(nyy, 1);

for i = 1:12
    cons(i, 1) = mean(varbShape(:, i));
    varbShape(:, i) = varbShape(:, i) - cons(i, 1); %  discrete point minus annual mean
    X = [x0, squeeze(timeShape(:, i))];
    [b, ~, r, ~, stats] = regress(varbShape(:, i), X);
    trendMon(i, 1) = b(2, 1); % k
    trendMon(i, 2) = b(1, 1); % b
    varbShape(:, i) = varbShape(:, i) - timeShape(:, i) .* trendMon(i, 1) - trendMon(i, 2); % detrend after dvars : discrete point minus (kx+b)
    varbShapecopy(:, i) = r; %
    pvalue(i, 1) = stats(1, 3);
end

varbShape = varbShape';
varbShapecopy = varbShapecopy';
varb1 = reshape(varbShape, ntime, 1);
ans1 = varbShapecopy - varbShape;
plot(time(1:18), varbShapecopy(2, :), time(1:18), varbShape(2, :))
% new method (precise query )
varb = dDeltaE;

period = time; % >= datenum(1988,03,01) & time <datenum(2017,03,1);
nperiod = length(time); % nperiod = length(time(period));
[yy, mm, dd] = datevec(time); % precise query
var_m = zeros(nlonf, nlatf, 12); % monthly averages
ntime = length(time);

% var = var(:,:,period);
for monthNum = 1:12
    month_index = mm == monthNum;
    var_m(:, :, monthNum) = nanmean(varb(:, :, month_index), 3); % this is the average of each month
    var_dividMon(:,:,:,monthNum)=varb(:, :, month_index);
    time_dividMon(:,monthNum)=time(month_index);
end

for i = 20:20 %1:nlon
    for j = 20:20 %1:nlat
        for monthNum = 1:12
            varbShape1=squeeze(squeeze(var_dividMon(i,j,:,:)));
            cons1(monthNum, 1) = mean(varbShape1(:, monthNum));
            varbShape1(:, monthNum) = varbShape1(:, monthNum) - cons1(monthNum, 1); %  discrete point minus annual mean
            X = [x0, squeeze(time_dividMon(:, monthNum))];
            [b, ~, r, ~, stats] = regress(varbShape1(:, monthNum), X);
            trendMon1(monthNum, 1) = b(2, 1); % k
            trendMon1(monthNum, 2) = b(1, 1); % b
            varbShape1(:, monthNum) = varbShape1(:, monthNum) - time_dividMon(:, monthNum) .* trendMon1(monthNum, 1) - trendMon1(monthNum, 2); % detrend after dvars : discrete point minus (kx+b)
            varbShapecopy1(:, monthNum) = r; %
            pvalue(monthNum, 1) = stats(1, 3);
        end        
    end
end





% % lab2: regular data
% time = linspace(1, 24, 24);
% ntime = length(time);

% varb = 20 * sin(time * pi / 6 + pi / 6);
