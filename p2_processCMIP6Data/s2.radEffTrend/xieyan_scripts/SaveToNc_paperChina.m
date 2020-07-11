%% datasets
clear

vector_var_str{1} = 'toa';
vector_var_str{2} = 'sfc';

vector_var_str2{1} = '';
vector_var_str2{2} = '_lw';
vector_var_str2{3} = '_sw';

vector_var_str3{1} = 'all';
vector_var_str3{2} = 'clr';

filename1 = '/data1/xieyan/sfc_all.nc';
% filename1 = '/data1/xieyan/EBAF_clear/sfc_all.nc';
time = ncread(filename1,'time'); ntime = length(time);
lon = ncread(filename1,'longitude'); nlon = length(lon);
lat = ncread(filename1,'latitude'); nlat = length(lat);
lonf = 70.25:0.5:140.25; nlonf = length(lonf);
latf = 55.25:-0.5:14.75; nlatf = length(latf);
time2 = 1:ntime; 
[Xlon,Ylat,Ttime] = meshgrid(lat,lon,time2);
[Xlonf,Ylatf,Ttimef] = meshgrid(latf,lonf,time2); 
%% EBAF & ERAi
% 2) dRx (time, lat, lon, budget, spectrum, sky, x): radiation anomaly, deseasoned and detrended
    % lat is above 70 degree north; 
    % budget: TOA, SFC;
    % spectrum: NET, LW, SW;
    % sky: all-sky, clear-sky;
    % x: total, temperature, water vapor, albedo, cloud, non-cloud,
    % diagnosed all, residual, ta(atmosphere), ts(surface)
    % datasets_radiation2018 / datasets_radiation2018_supplement
    dRx = NaN(5,2,3,2,nlonf,nlatf,216);                                        % deseasonalized and detrended
    dRx0 = NaN(5,2,3,2,nlonf,nlatf,216);                                        % only deseasonalized
    trend_Rx = NaN(5,2,3,2,nlonf,nlatf,12,2);
    p_Rx = NaN(5,2,3,2,nlonf,nlatf,12);
    cons_Rx = NaN(5,2,3,2,nlonf,nlatf,12);
    
    dRxs = NaN(5,2,3,2,nlonf,nlatf,216);                                        % deseasonalized and detrended
    dRxs0 = NaN(5,2,3,2,nlonf,nlatf,216);                                        % only deseasonalized
    trend_Rxs = NaN(5,2,3,2,nlonf,nlatf,12,2);
    p_Rxs = NaN(5,2,3,2,nlonf,nlatf,12);
    cons_Rxs = NaN(5,2,3,2,nlonf,nlatf,12);
    
    for ii = 1:2
        for zz = 1:2
            for jj = 1:3
                var = char(vector_var_str{ii});
                var3 = char(vector_var_str3{zz});
                var2 = char(vector_var_str2{jj});
                filename = strcat('/data1/xieyan/',var,'_',var3,'.nc');
%                 filename = strcat('/data1/xieyan/EBAF_clear/',var,'_',var3,'.nc');
                
                tot0 = ncread(filename,strcat('dR',var2)); tot0 = interp3(Xlon,Ylat,Ttime,tot0,Xlonf,Ylatf,Ttimef); tot = NaN(size(tot0));
                T0 = ncread(filename,'dR_T'); T0 = interp3(Xlon,Ylat,Ttime,T0,Xlonf,Ylatf,Ttimef); T = NaN(size(T0));
                ta0 = ncread(filename,'dR_ta'); ta0 = interp3(Xlon,Ylat,Ttime,ta0,Xlonf,Ylatf,Ttimef); ta = NaN(size(ta0));
                ts0 = ncread(filename,'dR_ts'); ts0 = interp3(Xlon,Ylat,Ttime,ts0,Xlonf,Ylatf,Ttimef); ts = NaN(size(ts0));
                q0 = ncread(filename,strcat('dR',var2,'_q')); q0 = interp3(Xlon,Ylat,Ttime,q0,Xlonf,Ylatf,Ttimef); q = NaN(size(q0));
                a0 = ncread(filename,strcat('dR_alb')); a0 = interp3(Xlon,Ylat,Ttime,a0,Xlonf,Ylatf,Ttimef); a = NaN(size(a0));
                nonc0 = ncread(filename,strcat('dR',var2,'_noncloud')); nonc0 = interp3(Xlon,Ylat,Ttime,nonc0,Xlonf,Ylatf,Ttimef); nonc = NaN(size(nonc0));
                res0 = ncread(filename,strcat('dR',var2,'_res')); res0 = interp3(Xlon,Ylat,Ttime,res0,Xlonf,Ylatf,Ttimef); res = NaN(size(res0));
                if strcmp(var3,'all')
                    c0 = ncread(filename,strcat('dR',var2,'_c')); c0 = interp3(Xlon,Ylat,Ttime,c0,Xlonf,Ylatf,Ttimef); c = NaN(size(c0));
                    diag0 = nonc0+c0; diag = NaN(size(diag0));
                elseif strcmp(var3,'clr')
                    diag0 = nonc0; diag = NaN(size(diag0));
                end
                
                % detrending
                trend_tot = zeros(nlonf,nlatf,12,2); trend_T = zeros(nlonf,nlatf,12,2); trend_q = zeros(nlonf,nlatf,12,2); trend_a = zeros(nlonf,nlatf,12,2); trend_c = zeros(nlonf,nlatf,12,2);
                trend_nonc = zeros(nlonf,nlatf,12,2); trend_diag = zeros(nlonf,nlatf,12,2); trend_res = zeros(nlonf,nlatf,12,2); trend_ta = zeros(nlonf,nlatf,12,2); trend_ts = zeros(nlonf,nlatf,12,2);
                
                p_tot = zeros(nlonf,nlatf,12); p_T = zeros(nlonf,nlatf,12); p_q = zeros(nlonf,nlatf,12); p_a = zeros(nlonf,nlatf,12); p_c = zeros(nlonf,nlatf,12);
                p_nonc = zeros(nlonf,nlatf,12); p_diag = zeros(nlonf,nlatf,12); p_res = zeros(nlonf,nlatf,12); p_ta = zeros(nlonf,nlatf,12); p_ts = zeros(nlonf,nlatf,12);
                
                cons_tot = zeros(nlonf,nlatf,12); cons_T = zeros(nlonf,nlatf,12); cons_q = zeros(nlonf,nlatf,12); cons_a = zeros(nlonf,nlatf,12); cons_c = zeros(nlonf,nlatf,12);
                cons_nonc = zeros(nlonf,nlatf,12); cons_diag = zeros(nlonf,nlatf,12); cons_res = zeros(nlonf,nlatf,12); cons_ta = zeros(nlonf,nlatf,12); cons_ts = zeros(nlonf,nlatf,12);

                for i = 1:nlonf
                    for j = 1:nlatf
                        [tot(i,j,:),trend_tot(i,j,:,:),cons_tot(i,j,:),p_tot(i,j,:)] = detrend_yan(tot0(i,j,:),time); 
                        [T(i,j,:),trend_T(i,j,:,:),cons_T(i,j,:),p_T(i,j,:)] = detrend_yan(T0(i,j,:),time);
                        [q(i,j,:),trend_q(i,j,:,:),cons_q(i,j,:),p_q(i,j,:)] = detrend_yan(q0(i,j,:),time);
                        [a(i,j,:),trend_a(i,j,:,:),cons_a(i,j,:),p_a(i,j,:)] = detrend_yan(a0(i,j,:),time);
                        [nonc(i,j,:),trend_nonc(i,j,:,:),cons_nonc(i,j,:),p_nonc(i,j,:)] = detrend_yan(nonc0(i,j,:),time); 
                        [diag(i,j,:),trend_diag(i,j,:,:),cons_diag(i,j,:),p_diag(i,j,:)] = detrend_yan(diag0(i,j,:),time);
                        [res(i,j,:),trend_res(i,j,:,:),cons_res(i,j,:),p_res(i,j,:)] = detrend_yan(res0(i,j,:),time);
                        [ta(i,j,:),trend_ta(i,j,:,:),cons_ta(i,j,:),p_ta(i,j,:)] = detrend_yan(ta0(i,j,:),time);
                        [ts(i,j,:),trend_ts(i,j,:,:),cons_ts(i,j,:),p_ts(i,j,:)] = detrend_yan(ts0(i,j,:),time);
                        if strcmp(var3,'all')
                            [c(i,j,:),trend_c(i,j,:,:),cons_c(i,j,:),p_c(i,j,:)] = detrend_yan(c0(i,j,:),time);
                        end
                    end
                end
                
                % save to dRx and trend_Rx
                dRx(1,zz,jj,ii,:,:,:) = tot; dRx0(1,zz,jj,ii,:,:,:) = tot0;trend_Rx(1,zz,jj,ii,:,:,:,:) = trend_tot; p_Rx(1,zz,jj,ii,:,:,:) = p_tot; cons_Rx(1,zz,jj,ii,:,:,:) = cons_tot;
                dRx(3,zz,jj,ii,:,:,:) = q; dRx0(3,zz,jj,ii,:,:,:) = q0;  trend_Rx(3,zz,jj,ii,:,:,:,:) = trend_q; p_Rx(3,zz,jj,ii,:,:,:) = p_q; cons_Rx(3,zz,jj,ii,:,:,:) = cons_q;
                dRxs(1,zz,jj,ii,:,:,:) = nonc; dRxs0(1,zz,jj,ii,:,:,:) = nonc0;  trend_Rxs(1,zz,jj,ii,:,:,:,:) = trend_nonc; p_Rxs(1,zz,jj,ii,:,:,:) = p_nonc; cons_Rxs(1,zz,jj,ii,:,:,:) = cons_nonc;
                dRxs(2,zz,jj,ii,:,:,:) = diag; dRxs0(2,zz,jj,ii,:,:,:) = diag0;  trend_Rxs(2,zz,jj,ii,:,:,:,:) = trend_diag; p_Rxs(2,zz,jj,ii,:,:,:) = p_diag; cons_Rxs(2,zz,jj,ii,:,:,:) = cons_diag;
                dRxs(3,zz,jj,ii,:,:,:) = res; dRxs0(3,zz,jj,ii,:,:,:) = res0;  trend_Rxs(3,zz,jj,ii,:,:,:,:) = trend_res; p_Rxs(3,zz,jj,ii,:,:,:) = p_res; cons_Rxs(3,zz,jj,ii,:,:,:) = cons_res;
                if zz == 1
                    dRx(5,zz,jj,ii,:,:,:) = c; dRx0(5,zz,jj,ii,:,:,:) = c0;trend_Rx(5,zz,jj,ii,:,:,:,:) = trend_c; p_Rx(5,zz,jj,ii,:,:,:) = p_c; cons_Rx(5,zz,jj,ii,:,:,:) = cons_c;
                end
                if jj ==1
                    dRx(2,zz,jj,ii,:,:,:) = T; dRx0(2,zz,jj,ii,:,:,:) = T0;trend_Rx(2,zz,jj,ii,:,:,:,:) = trend_T; p_Rx(2,zz,jj,ii,:,:,:) = p_T; cons_Rx(2,zz,jj,ii,:,:,:) = cons_T;
                    dRx(4,zz,jj,ii,:,:,:) = a; dRx0(4,zz,jj,ii,:,:,:) = a0; trend_Rx(4,zz,jj,ii,:,:,:,:) = trend_a; p_Rx(4,zz,jj,ii,:,:,:) = p_a; cons_Rx(4,zz,jj,ii,:,:,:) = cons_a;
                    dRxs(4,zz,jj,ii,:,:,:) = ta; dRxs0(4,zz,jj,ii,:,:,:) = ta0;trend_Rxs(4,zz,jj,ii,:,:,:,:) = trend_ta; p_Rxs(4,zz,jj,ii,:,:,:) = p_ta; cons_Rxs(4,zz,jj,ii,:,:,:) = cons_ta;
                    dRxs(5,zz,jj,ii,:,:,:) = ts; dRxs0(5,zz,jj,ii,:,:,:) = ts0;trend_Rxs(5,zz,jj,ii,:,:,:,:) = trend_ts; p_Rxs(5,zz,jj,ii,:,:,:) = p_ts; cons_Rxs(5,zz,jj,ii,:,:,:) = cons_ts;
                elseif jj == 2
                    dRx(2,zz,jj,ii,:,:,:) = T; dRx0(2,zz,jj,ii,:,:,:) = T0;trend_Rx(2,zz,jj,ii,:,:,:,:) = trend_T; p_Rx(2,zz,jj,ii,:,:,:) = p_T; cons_Rx(2,zz,jj,ii,:,:,:) = cons_T;
                    dRxs(4,zz,jj,ii,:,:,:) = ta; dRxs0(4,zz,jj,ii,:,:,:) = ta0;trend_Rxs(4,zz,jj,ii,:,:,:,:) = trend_ta; p_Rxs(4,zz,jj,ii,:,:,:) = p_ta; cons_Rxs(4,zz,jj,ii,:,:,:) = cons_ta;
                    dRxs(5,zz,jj,ii,:,:,:) = ts; dRxs0(5,zz,jj,ii,:,:,:) = ts0;trend_Rxs(5,zz,jj,ii,:,:,:,:) = trend_ts; p_Rxs(5,zz,jj,ii,:,:,:) = p_ts; cons_Rxs(5,zz,jj,ii,:,:,:) = cons_ts;
                elseif jj == 3
                    dRx(4,zz,jj,ii,:,:,:) = a; dRx0(4,zz,jj,ii,:,:,:) = a0; trend_Rx(4,zz,jj,ii,:,:,:,:) = trend_a; p_Rx(4,zz,jj,ii,:,:,:) = p_a; cons_Rx(4,zz,jj,ii,:,:,:) = cons_a;
                end
            end
        end
    end
    
    % 5) R (month, lat, lon, budget, spectrum, sky): R monthly climatology
% month = 1:12; 
% [Xlon,Ylat,Tmonth] = meshgrid(lat,lon,month);
% [Xlonf,Ylatf,Tomnthf] = meshgrid(latf,lonf,month);    
%    R = NaN(2,3,2,nlonf,nlatf,12);
%     for ii = 1:2
%         for jj = 1:3
%             for zz = 1:2
%                 var = char(vector_var_str{ii});
%                 var2 = char(vector_var_str2{jj});
%                 var3 = char(vector_var_str3{zz});
%                  filename = strcat('/data1/xieyan/',var,'_',var3,'.nc');
% %                 filename = strcat('/data1/xieyan/EBAF_clear/',var,'_',var3,'.nc');
%                 
%                 if strcmp(var2,'_lw') || strcmp(var2,'_sw')
%                     temp = ncread(filename,strcat('Clim_dR',var2));
%                     temp = interp3(Xlon,Ylat,Tmonth,temp,Xlonf,Ylatf,Tmonthf); 
%                 else
%                     temp = ncread(filename,'Clim_dR');
%                     temp = interp3(Xlon,Ylat,Tmonth,temp,Xlonf,Ylatf,Tmonthf); 
%                 end
%                 % save to Rx
%                 R(zz,jj,ii,:,:,:) = temp;
%             end
%         end
%     end
% China Land Map - Only trend
map1 = '/data1/xieyan/maps/worldcountry.shp';
map2 = '/data1/xieyan/maps/Chinaregion.shp';
readmap1 = shaperead(map1);
readmap2 = shaperead(map2);
[Xlonch,Ylatch]=meshgrid(lonf,latf);  
isin1 = inpolygon(Xlonch,Ylatch,readmap1(49,1).X,readmap1(49,1).Y);         % mainland
isin2 = inpolygon(Xlonch,Ylatch,readmap2(6,1).X,readmap2(6,1).Y);           % Taiwan Province
isin3 = inpolygon(Xlonch,Ylatch,readmap2(5,1).X,readmap2(5,1).Y);           % Xizang Province
isin12 = isin1|isin2; 
isinch = isin12|isin3;
isinch = isinch'; 

for xx = 1:5
   for sky = 1:2
      for spec = 1:3 
         for budg = 1:2
             for month = 1:12
                    for coef = 1:2
                        temp = squeeze(trend_Rx(xx,sky,spec,budg,:,:,month,coef));
                        temp(~isinch) = NaN; 
                        trend_Rx(xx,sky,spec,budg,:,:,month,coef) = temp;
                        temp = squeeze(trend_Rxs(xx,sky,spec,budg,:,:,month,coef));
                        temp(~isinch) = NaN; 
                        trend_Rxs(xx,sky,spec,budg,:,:,month,coef) = temp;
                    end
                    temp = squeeze(p_Rx(xx,sky,spec,budg,:,:,month));
                    temp(~isinch) = NaN; 
                    p_Rx(xx,sky,spec,budg,:,:,month) = temp;
                    temp = squeeze(p_Rxs(xx,sky,spec,budg,:,:,month));
                    temp(~isinch) = NaN; 
                    p_Rxs(xx,sky,spec,budg,:,:,month) = temp;
              end
         end
      end
   end
end
month = 1:12;
two = 1:2;
three = 1:3;
five = 1:5;

%% 2:Radiation Detrended
% ncid = netcdf.create('/data1/xieyan/EBAF_clear/Trendrad_China.nc','NC_WRITE');
ncid = netcdf.create('/data1/xieyan/Trendrad_China.nc','NC_WRITE');
%Define the dimensions
dimidlon = netcdf.defDim(ncid,'longitude',nlonf);
dimidlat = netcdf.defDim(ncid,'latitude',nlatf);
dimidmonth = netcdf.defDim(ncid,'month',12);
dimidbudget = netcdf.defDim(ncid,'budget',2);
dimidspec = netcdf.defDim(ncid,'spectrum',3);
dimidsky = netcdf.defDim(ncid,'sky',2);
dimidx = netcdf.defDim(ncid,'x',5);
dimidcoef = netcdf.defDim(ncid,'coef',2);
%Define IDs for the dimension variables
longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
month_ID=netcdf.defVar(ncid,'month','double',dimidmonth);
budget_ID =netcdf.defVar(ncid,'budget','int',dimidbudget);
spec_ID =netcdf.defVar(ncid,'spectrum','int',dimidspec);
sky_ID =netcdf.defVar(ncid,'sky','int',dimidsky);
x_ID =netcdf.defVar(ncid,'x','int',dimidx);
coef_ID = netcdf.defVar(ncid,'coef','int',dimidcoef);
%Define the main variable ()
trend_Rx_ID = netcdf.defVar(ncid,'trend_Rx','double',[dimidx dimidsky dimidspec dimidbudget dimidlon dimidlat dimidmonth dimidcoef]);
p_Rx_ID = netcdf.defVar(ncid,'p_Rx','double',[dimidx dimidsky dimidspec dimidbudget dimidlon dimidlat dimidmonth]);
% cons_Rx_ID = netcdf.defVar(ncid,'cons_Rx','double',[dimidx dimidsky dimidspec dimidbudget dimidlon dimidlat dimidmonth]);
%Define units for the dimension variables
netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
netcdf.putAtt(ncid,month_ID,'long_name','month');
netcdf.putAtt(ncid,month_ID,'units','Month (JFMAMJJASOND)');
netcdf.putAtt(ncid,trend_Rx_ID,'long_name','Coefficients to Determine Radiation Linear Trend Component: trend = a*time + b');
netcdf.putAtt(ncid,trend_Rx_ID,'units','a:W m-2 d-1; b:W m-2');
netcdf.putAtt(ncid,p_Rx_ID,'long_name','p value of trending');
netcdf.putAtt(ncid,p_Rx_ID,'units','1');
% netcdf.putAtt(ncid,cons_Rx_ID,'long_name','Mean value of deseasonalized Radiation');
% netcdf.putAtt(ncid,cons_Rx_ID,'units','W m-2');
netcdf.putAtt(ncid,budget_ID,'standard_name','budget');
netcdf.putAtt(ncid,budget_ID,'units','TOA; SFC');
netcdf.putAtt(ncid,spec_ID,'standard_name','spectrum');
netcdf.putAtt(ncid,spec_ID,'units','NET; LW; SW');
netcdf.putAtt(ncid,sky_ID,'standard_name','sky');
netcdf.putAtt(ncid,sky_ID,'units','all-sky; clear-sky');
netcdf.putAtt(ncid,x_ID,'standard_name','x');
netcdf.putAtt(ncid,x_ID,'units','total, temperature, water vapor, albedo, cloud');
netcdf.putAtt(ncid,coef_ID,'units','a; b');
%We are done defining the NetCdf
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid,longitude_ID,lonf);
netcdf.putVar(ncid,latitude_ID,latf);
netcdf.putVar(ncid,month_ID,month);
netcdf.putVar(ncid,budget_ID,two);
netcdf.putVar(ncid,spec_ID,three);
netcdf.putVar(ncid,sky_ID,two);
netcdf.putVar(ncid,x_ID,five);
netcdf.putVar(ncid,coef_ID,two);
%Then store main variable
netcdf.putVar(ncid,trend_Rx_ID,trend_Rx);
netcdf.putVar(ncid,p_Rx_ID,p_Rx);
% netcdf.putVar(ncid,cons_Rx_ID,cons_Rx);
%We're done, close the netcdf
netcdf.close(ncid)

%% Supplementary detrended radiation
% ncid = netcdf.create('/data1/xieyan/EBAF_clear/Trendrad_Chinasupp.nc','NC_WRITE');
ncid = netcdf.create('/data1/xieyan/Trendrad_Chinasupp.nc','NC_WRITE');
%Define the dimensions
dimidlon = netcdf.defDim(ncid,'longitude',nlonf);
dimidlat = netcdf.defDim(ncid,'latitude',nlatf);
dimidmonth = netcdf.defDim(ncid,'month',12);
dimidbudget = netcdf.defDim(ncid,'budget',2);
dimidspec = netcdf.defDim(ncid,'spectrum',3);
dimidsky = netcdf.defDim(ncid,'sky',2);
dimidx = netcdf.defDim(ncid,'x',5);
dimidcoef = netcdf.defDim(ncid,'coef',2);
%Define IDs for the dimension variables
longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
month_ID=netcdf.defVar(ncid,'month','double',dimidmonth);
budget_ID =netcdf.defVar(ncid,'budget','int',dimidbudget);
spec_ID =netcdf.defVar(ncid,'spectrum','int',dimidspec);
sky_ID =netcdf.defVar(ncid,'sky','int',dimidsky);
x_ID =netcdf.defVar(ncid,'x','int',dimidx);
coef_ID = netcdf.defVar(ncid,'coef','int',dimidcoef);
%Define the main variable ()
trend_Rxs_ID = netcdf.defVar(ncid,'trend_Rxs','double',[dimidx dimidsky dimidspec dimidbudget dimidlon dimidlat dimidmonth dimidcoef]);
p_Rxs_ID = netcdf.defVar(ncid,'p_Rxs','double',[dimidx dimidsky dimidspec dimidbudget dimidlon dimidlat dimidmonth]);
% cons_Rxs_ID = netcdf.defVar(ncid,'cons_Rxs','double',[dimidx dimidsky dimidspec dimidbudget dimidlon dimidlat dimidmonth]);
%Define units for the dimension variables
netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
netcdf.putAtt(ncid,month_ID,'long_name','month');
netcdf.putAtt(ncid,month_ID,'units','Month (JFMAMJJASOND)');
netcdf.putAtt(ncid,trend_Rxs_ID,'long_name','Coefficients to Determine Radiation Linear Trend Component: trend = a*time + b');
netcdf.putAtt(ncid,trend_Rxs_ID,'units','a:W m-2 d-1; b:W m-2');
netcdf.putAtt(ncid,p_Rxs_ID,'long_name','p value of trending');
netcdf.putAtt(ncid,p_Rxs_ID,'units','1');
% netcdf.putAtt(ncid,cons_Rxs_ID,'long_name','Mean value of deseasonalized Radiation');
% netcdf.putAtt(ncid,cons_Rxs_ID,'units','W m-2');
netcdf.putAtt(ncid,budget_ID,'standard_name','budget');
netcdf.putAtt(ncid,budget_ID,'units','TOA; SFC');
netcdf.putAtt(ncid,spec_ID,'standard_name','spectrum');
netcdf.putAtt(ncid,spec_ID,'units','NET; LW; SW');
netcdf.putAtt(ncid,sky_ID,'standard_name','sky');
netcdf.putAtt(ncid,sky_ID,'units','all-sky; clear-sky');
netcdf.putAtt(ncid,x_ID,'standard_name','x');
netcdf.putAtt(ncid,x_ID,'units','non-cloud, diagnosed all, residual, ta, ts');
%We are done defining the NetCdf
netcdf.endDef(ncid);
%Then store the dimension variables in
netcdf.putVar(ncid,longitude_ID,lonf);
netcdf.putVar(ncid,latitude_ID,latf);
netcdf.putVar(ncid,month_ID,month);
netcdf.putVar(ncid,budget_ID,two);
netcdf.putVar(ncid,spec_ID,three);
netcdf.putVar(ncid,sky_ID,two);
netcdf.putVar(ncid,x_ID,five);
netcdf.putVar(ncid,coef_ID,two);
%Then store main variable
netcdf.putVar(ncid,trend_Rxs_ID,trend_Rxs);
netcdf.putVar(ncid,p_Rxs_ID,p_Rxs);
% netcdf.putVar(ncid,cons_Rxs_ID,cons_Rxs);
%We're done, close the netcdf
netcdf.close(ncid)
