% modify path first
% output:

clc; clear; tic;
Rv = 487.5;
Lv = 2.5e6;
varLevel1{1} = 'toa'; varLevel1{2} = 'sfc';
varLevel2{1} = 'all'; varLevel2{2} = 'clr';

%% different time series, 1mean 2000-03 to 2018-02(18*12). 2 mean 200207-201706(15*12)
for p_1 = 1:2
    [readme, tLin] = observeParameters(p_1); % fuction% readme,  tLin, tLin.time tLin.start tLin.inter tLin.startmonth 
    %modify path first
    kernelAccum_path = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/kernelsCal/kernels');%kernel's path
    inputpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/anomaly/');
    inputpath1 = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/rawdata_regrid/');
    outpath = fullfile('/data1/liuyincheng/Observe-process', tLin.time{p_1}, 'ERAi/radEffect/');
    auto_mkdir(outpath)

    %% Read in ERAi anomalies and temperature Climatology
    load([inputpath, 'dvars.mat']); % Clim_alb, Clim_q, Clim_ta, Clim_ts, dalb, dq, dta, dts, latf, lonf, plevf, readme, time
    load([inputpath1, 'vars.mat']); % alb, q, ta, ts, latf, lonf, plevf, readme, time

    nlatf = length(latf); nlonf = length(lonf);ntime=length(time);
    time2 = 1:ntime;

    %% define qAnom properly
    %because kernel q's unit is W/m2/K/100mb
    dhus2 = zeros(144, 73, 24, ntime);
    p_hus=tLin.startmonth{p_1}-2;% the first month cal hus 
    for monNum = 1:ntime
        dhus2(:, :, :, monNum) = (log(q(:, :, :, monNum)) - log(Clim_q(:, :, :, mod(monNum + p_hus, 12) + 1))) .* Clim_ta(:, :, :, mod(monNum + p_hus, 12) + 1).^2 * Rv / Lv;
    end

    dhus2(isnan(dhus2)) = 0;

    %% Read in Kernels
    % kernel data please connect McGill University Prof.Yi Huang Email: yi.huang@mcgill.ca
    %lwkernel is 144*73*24*12(lon*lat*level*month)
    % kernelAccum_path = '/home/lyc/repeat_Kernel_experiment/testdata/kernels_cal/kernel_accumulate_'; %kernel's path
    level = 24;
    t_scflevel = 25;
    startMonth = tLin.startmonth{p_1};

    for ii = 1:2
        var1 = char(varLevel1{ii}); % 1.toa, 2.sfc

        for jj = 1:2

            var2 = char(varLevel2{jj}); % 1.cld 2.clr

            wv_lwkernel = ncread([kernelAccum_path, '_', var1, '_', var2, '.nc'], 'wv_lwkernel');
            wv_swkernel = ncread([kernelAccum_path, '_', var1, '_', var2, '.nc'], 'wv_swkernel');
            t_lwkernel = ncread([kernelAccum_path, '_', var1, '_', var2, '.nc'], 't_lwkernel');
            ts_lwkernel = ncread([kernelAccum_path, '_', var1, '_', var2, '.nc'], 'ts_lwkernel');
            alb_swkernel = ncread([kernelAccum_path, '_', var1, '_', var2, '.nc'], 'alb_swkernel');

            if strcmp(var1, 'toa')
                %% Radiative Effect
                wvlwEffect = zeros(144, 73, level, ntime);
                wvswEffect = zeros(144, 73, level, ntime);
                taEffect = zeros(144, 73, level, ntime);
                tsEffect = zeros(144, 73, ntime);
                albEffect = zeros(144, 73, ntime);
                %need loop because Kernels are from bottom to top and nc data is from top to bottom also the mod function is like that because the data starts in March thus it goes in a cycle from 1-12 starting at 3
                startMonth = tLin.startmonth{p_1};
                for monNum = 1:ntime

                    for ll = 1:level
                        wvlwEffect(:, :, ll, monNum) = wv_lwkernel(:, :, 25 - ll, mod(monNum + startMonth -2, 12) + 1) .* dhus2(:, :, ll, monNum);
                        wvswEffect(:, :, ll, monNum) = wv_swkernel(:, :, 25 - ll, mod(monNum + startMonth -2, 12) + 1) .* dhus2(:, :, ll, monNum);
                        taEffect(:, :, ll, monNum) = t_lwkernel(:, :, 25 - ll, mod(monNum + startMonth -2, 12) + 1) .* dta(:, :, ll, monNum);
                    end

                    tsEffect(:, :, monNum) = ts_lwkernel(:, :, mod(monNum + startMonth -2, 12) + 1) .* dts(:, :, monNum);
                    albEffect(:, :, monNum) = alb_swkernel(:, :, mod(monNum + startMonth -2, 12) + 1) .* dalb(:, :, monNum) * 100;
                end

                %sum over all the pressure levels to get the total radiative effect for each individual month.
                wvlwEffect = squeeze(sum(wvlwEffect(:, :, :, :), 3));
                wvswEffect = squeeze(sum(wvswEffect(:, :, :, :), 3));
                taEffect = squeeze(sum(taEffect(:, :, :, :), 3));
            elseif strcmp(var1, 'sfc')
                %% Radiative Effect
                wvlwEffect = zeros(144, 73, level, ntime);
                wvswEffect = zeros(144, 73, level, ntime);
                taEffect = zeros(144, 73, t_scflevel, ntime);
                t0Effect = zeros(144, 73, ntime);
                tsEffect = zeros(144, 73, ntime);
                albEffect = zeros(144, 73, ntime);
                %need loop because Kernels are from bottom to top and ncl data is from top to bottom also the mod function is like that because the data starts in March thus it goes in a cycle from 1-12 starting at 3
                for monNum = 1:ntime

                    for ll = 1:level
                        wvlwEffect(:, :, ll, monNum) = wv_lwkernel(:, :, 25 - ll, mod(monNum + startMonth -2, 12) + 1) .* dhus2(:, :, ll, monNum);
                        wvswEffect(:, :, ll, monNum) = wv_swkernel(:, :, 25 - ll, mod(monNum + startMonth -2, 12) + 1) .* dhus2(:, :, ll, monNum);
                        taEffect(:, :, ll, monNum) = t_lwkernel(:, :, 26 - ll, mod(monNum + startMonth -2, 12) + 1) .* dta(:, :, ll, monNum);
                    end

                    t0Effect(:, :, monNum) = squeeze(t_lwkernel(:, :, 1, mod(monNum + startMonth -2, 12) + 1)) .* dts(:, :, monNum);
                    tsEffect(:, :, monNum) = ts_lwkernel(:, :, mod(monNum + startMonth -2, 12) + 1) .* dts(:, :, monNum);
                    albEffect(:, :, monNum) = alb_swkernel(:, :, mod(monNum + startMonth -2, 12) + 1) .* dalb(:, :, monNum) * 100;
                end

                %sum over cld the pressure levels to get the total radiative effect for each individual month.
                wvlwEffect = squeeze(sum(wvlwEffect(:, :, :, :), 3));
                wvswEffect = squeeze(sum(wvswEffect(:, :, :, :), 3));
                taEffect = squeeze(sum(taEffect(:, :, :, :), 3));
                taEffect = taEffect + t0Effect;
            end

            %% Regrid the latitude
            latPlot = 88.75:-2.5:-88.75; nlatPlot = length(latPlot);
            lonPlot = lonf; nlonPlot = nlonf;
            [Xlonf, Ylatf, Ttimef] = meshgrid(latf, lonf, time2);
            [Xlonff, Ylatff, Ttimeff] = meshgrid(latPlot, lonf, time2);

            wvlwEffect = interp3(Xlonf, Ylatf, Ttimef, wvlwEffect, Xlonff, Ylatff, Ttimeff);
            wvswEffect = interp3(Xlonf, Ylatf, Ttimef, wvswEffect, Xlonff, Ylatff, Ttimeff);
            taEffect = interp3(Xlonf, Ylatf, Ttimef, taEffect, Xlonff, Ylatff, Ttimeff);
            tsEffect = interp3(Xlonf, Ylatf, Ttimef, tsEffect, Xlonff, Ylatff, Ttimeff);
            albEffect = interp3(Xlonf, Ylatf, Ttimef, albEffect, Xlonff, Ylatff, Ttimeff);

            totalEffect = wvlwEffect + wvswEffect + taEffect + tsEffect + albEffect;

            %% save as nc file
            saveto = strcat(outpath, 'dradEffect_', var1, '_', var2, '.nc');
            ncid = netcdf.create(saveto, 'NC_WRITE');

            %Define the dimensions
            dimidlon = netcdf.defDim(ncid, 'longitude', nlonPlot);
            dimidlat = netcdf.defDim(ncid, 'latitude', nlatPlot);
            dimidtime = netcdf.defDim(ncid, 'time', ntime);

            %Define IDs for the dimension variables
            longitude_ID = netcdf.defVar(ncid, 'longitude', 'double', dimidlon);
            latitude_ID = netcdf.defVar(ncid, 'latitude', 'double', dimidlat);
            time_ID = netcdf.defVar(ncid, 'time', 'double', dimidtime);

            %Define the main variable
            t_ID = netcdf.defVar(ncid, 'taRadEff', 'double', [dimidlon dimidlat dimidtime]);
            wvlw_ID = netcdf.defVar(ncid, 'wvlwRadEff', 'double', [dimidlon dimidlat dimidtime]);
            wvsw_ID = netcdf.defVar(ncid, 'wvswRadEff', 'double', [dimidlon dimidlat dimidtime]);
            ts_ID = netcdf.defVar(ncid, 'tsRadEff', 'double', [dimidlon dimidlat dimidtime]);
            alb_ID = netcdf.defVar(ncid, 'albRadEff', 'double', [dimidlon dimidlat dimidtime]);
            total_ID = netcdf.defVar(ncid, 'totalRadEff', 'double', [dimidlon dimidlat dimidtime]);

            %Define units for the dimension variables
            netcdf.putAtt(ncid, longitude_ID, 'standard_name', 'longitude');
            netcdf.putAtt(ncid, latitude_ID, 'standard_name', 'latitude');

            % netcdf.putAtt(ncid,pressure_ID,'standard_name','pressure');
            netcdf.putAtt(ncid, time_ID, 'long_name', 'time');
            netcdf.putAtt(ncid, time_ID, 'units', 'days since 0000-00-00 00:00:0.0');
            netcdf.putAtt(ncid, time_ID, 'calendar', 'gregorian');

            netcdf.putAtt(ncid, t_ID, 'long_name', 'Radiation Effect due to Temperature Anomaly (longwave)');
            netcdf.putAtt(ncid, t_ID, 'units', 'W m-2');

            netcdf.putAtt(ncid, wvlw_ID, 'long_name', 'Radiation Effect due to Specific Humidity Anomaly (longwave)');
            netcdf.putAtt(ncid, wvlw_ID, 'units', 'W m-2');

            netcdf.putAtt(ncid, wvsw_ID, 'long_name', 'Radiation Effect due to Specific Humidity Anomaly (shortwave');
            netcdf.putAtt(ncid, wvsw_ID, 'units', 'W m-2');

            netcdf.putAtt(ncid, ts_ID, 'long_name', 'Radiation Effect due to Surface Temperature Anomaly (longwave)');
            netcdf.putAtt(ncid, ts_ID, 'units', 'W m-2');

            netcdf.putAtt(ncid, alb_ID, 'long_name', 'Radiation Effect due to Albedo Anomaly (shortwave)');
            netcdf.putAtt(ncid, alb_ID, 'units', 'W m-2');

            netcdf.putAtt(ncid, total_ID, 'long_name', 'Total Radiation Effect');
            netcdf.putAtt(ncid, total_ID, 'units', 'W m-2');

            %We are done defining the NetCdf
            netcdf.endDef(ncid);

            %Then store the dimension variables in
            netcdf.putVar(ncid, longitude_ID, lonPlot);
            netcdf.putVar(ncid, latitude_ID, latPlot);
            netcdf.putVar(ncid, time_ID, time);

            %Then store my main variable

            netcdf.putVar(ncid, t_ID, taEffect);
            netcdf.putVar(ncid, wvlw_ID, wvlwEffect);
            netcdf.putVar(ncid, wvsw_ID, wvswEffect);
            netcdf.putVar(ncid, ts_ID, tsEffect);
            netcdf.putVar(ncid, alb_ID, albEffect);
            netcdf.putVar(ncid, total_ID, totalEffect);

            %We're done, close the netcdf
            netcdf.close(ncid)
        end

    end

end

t = toc;
disp(t)
