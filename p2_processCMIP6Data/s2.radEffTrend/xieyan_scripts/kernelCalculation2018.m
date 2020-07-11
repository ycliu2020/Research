clear
close all
Rv = 487.5;
Lv = 2.5e6;
vector_var_str{1} = 'toa';
vector_var_str{2} = 'sfc';
vector_var_str2{1} = 'all';
vector_var_str2{2} = 'clr';

%% Read in ERAi anomalies and temperature Climatology
filename = '/data1/xieyan/ERAi/ERAi_variAnom2018.nc';
latf = ncread(filename,'latitude'); nlatf = length(latf);
lonf = ncread(filename,'longitude'); nlonf = length(lonf);
time2 = 1:216;
q = ncread(filename,'q');
time = ncread(filename,'time'); ntime = length(time);       
qAnom = ncread(filename,'Anom_q');
tAnom = ncread(filename,'Anom_t');
tsAnom = ncread(filename,'Anom_ts');
albAnom = ncread(filename,'Anom_alb');
%need the Climatology of t to properly calculate q with the kernels
tClim = ncread(filename,'Clim_t');
qClim = ncread(filename,'Clim_q');
%now define qAnom properly
qAnom2 = zeros(144,73,19,216);
for kk = 1:216
    qAnom2(:,:,:,kk) = (log(q(:,:,:,kk))-log(qClim(:,:,:,mod(kk+1,12)+1))).*tClim(:,:,:,mod(kk+1,12)+1).^2*Rv/Lv;
end
        
%% Read in Kernels
% kernel数据请联系加拿大McGill University Prof.Yi Huang Email: yi.huang@mcgill.ca
for ii = 1:2
    for jj = 1:2
        var = char(vector_var_str{ii});
        var2 = char(vector_var_str2{jj});
        if strcmp(var,'toa')
            if strcmp(var2,'all')
                wvlwCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_wv_toa_cld_nodp.nc','lwkernel');
                wvswCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_wv_toa_cld_nodp.nc','swkernel');
                tCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_t_toa_cld_nodp.nc','lwkernel');
                tsCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_ts_toa_cld_nodp.nc','lwkernel');
                albCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_alb_toa_cld_nodp.nc','swkernel');
            else
                wvlwClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_wv_toa_clr_nodp.nc','lwkernel');
                wvswClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_wv_toa_clr_nodp.nc','swkernel');
                tClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_t_toa_clr_nodp.nc','lwkernel');
                tsClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_ts_toa_clr_nodp.nc','lwkernel');
                albClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/toa/RRTMG_alb_toa_clr_nodp.nc','swkernel');
            end
        elseif strcmp(var,'sfc')
            %note that this has an extra dimension which holds upward-downward-net
            %fluxes (we should probs use net)
            if strcmp(var2,'all')
                wvlwCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_wv_surface_cld_nodp.nc','lwkernel');
                wvswCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_wv_surface_cld_nodp.nc','swkernel');
                tCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_t_surface_cld_nodp.nc','lwkernel');
                tsCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_ts_surface_cld_nodp.nc','lwkernel');
                albCldKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_alb_surface_cld_nodp.nc','swkernel');
                
                wvlwCldKernel = wvlwCldKernel(:,:,:,:,3);
                wvswCldKernel = wvswCldKernel(:,:,:,:,3);
                tCldKernel = tCldKernel(:,:,:,:,3);
                tsCldKernel = tsCldKernel(:,:,:,3);
                albCldKernel = albCldKernel(:,:,:,3); 
            else
                wvlwClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_wv_surface_clr_nodp.nc','lwkernel');
                wvswClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_wv_surface_clr_nodp.nc','swkernel');
                tClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_t_surface_clr_nodp.nc','lwkernel');
                tsClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_ts_surface_clr_nodp.nc','lwkernel');
                albClrKernel = ncread('/gs/project/ktg-565-aa/xxtan/data/kernel/rrtmgkernel/old/kernel-nodp/surface/RRTMG_alb_surface_clr_nodp.nc','swkernel');
                
                wvlwClrKernel = wvlwClrKernel(:,:,:,:,3);
                wvswClrKernel = wvswClrKernel(:,:,:,:,3);
                tClrKernel = tClrKernel(:,:,:,:,3);
                tsClrKernel = tsClrKernel(:,:,:,3);
                albClrKernel = albClrKernel(:,:,:,3);
            end
        end     
  
        %% Radiative Effect
        wvlwEffect = zeros(144,73,19,216);
        wvswEffect = zeros(144,73,19,216);
        tEffect = zeros(144,73,19,216);
        tsEffect = zeros(144,73,216);
        albEffect = zeros(144,73,216);
        %need loop because Kernels are from bottom to top and ncl data is from top
        %to bottom
        %also the mod function is like that because the data starts in March thus
        %it goes in a cycle from 1-12 starting at 3
        if strcmp(var2,'all')
            for kk = 1:216
                for ll = 1:19
                    wvlwEffect(:,:,ll,kk) = wvlwCldKernel(:,:,20-ll,mod(kk+1,12)+1).*qAnom2(:,:,ll,kk);
                    wvswEffect(:,:,ll,kk) = wvswCldKernel(:,:,20-ll,mod(kk+1,12)+1).*qAnom2(:,:,ll,kk);
                    tEffect(:,:,ll,kk) = tCldKernel(:,:,20-ll,mod(kk+1,12)+1).*tAnom(:,:,ll,kk);
                end
                tsEffect(:,:,kk) = tsCldKernel(:,:,mod(kk+1,12)+1).*tsAnom(:,:,kk);
                albEffect(:,:,kk) = albCldKernel(:,:,mod(kk+1,12)+1).*albAnom(:,:,kk)*100;
            end
        else
            for kk = 1:216
                for ll = 1:19
                    wvlwEffect(:,:,ll,kk) = wvlwClrKernel(:,:,20-ll,mod(kk+1,12)+1).*qAnom2(:,:,ll,kk);
                    wvswEffect(:,:,ll,kk) = wvswClrKernel(:,:,20-ll,mod(kk+1,12)+1).*qAnom2(:,:,ll,kk);
                    tEffect(:,:,ll,kk) = tClrKernel(:,:,20-ll,mod(kk+1,12)+1).*tAnom(:,:,ll,kk);
                end
                tsEffect(:,:,kk) = tsClrKernel(:,:,mod(kk+1,12)+1).*tsAnom(:,:,kk);
                albEffect(:,:,kk) = albClrKernel(:,:,mod(kk+1,12)+1).*albAnom(:,:,kk)*100;
            end
            
        end
        %%
        
        %sum over all the pressure levels to get the total radiative effect for
        %each individual month.
        wvlwEffect = sum(wvlwEffect(:,:,:,:),3);
        wvswEffect = sum(wvswEffect(:,:,:,:),3);
        tEffect = sum(tEffect(:,:,:,:),3);
        
        wvlwEffect = squeeze(wvlwEffect);
        wvswEffect = squeeze(wvswEffect);
        tEffect = squeeze(tEffect);
        
        %% Regrid the latitude
        latff = 88.75:-2.5:-88.75; nlatff = length(latff);
        [Xlonf, Ylatf, Ttimef] = meshgrid(latf,lonf,time2);
        [Xlonff, Ylatff, Ttimeff] = meshgrid(latff,lonf,time2);
        
        wvlwEffect = interp3(Xlonf,Ylatf,Ttimef,wvlwEffect,Xlonff,Ylatff,Ttimeff);
        wvswEffect = interp3(Xlonf,Ylatf,Ttimef,wvswEffect,Xlonff,Ylatff,Ttimeff);
        tEffect = interp3(Xlonf,Ylatf,Ttimef,tEffect,Xlonff,Ylatff,Ttimeff);
        tsEffect = interp3(Xlonf,Ylatf,Ttimef,tsEffect,Xlonff,Ylatff,Ttimeff);
        albEffect = interp3(Xlonf,Ylatf,Ttimef,albEffect,Xlonff,Ylatff,Ttimeff);
        
        totalEffect = wvlwEffect + wvswEffect + tEffect + tsEffect + albEffect;
        
        %% create a .nc file from the regridded data (reference: p-martineau.com/saving-netcdf-in-matlab/)
        directory_name = '/data1/xieyan/ERAi/';
%         saveto = strcat(directory_name,'eraiRadAnom_1988',upper(var),var2,'sky.nc');
        saveto = strcat(directory_name,'eraiRadAnom_2018',upper(var),var2,'sky.nc');
        ncid = netcdf.create(saveto,'NC_WRITE');
        
        %Define the dimensions
        dimidlon = netcdf.defDim(ncid,'longitude',nlonf);
        dimidlat = netcdf.defDim(ncid,'latitude',nlatff);
        dimidtime = netcdf.defDim(ncid,'time',216);
        
        %Define IDs for the dimension variables
        longitude_ID=netcdf.defVar(ncid,'longitude','double',dimidlon);
        latitude_ID=netcdf.defVar(ncid,'latitude','double',dimidlat);
        time_ID=netcdf.defVar(ncid,'time','double',dimidtime);
        
        %Define the main variable ()
        t_ID = netcdf.defVar(ncid,'tRadEff','double',[dimidlon dimidlat dimidtime]);
        wvlw_ID = netcdf.defVar(ncid,'wvlwRadEff','double',[dimidlon dimidlat dimidtime]);
        wvsw_ID = netcdf.defVar(ncid,'wvswRadEff','double',[dimidlon dimidlat dimidtime]);
        ts_ID = netcdf.defVar(ncid,'tsRadEff','double',[dimidlon dimidlat dimidtime]);
        alb_ID = netcdf.defVar(ncid,'albRadEff','double',[dimidlon dimidlat dimidtime]);
        total_ID = netcdf.defVar(ncid,'totalRadEff','double',[dimidlon dimidlat dimidtime]);
        
        %Define units for the dimension variables
        netcdf.putAtt(ncid,longitude_ID,'standard_name','longitude');
        netcdf.putAtt(ncid,latitude_ID,'standard_name','latitude');
        
        % netcdf.putAtt(ncid,pressure_ID,'standard_name','pressure');
        netcdf.putAtt(ncid,time_ID,'long_name','time');
        netcdf.putAtt(ncid,time_ID,'units','days since 0000-00-00 00:00:0.0');
        netcdf.putAtt(ncid,time_ID,'calendar','gregorian');
        
        netcdf.putAtt(ncid, t_ID,'long_name','Radiation Effect due to Temperature Anomaly (longwave)');
        netcdf.putAtt(ncid, t_ID,'units', 'W m-2');
        
        netcdf.putAtt(ncid,wvlw_ID,'long_name','Radiation Effect due to Specific Humidity Anomaly (longwave)');
        netcdf.putAtt(ncid,wvlw_ID,'units', 'W m-2');
        
        netcdf.putAtt(ncid, wvsw_ID,'long_name','Radiation Effect due to Specific Humidity Anomaly (shortwave');
        netcdf.putAtt(ncid, wvsw_ID,'units', 'W m-2');
        
        netcdf.putAtt(ncid,ts_ID,'long_name','Radiation Effect due to Surface Temperature Anomaly (longwave)');
        netcdf.putAtt(ncid, ts_ID,'units', 'W m-2');
        
        netcdf.putAtt(ncid, alb_ID,'long_name','Radiation Effect due to Albedo Anomaly (shortwave)');
        netcdf.putAtt(ncid, alb_ID,'units', 'W m-2');
        
        netcdf.putAtt(ncid,total_ID,'long_name','Total Radiation Effect');
        netcdf.putAtt(ncid, total_ID,'units', 'W m-2');
        
        %We are done defining the NetCdf
        netcdf.endDef(ncid);
        
        %Then store the dimension variables in
        netcdf.putVar(ncid,longitude_ID,lonf);
        netcdf.putVar(ncid,latitude_ID,latff);
        netcdf.putVar(ncid,time_ID,time);
        
        %Then store my main variable
        
        netcdf.putVar(ncid,t_ID,tEffect);
        netcdf.putVar(ncid,wvlw_ID,wvlwEffect);
        netcdf.putVar(ncid,wvsw_ID,wvswEffect);
        netcdf.putVar(ncid,ts_ID,tsEffect);
        netcdf.putVar(ncid,alb_ID,albEffect);
        netcdf.putVar(ncid,total_ID,totalEffect);
        
        %We're done, close the netcdf
        netcdf.close(ncid)
        
    end
end

