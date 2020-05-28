% cal real rad
% Sfc_cld/clr net thermal radiation,
% Sfc_cld/clr net short radiation,
% Toa_cld/clr net thermal radiation,
% Toa_cld/clr net short radiation,
for p_1 = 5:5%1 mean amip 2000; 2 mean amip 1980;3 means ssp245, 4 means ssp370; 5 mean amip-hist 2000; 6 mean amip-hist 1980
    nowpath = pwd;
    [readme, Experiment, level, tLin, mPlev, vars] = modelParameters(p_1);

    % inputPath
    inputPath = ['/data1/liuyincheng/cmip6-process/', level.time1{p_1}]; %/data1/liuyincheng/cmip6-process/2000-2014/
    % model loop

    for level1 = 1:length(level.model2)
        % load dvars
        varsPath = [inputPath, level.model2{level1}, '/', level.process3{1}]; %/data1/liuyincheng/cmip6-process/2000-2014/MRI-ESM2-0/rawdata
        % 1.sfc dvars(includ cld and clr sky)
        load([varsPath, 'rlus.mat'])% surface_upwelling_longwave_flux_in_air
        load([varsPath, 'rsus.mat']); load([varsPath, 'rsuscs.mat'])% surface_upwelling_shortwave_flux_in_air
        load([varsPath, 'rsds.mat']); load([varsPath, 'rsdscs.mat'])% surface_downwelling_shortwave_flux_in_air
        load([varsPath, 'rlds.mat']); load([varsPath, 'rldscs.mat'])% surface_downwelling_longwave_flux_in_air
        % 2.toa dvars(includ cld and clr sky)
        load([varsPath, 'rsdt.mat'])% toa_incoming_shortwave_flux
        load([varsPath, 'rsut.mat']); load([varsPath, 'rsutcs.mat'])% toa_outgoing_shortwave_flux
        load([varsPath,  'rlut.mat']); load([varsPath,  'rlutcs.mat'])% toa_outgoing_longwave_flux

        % real rad( 1sfc cld,2sfc clr,3toa cld,4toa clr)
        ll_rad(:, :, :, 1) = rlds - rlus; % Sfc_cld net thermal radiation,
        ll_rad(:, :, :, 2) = rldscs - rlus; % Sfc_clr net thermal radiation,
        ll_rad(:, :, :, 3) = -rlut; % toa_cld net thermal radiation
        ll_rad(:, :, :, 4) = -rlutcs; % toa_clr net thermal radiation

        ss_rad(:, :, :, 1) = rsds - rsus; % Sfc_cld net solar radiation,
        ss_rad(:, :, :, 2) = rsdscs - rsuscs; % Sfc_clr net solar radiation,
        ss_rad(:, :, :, 3) = rsdt - rsut; % toa_cld net solar radiation,
        ss_rad(:, :, :, 4) = rsdt - rsutcs; % toa_clr net solar radiation,
        if all(isnan(ll_rad(:))==0)==0
            disp('warning: ll_rad have NaN data !')
        elseif all(isnan(ss_rad(:))==0)==0
            disp('warning: ss_rad have NaN data !')
        end
        

    end


end