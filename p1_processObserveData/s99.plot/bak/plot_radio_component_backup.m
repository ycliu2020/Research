% multiply 10a 
for
trend_dR_ts_m = squeeze(mon_trend_dR_ts(:,:,:,1))*365*10;
trend_rhs_m = squeeze(mon_trend_rhs(:,:,:,1))*365*10;
trend_dR_q_m = squeeze(mon_trend_dR_q(:,:,:,1))*365*10;
trend_dR_alb_m = squeeze(mon_trend_dR_alb(:,:,:,1))*365*10;
trend_dR_t_m = squeeze(mon_trend_dR_t(:,:,:,1))*365*10;
trend_dR_cloud_m = squeeze(mon_trend_dR_cloud(:,:,:,1))*365*10;

%season
trend_dR_ts_s = sea_trend_dR_ts(:,:,:)*365*10;
trend_rhs_s = sea_trend_rhs(:,:,:)*365*10;
trend_dR_q_s = sea_trend_dR_q(:,:,:)*365*10;
trend_dR_alb_s = sea_trend_dR_alb(:,:,:)*365*10;
trend_dR_t_s = sea_trend_dR_t(:,:,:)*365*10;
trend_dR_cloud_s = sea_trend_dR_cloud(:,:,:)*365*10;

%year
trend_dR_ts_y = yr_trend_dR_ts(:,:)*365*10;
trend_rhs_y = yr_trend_rhs(:,:)*365*10;
trend_dR_q_y = yr_trend_dR_q(:,:)*365*10;
trend_dR_alb_y = yr_trend_dR_alb(:,:)*365*10;
trend_dR_t_y = yr_trend_dR_t(:,:)*365*10;
trend_dR_cloud_y = yr_trend_dR_cloud(:,:)*365*10;