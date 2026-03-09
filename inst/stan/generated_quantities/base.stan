  vector[T] effect = y - synthetic_y;
  real avg_effect_pre = mean(head(effect, T_pre));
  real avg_effect_post = mean(tail(effect, T - T_pre));
  vector[T - T_pre] avg_effect_post_cumulative = 
  cumulative_sum(tail(effect, T - T_pre)) ./ tt;
  vector[T - T_pre] relative_change = 
  tail(synthetic_y, T - T_pre) ./ synthetic_y[T_pre];
  
  real RMSE_pre = sqrt(mean(head(effect, T_pre)^2));
  real RMSE_post = sqrt(mean(tail(effect, T - T_pre)^2));
  real RMSE_ratio = RMSE_post / RMSE_pre;
  real effective_donors = inv(sum(square(omega)));
