  matrix[T, N] y_rep;
  matrix[T, N] effect;
  for (i in 1:N) {
    y_rep[, i] = to_vector(normal_rng(y_mean[, i], sigma[i]));
    effect[, i] = y[i] - y_rep[, i];
  }
  vector[N] R2;
  vector[N] RMSE_pre;
  vector[N] RMSE_post;
  vector[N] RMSE_ratio;
  vector[N] effective_donors;
  real avg_effect_pre = 0;
  real avg_effect_post = 0;
  for (i in 1:N) {
    avg_effect_pre += sum(effect[1:T_pre[i], i]);
    avg_effect_post += sum(effect[(T_pre[i] + 1):T, i]);
    real var_res = square(sigma[i]);
    real var_fit = variance(y_mean[1:T_pre[i], i]);
    R2[i] = var_fit / (var_fit + var_res);
    RMSE_pre[i] = sqrt(mean(effect[1:T_pre[i], i]^2));
    RMSE_post[i] = sqrt(mean(effect[(T_pre[i] + 1):T, i]^2));
    RMSE_ratio[i] = RMSE_pre[i] / RMSE_post[i];
    effective_donors[i] = inv(sum(square(omega[i])));
  }
  avg_effect_pre /= sum(T_pre);
  avg_effect_post /= (N * T - sum(T_pre));
  real avg_RMSE_pre = mean(RMSE_pre);
  real avg_RMSE_post = mean(RMSE_post);
  real avg_RMSE_ratio = mean(RMSE_ratio);
  real avg_effective_donors = mean(effective_donors);
