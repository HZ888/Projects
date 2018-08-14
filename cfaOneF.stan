data
{
  int ny; //# of variables
  int nsub; //# of subjects
  matrix[nsub,ny] y; // matrix of outcome
  int nfac;
  int<lower=1,upper=nfac> y_vs_fac[ny]; // y_vs_f: list of factors associated with outcome
}
transformed data
{
  matrix [nsub,ny] stand_y; //standardized outcome
  for (i in 1:ny)
  {
    stand_y[,i] = (y[,i]-mean(y[,i])) / sd(y[,i]);
  }
}
parameters
{
  matrix[nfac,nsub] normal01; //
  cholesky_factor_corr[nfac] fac_cor_helper; //correlations among factors
  vector[ny] stand_y_mean;
  vector<lower=0>[ny] stand_y_error;
  vector<lower=0, upper=1>[ny] betas;
}
transformed parameters
{
  matrix[nsub,nfac] sub_fac;
  sub_fac=transpose(fac_cor_helper*normal01);
}
model
{
  to_vector(normal01) ~ normal(0,1);
  fac_cor_helper ~ lkj_corr_cholesky(1); //flat prior on correlations among factors
  stand_y_mean ~ normal(0,1);
  stand_y_error ~ weibull(2,1);
  for (i in 1:ny)
  {
    stand_y[,i] ~ normal(stand_y_mean[i] + sub_fac[,y_vs_fac[i]] * betas[i], stand_y_error[i]);
    betas[i] ~ beta(700,200);
    //betas[i] ~ normal(0.2,0.01);
  }
}
generated quantities
{
  corr_matrix[nfac] fac_cor;
  vector[ny] y_means;
  vector[ny] y_error;
  fac_cor = multiply_lower_tri_self_transpose(fac_cor_helper);
  for (i in 1:ny)
  {
    y_means[i] = stand_y_mean[i] * sd(y[,i])+mean(y[,i]);
    y_error[i] = stand_y_error[i] * sd(y[,i]);
  }
}
