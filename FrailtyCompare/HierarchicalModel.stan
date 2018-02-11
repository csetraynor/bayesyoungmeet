//Example Hierarchical 
data {
  int<lower=0> J; //number of clusters
  real y[J] ; //estimated treatment effects
  real<lower=0> sigma[J]; //s.e. treatment effects
}
parameters {
  real mu;
  real<lower=0> tau;
  real eta[J];
}
transformed parameters{
  real theta[J];
  for (j in 1:J)
    theta[j] = mu + tau * eta[j];
}
model {
  target += normal_lpdf(eta | 0, 1);
  target += normal_lpdf(y | theta, sigma);
}