set.seed(342)

df_draws = tbl_df(as.data.frame(simulated_fit))


beta = df_draws %>% select(starts_with("beta") )
lambda = df_draws %>% select(starts_with("lambda"))

pp_newdata <- mapply(function(x, y)  { pem_sim_data( beta = x,
                                  lambda = y,
                                  n = test_n,
                                  tau = seq(0, 160, by = 20), 
                                  X = matrix(c(rnorm(100),
                                               sample(c(0,1), 100, replace = T))
                                             , ncol=2)
)}, x = as.list(beta), y = as.list(lambda)
  )
