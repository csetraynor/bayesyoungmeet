##Parameterization
require(ggplot2)
require(scales)

#Consider the resulting distribution of alpha over a range of values for alpha_raw sampled from our normal(0, 1) prior:
alpha_raw <- rnorm(1000, 0, 1)
tau_al <- 10
log_alpha <- alpha_raw * tau_al
alpha <- exp(log_alpha)
ggplot(data.frame(alpha = alpha, alpha_raw = alpha_raw), 
       aes(x = alpha)) + 
  geom_density() + 
  scale_x_log10(labels = scientific)

#This distribution is centered at 0 and has more consistent behavior throughout its range of values.
ggplot(data.frame(alpha = alpha, alpha_raw = alpha_raw), 
       aes(x = alpha, y = alpha_raw)) + 
  geom_density2d() + 
  scale_x_log10(labels = scientific)
