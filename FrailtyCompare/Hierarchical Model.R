#create data J = cluster , y = response, sigma = known standard deviation
cluster_linear_data <- list(J = 8, 
                            y = c(28, 8, -3, 7, -1, 1, 18, 12),
                            sigma = c(15, 10, 16, 11, 9, 11, 10, 18))

fit <- stan(file = 'HierarchicalModel.stan', data = cluster_linear_data,
            iter = 1000, chains = 4)
print(fit)
plot(fit)
pairs(fit, pars = c('mu', 'tau', 'lp__'))

la <- extract(fit, permuted = TRUE)
mu <- la$mu

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit, permuted = FALSE) 

### use S3 functions as.array (or as.matrix) on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)

print(fit, digits = 1)


##Rats example
y <- as.matrix(read.table('https://raw.github.com/wiki/stan-dev/rstan/rats.txt', header = TRUE))
x <- c(8, 15, 22, 29, 36)
xbar <- mean(x)
N <- nrow(y)
T <- ncol(y)
rats_fit <- stan(file = 'https://raw.githubusercontent.com/stan-dev/example-models/master/bugs_examples/vol1/rats/rats.stan')
