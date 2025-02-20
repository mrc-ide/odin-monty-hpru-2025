library(monty)

m <- monty_example("banana", sigma = 0.5)

a <- seq(-2, 6, length.out = 1000)
b <- seq(-2.5, 2.5, length.out = 1000)
z <- outer(a, b, function(alpha, beta) {
  exp(monty_model_density(m, rbind(alpha, beta)))
})

point1 <- c(0,0)
point2 <- c(2,1.35)

col1 <- "#33D1FF44"
col2 <- "#F817F144"

image(a, b, z, xlab = "alpha", ylab = "beta", main = "HMC proposal")
points(point1[1], point1[2], col = col1, pch = 19)
points(point2[1], point2[2], col = col2, pch = 19)

n_samples <- 1000

sampler_hmc <- monty_sampler_hmc(epsilon = 0.01, n_integration_steps = 640, debug = T)
sampler_rw <- monty_sampler_random_walk(vcv = diag(2) * 0.3)

image(a, b, z, xlab = "alpha", ylab = "beta", main = "HMC proposal")
res_hmc <- monty_sample(m, sampler_hmc, 1, point1, n_samples)
1-sum(res_hmc$pars["alpha",1,]==point1[1])/n_samples
points(t(res_hmc$pars[,1,]), col = col1, pch = 19)

res_hmc <- monty_sample(m, sampler_hmc, 1, point2, n_samples)
1-sum(res_hmc$pars["alpha",1,]==point2[1])/n_samples
points(t(res_hmc$pars[,1,]), col = col2, pch = 19)

res_sampler_hmc2 <- monty_sampler_hmc(epsilon = 0.2,
                                  n_integration_steps = 5,
                                  debug = T)

res <- monty_sample(m, sampler_hmc2, n_steps = 1, initial = c(0,0))

lines(res$details[[1]]$pars[1,,], res$details[[1]]$pars[2,,], col = "black")



image(a, b, z, xlab = "alpha", ylab = "beta", main = "RW proposal")
res_rw <- monty_sample(m, sampler_rw, 1, point1, n_samples)
1-sum(res_rw$pars["alpha",1,]==point1[1])/n_samples
points(t(res_rw$pars[,1,]), col = col1, pch = 19)

res_rw <- monty_sample(m, sampler_rw, 1, point2, n_samples)
1-sum(res_rw$pars["alpha",1,]==point2[1])/n_samples
points(t(res_rw$pars[,1,]), col = col2, pch = 19)

m <- monty_dsl({
  x ~ Normal(0, 1)
})

sampler_hmc <- monty_sampler_hmc(epsilon = 0.01,
                                 n_integration_steps = 1000,
                                 debug = T)

res <- monty_sample(m, sampler_hmc, n_steps = 1, initial = 0)

plot(res$details[[1]]$pars, type = "l")

library(monty)

set.seed(42)
m <- monty_example("banana", sigma = 0.5)

n_samples <- 1000
point2 <- c(2,1.35)
sampler_hmc <- monty_sampler_hmc(epsilon = 0.1, n_integration_steps = 700, debug = T)
res_hmc <- monty_sample(m, sampler_hmc, 1, point2, n_samples)

#res_hmc <- monty_sample_continue(res_hmc, 2)

plot(NULL, xlim = c(-2, 15), ylim = c(-3, 3))
for(i in 1:n_samples){
  lines(t(res_hmc$details[[i]]$pars[, , 1]))
}

theta <- seq(0, 2 * pi, length.out = 10000)
z95 <- local({
  sigma <- 0.5
  r <- sqrt(qchisq(.95, df = 2))
  x <- r * cos(theta)
  y <- r * sin(theta)
  cbind(x^2 + y * sigma, x)
})

lines(z95[, 1], z95[, 2], col = "red")

beta <- 1.35
alpha <- 2
sigma <- 0.5

-1/2*(beta^2 + ((alpha-beta^2)/sigma)^2) - log(2*pi)
