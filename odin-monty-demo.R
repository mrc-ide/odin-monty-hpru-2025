## Read in incidence data; see data/incidence.csv for original
d <- data.frame(
  time = 1:20,
  cases = c(12, 23, 25, 36, 30, 57, 59, 62, 47, 52, 56, 33, 34, 19, 27,
            25, 15, 20, 11, 7))

## system of equations written in odin, compile into dust
## Note that dt (the step size) is a special parameter
sir <- odin2::odin({
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(incidence) <- incidence + n_SI
  
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  
  initial(S) <- N - I0
  initial(I) <- I0
  initial(R) <- 0
  initial(incidence, zero_every = 1) <- 0
  
  N <- parameter(1000)
  I0 <- parameter(10)
  beta <- parameter(0.2)
  gamma <- parameter(0.1)
  
  cases <- data()
  cases ~ Poisson(incidence)
})

## Run a single simulation
sys <- dust2::dust_system_create(sir, pars = list(), dt = 0.25)
dust2::dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust2::dust_system_simulate(sys, t)
y <- dust2::dust_unpack_state(sys, y)

plot(t, y$incidence, type = "l", xlab = "Time", ylab = "Infection incidence")

## Run multiple simulations
sys <- dust2::dust_system_create(sir, pars = list(), n_particles = 50, 
                                 dt = 0.25)
dust2::dust_system_set_state_initial(sys)
t <- seq(0, 100)
y <- dust2::dust_system_simulate(sys, t)
y <- dust2::dust_unpack_state(sys, y)
matplot(t, t(y$incidence), type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infected population")



## Create a particle filter to fit to the data
filter <- dust2::dust_filter_create(sir, time_start = 0,
                                    data = d, dt = 0.25,
                                    n_particles = 200)
pars <- list(beta = 1, gamma = 0.6)
## Each time the filter is run it produces a different estimate of
## the marginal log-likelihood
replicate(10, dust2::dust_likelihood_run(filter, pars))


## Create a packer to "unpack" from a vector of fitted parameters
## to a list of parameters for input in our dust system
packer <- monty::monty_packer(scalar = c("beta", "gamma"),
                              fixed = list(I0 = 10, N = 1000))
## Create a monty model using our filter
likelihood <- dust2::dust_likelihood_monty(filter, packer,
                                           save_trajectories = TRUE)

## Create a monty model for prior distributions via the monty DSL
prior <- monty::monty_dsl({
  beta ~ Exponential(mean = 1)
  gamma ~ Exponential(mean = 0.5)
})

## Combine the two models (on log-scale so we add them) into one monty model
posterior <- likelihood + prior

## Create a Metropolis-Hastings random walk sampler
vcv <- diag(2) * 0.01
sampler <- monty::monty_sampler_random_walk(vcv)

## Run PMCMC
samples <- monty::monty_sample(posterior, sampler, n_steps = 1000,
                               initial = c(0.3, 0.1),
                               n_chains = 4)

## Can get diagnostics via the posterior package
samples_df <- posterior::as_draws_df(samples)
posterior::summarise_draws(samples_df)

## Thin samples, and then plot fitted trajectories against the data
samples <- monty::monty_samples_thin(samples, thinning_factor = 4, burnin = 500)
y <- dust2::dust_unpack_state(filter,
                              samples$observations$trajectories)
incidence <- array(y$incidence, c(20, 500))
matplot(d$time, incidence, type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infected population", ylim = c(0, 75))
points(cases ~ time, d, pch = 19, col = "red")





#### Fitting a deterministic "expectation" version of the model
#### (All random draws in odin code replaced by mean of the distribution)

## Use "unfilter" - runs only a single "particle"
filter <- dust2::dust_unfilter_create(sir, time_start = 0,
                                      data = d, dt = 0.25)
pars <- list(beta = 1, gamma = 0.6)
## Each run produces the same marginal log-likelihood
replicate(10, dust2::dust_likelihood_run(filter, pars))


## Remaining setup the same as before
packer <- monty::monty_packer(scalar = c("beta", "gamma"),
                              fixed = list(I0 = 10, N = 1000))
likelihood <- dust2::dust_likelihood_monty(filter, packer,
                                           save_trajectories = TRUE)

prior <- monty::monty_dsl({
  beta ~ Exponential(mean = 1)
  gamma ~ Exponential(mean = 0.5)
})

posterior <- likelihood + prior
## We will use smaller variance here
vcv <- diag(2) * 0.001

sampler <- monty::monty_sampler_random_walk(vcv)
samples_det <- monty::monty_sample(posterior, sampler, n_steps = 1000,
                                   initial = c(0.3, 0.1),
                                    n_chains = 4)
samples_df <- posterior::as_draws_df(samples_det)
posterior::summarise_draws(samples_df)

samples_det<- monty::monty_samples_thin(samples_det, thinning_factor = 4,
                                        burnin = 500)
y <- dust2::dust_unpack_state(filter,
                              samples_det$observations$trajectories)
incidence <- array(y$incidence, c(20, 500))
matplot(d$time, incidence, type = "l", lty = 1, col = "#00000044",
        xlab = "Time", ylab = "Infected population", ylim = c(0, 75))
points(cases ~ time, d, pch = 19, col = "red")


pars_stochastic <- array(samples$pars, c(2, 500))
pars_deterministic <- array(samples_det$pars, c(2, 500))
plot(pars_stochastic[1, ], pars_stochastic[2, ], ylab = "gamma", xlab = "beta",
     pch = 19, col = "blue")
points(pars_deterministic[1, ], pars_deterministic[2, ], pch = 19, col = "red")
legend("bottomright", c("stochastic fit", "deterministic fit"), pch = c(19, 19), 
       col = c("blue", "red"))
