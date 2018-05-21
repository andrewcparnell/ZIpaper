# Example 1 from ZI paper 2018. Fish data set. 

# Boiler plate code. Clear workspace and load in packages
rm(list = ls())

# Install auf if not already installed
# devtools::install_github('andrewcparnell/auf)
library(auf)
packages('tidyverse', 'boot', 'rstan', 'gridExtra',
         'ggpubr', 'bayesplot', 'mgcv', 'loo')
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Load in the data
zinb = read.csv("fish.csv")[-c(89, 160),]

# Explanatory plots -------------------------------------------------------

# Create an explanatory plot
p1 = ggplot(zinb, aes(x = count)) + 
  theme_bw() +
  geom_bar(aes(fill = log(..count..))) + 
  #geom_histogram(aes(fill = log(..density..))) + 
  theme(legend.position = 'none') + 
  labs(x = 'Count of fish caught',
       y = 'Number of observations')

p2 = ggplot(zinb, aes(x = as.factor(persons), y = count, fill = (persons))) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_boxplot() + 
  labs(y = 'Count of fish caught',
       x = 'Number of people on trip')

p3 = ggplot(zinb, aes(x = as.factor(camper), 
                      y = count, fill = camper)) +
  theme_bw() +
  theme(legend.position = 'none') + 
  geom_boxplot() + 
  labs(y = 'Count of fish caught',
       x = 'Camper van present')

ggarrange(p1, p2, p3, ncol=3, nrow=1)
ggsave(file = 'fish_fig_1.pdf', width = 8, height = 4)

# First stan run ----------------------------------------------------------

# Run the first stan model - compare the three versions of hurdle, additive mixture
zi1 = stanc(file = 'model1.stan') # Check Stan file
zi1_run = stan_model(stanc_ret = zi1) # Compile Stan code
# Prep data
y = zinb$count
x = with(zinb, cbind(1, 
                     persons - mean(persons),
                     camper - mean(camper)))
K = ncol(x)
x_cov = uniquecombs(x)
N_cov = nrow(x_cov)
fit_mult = sampling(zi1_run,
                    data = list(y = y,
                                tau1 = 0,
                                tau2 = 0,
                                x = x,
                                K = K,
                                N = nrow(x),
                                x_cov = x_cov,
                                N_cov = N_cov))
fit_add = sampling(zi1_run,
                   data = list(y = y,
                               tau1 = 0,
                               tau2 = -1,
                               x = x,
                               K = K,
                               N = nrow(x),
                               x_cov = x_cov,
                               N_cov = N_cov))
fit_hurd = sampling(zi1_run,
                    data = list(y = y,
                                tau1 = 1,# Hurdle is (1,-1) not (-1,-1)
                                tau2 = -1,
                                x = x,
                                K = K,
                                N = nrow(x),
                                x_cov = x_cov,
                                N_cov = N_cov))

# Diagnostics for model runs
important_pars = c('alpha', 'beta')
post_mult = as.array(fit_mult)
post_add = as.array(fit_add)
post_hurd = as.array(fit_hurd)
mcmc_pairs(post_mult, regex_pars = important_pars,
           off_diag_args = list(size = 1, alpha = 0.5))
mcmc_pairs(post_add, regex_pars = important_pars,
           off_diag_args = list(size = 1, alpha = 0.5))
mcmc_pairs(post_hurd, regex_pars = important_pars,
           off_diag_args = list(size = 1, alpha = 0.5))

# Now use loo package to get WAIC values
log_lik1 = extract_log_lik(fit_mult)
log_lik2 = extract_log_lik(fit_add)
log_lik3 = extract_log_lik(fit_hurd)
waic_mult = waic(log_lik1)
waic_add = waic(log_lik2)
waic_hurd = waic(log_lik3)
print(compare(waic_mult, waic_add, waic_hurd), digits = 2)

# Comparison plot of pi0 vs pit0
lpi0 = rstan::extract(fit_hurd, pars = 'lpi0')$lpi0
lpit0 = rstan::extract(fit_hurd, pars = 'lpit0')$lpit0
lpi0_quantiles = apply(lpi0, 2, 'quantile', probs = c(25,50,75)/100)
lpit0_quantiles = apply(lpit0, 2, 'quantile', probs = c(25,50,75)/100)

df = data.frame(pi0_50 = exp(lpi0_quantiles[2,]),
                pi0_25 = exp(lpi0_quantiles[1,]), 
                pi0_75 = exp(lpi0_quantiles[3,]), 
                pit0_50 = exp(lpit0_quantiles[2,]),
                pit0_25 = exp(lpit0_quantiles[1,]), 
                pit0_75 = exp(lpit0_quantiles[3,]), 
                persons = as.factor(zinb$persons),
                camper = as.factor(zinb$camper))
df2 = data.frame(logit_pi0_50 = logit(exp(lpi0_quantiles[2,])),
                 logit_pi0_25 = logit(exp(lpi0_quantiles[1,])), 
                 logit_pi0_75 = logit(exp(lpi0_quantiles[3,])), 
                 logit_pit0_50 = logit(exp(lpit0_quantiles[2,])),
                 logit_pit0_25 = logit(exp(lpit0_quantiles[1,])), 
                 logit_pit0_75 = logit(exp(lpit0_quantiles[3,])), 
                 persons = as.factor(zinb$persons),
                 camper = as.factor(zinb$camper))
p1 = ggplot(df, aes(pi0_50, pit0_50, 
                    colour = persons, shape = camper)) + 
  ylim(0, 1) + 
  geom_point() +
  theme_bw() +
  geom_errorbar(aes(ymin = pit0_25, ymax = pit0_75)) + 
  geom_errorbarh(aes(xmin = pi0_25, xmax = pi0_75)) + 
  labs(x = expression(pi["0"]),
       y = expression(tilde(pi)["0"]))

p2 = ggplot(df2, aes(x = logit_pi0_50, logit_pit0_50,
                     colour = persons, shape = camper)) + 
  geom_point() +
  theme_bw() +
  geom_errorbar(aes(ymin = logit_pit0_25, ymax = logit_pit0_75)) + 
  geom_errorbarh(aes(xmin = logit_pi0_25, xmax = logit_pi0_75)) + 
  labs(x = expression(paste(logit,"(" , pi["0"],")")),
       y = expression(paste(logit,"(" ,tilde(pi)["0"],")")))

ggarrange(p1, p2, ncol=2, nrow=1, 
          common.legend = TRUE, legend="bottom")
ggsave(file = 'fish_fig_2.pdf', width = 8, height = 4)

# Compare with empirical --------------------------------------------------

# Calculate empirical probabilities of zero
# Need the proportions of zeroes for each combination of persons and camper
ans = aggregate(zinb$count, by = list(zinb$persons, zinb$camper), function(x) length(x[x==0])/length(x))
colnames(ans) = c('persons', 'camper', 'p0')
ans$persons_std = ans$persons - mean(zinb$persons)
ans$camper_std = ans$camper - mean(zinb$camper)

ans2 = as.data.frame(x_cov)
ans2$pi0_est = exp(lpi0_quantiles[2,])
ans2$pit0_est = exp(lpit0_quantiles[2,])
colnames(ans2) = c('const', 'persons_std', 'camper_std', 'pi0', 'pit0')

ans3 = merge(ans, ans2)

p3 = ggplot(ans3, aes(x = pi0, y = p0)) + geom_point() + 
  xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1)
p4 = ggplot(ans3, aes(x = pit0, y = p0)) + geom_point() + 
  xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1)
ggarrange(p3, p4, ncol=2, nrow=1)

lpi0 = rstan::extract(fit_hurd, pars = 'lpi0')$lpi0
lpit0 = rstan::extract(fit_hurd, pars = 'lpit0')$lpit0
lpi0_hurd_quantiles = apply(lpi0, 2, 'quantile', probs = c(25,50,75)/100)
lpit0_hurd_quantiles = apply(lpit0, 2, 'quantile', probs = c(25,50,75)/100)
ans4 = as.data.frame(x_cov)
ans4$pi0_est = exp(lpi0_hurd_quantiles[2,])
ans4$pit0_est = exp(lpit0_hurd_quantiles[2,])
colnames(ans4) = c('const', 'persons_std', 'camper_std', 'pi0', 'pit0')

ans5 = merge(ans, ans4)


p5 = ggplot(ans3, aes(x = pi0, y = p0)) + geom_point() + 
  geom_point(data = ans5, aes(x = pi0, y = p0, colour = 'red')) +
  xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1)
p6 = ggplot(ans3, aes(x = pit0, y = p0)) + geom_point() + 
  geom_point(data = ans5, aes(x = pit0, y = p0, colour = 'red')) +
  xlim(0, 1) + ylim(0, 1) + geom_abline(intercept = 0, slope = 1)
ggarrange(p5, p6, ncol=2, nrow=1)



# NOTE: CODE FROM HEREON EXPERIMENTAL -------------------------------------
# 
# # Second example where tau types are learned per covariate ----------------
# 
# zi4 = stanc(file = 'model4.stan') # Check Stan file
# zi4_run = stan_model(stanc_ret = zi4) # Compile Stan code
# 
# # Sort out data
# y = zinb$count
# x = with(zinb, cbind(1, 
#                      persons - mean(persons),
#                      camper - mean(camper)))
# K = ncol(x)
# dat = list(y = y,
#            x = x,
#            K = K,
#            N = nrow(x))
# my_init = function() {
#   list(tau1p1 = 0.5, tau2p1 = 0.5, alpha = 0, beta = c(0, 0, 0))
# } 
# fit = sampling(zi4_run,
#                data = dat,
#                init = my_init,
#                chain = 1)
# 
# # Second model with multiple types ----------------------------------------
# 
# # Set up stan
# zi2 = stanc(file = 'model2.stan') # Check Stan file
# zi2_run = stan_model(stanc_ret = zi2) # Compile Stan code
# 
# # Get all models remember it's (0,0), (0, -1) and (-1, -1) 
# # for mult, add, and hurdle respectively
# # so if there are four levels of variable persons the combinations are
# # 3^4 = 81 different models!
# all_models = expand.grid(1:3,1:3,1:3,1:3)
# waic_store = matrix(NA, ncol = 2, nrow = nrow(all_models))
# # Let 1 = mult, 2 = add, 3 = hurd
# y = zinb$count
# x = with(zinb, cbind(1, 
#                      persons - mean(persons),
#                      camper - mean(camper)))
# K = ncol(x)
# dat = list(y = y,
#            x = x,
#            K = K,
#            N = nrow(x))
# tau1_all = c(0, 0, -1)
# tau2_all = c(0, -1, -1)
# for(i in 1:nrow(all_models)) {
#   print(i)
#   # Prep data
#   curr_dat = dat
#   curr_dat$tau1 = rep(NA, length(dat$y))
#   curr_dat$tau2 = rep(NA, length(dat$y))
#   for(j in 1:length(dat$y)) {
#     curr_dat$tau1[j] = tau1_all[all_models[i,zinb$persons[j]]]
#     curr_dat$tau2[j] = tau2_all[all_models[i,zinb$persons[j]]]
#   }
#   
#   curr_fit = sampling(zi2_run,
#                       data = curr_dat)
#   curr_log_lik = extract_log_lik(curr_fit)
#   curr_waic = waic(curr_log_lik)
#   waic_store[i,] = c(curr_waic$waic, curr_waic$se_waic)
# }
# 
# write.csv(waic_store, file = 'waic_store.csv', 
#           quote = FALSE, row.names = FALSE)