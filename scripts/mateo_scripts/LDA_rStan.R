library(tm)
library(rstan)
library(slam)
library(shinystan)
library(topicmodels)

pkgbuild::has_build_tools(debug = TRUE)

reordered_ascat_for_modeling <- readRDS("STAN/reordered_ascat_for_modeling.rds")

# subseting the data to work with smaller samples sizes
# ids <- unique(tmp_data$ID)[1:20] # take only first 10 sample ids # nolint
# tmp_data <- tmp_data[which(tmp_data$ID %in% ids), ] # nolint

unique(reordered_ascat_for_modeling$distToNearestCNVQ)
tmp_data <- reordered_ascat_for_modeling

num_sigs <- 5
chains <- 4
niter <- 2000
# read and run stan code
fit <- stan(file = "STAN/LDA.stan", data = list(
  K = num_sigs,
  V_f1 = 7,
  V_f2 = 2,
  V_f3 = 6,
  V_f4 = 7,
  V_f5 = 7,
  V_f6 = 7,
  V_f7 = 3,
  V_f8 = 5,
  M = length(unique(tmp_data$ID)),
  N = nrow(tmp_data),
  doc = as.integer(as.factor(tmp_data$ID)),
  feature1 = as.integer(tmp_data$segmentSizeQ),
  feature2 = as.integer(tmp_data$isLOH),
  feature3 = as.integer(tmp_data$cnQ),
  feature4 = as.integer(tmp_data$distToCentromereQ),
  feature5 = as.integer(tmp_data$distToTelomereQ),
  feature6 = as.integer(tmp_data$changepointCNQ),
  feature7 = as.integer(tmp_data$sizeDiploidSegQ),
  feature8 = as.integer(tmp_data$distToNearestCNVQ),
  alpha = rep(1, num_sigs),
  beta_f1 = rep(1, 7),
  beta_f2 = rep(1, 2),
  beta_f3 = rep(1, 6),
  beta_f4 = rep(1, 7),
  beta_f5 = rep(1, 7),
  beta_f6 = rep(1, 7),
  beta_f7 = rep(1, 3),
  beta_f8 = rep(1, 5)
), chains = chains, iter = niter, verbose = T, cores = 2000)

saveRDS(fit, "STAN/stan_fit_6sigs_2kIter_4chains.rds")
