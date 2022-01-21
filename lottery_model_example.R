# load packages -----------------------------------------------------------
library(tidyverse)
library(shades)
library(latex2exp)
library(extrafont)

# simulation step function ------------------------------------------------

step <- function(n_mat, time) {
  # performs one time step of the spatiotemporal lottery model
  #
  # returns: a list with species densities, environmental parameters, a shared competition parameters.
  
  # get environmental and competition params for each species
  num_spp <- NCOL(n_mat) # third dimension indexes species
  #E_mat <- E_arr[i,,]
  #E_mat <- sapply(seq_len(num_spp), function(i) get_E())
  E_mat <-t(E_arr[,,time])
  C_vec <- log( rowSums(exp(E_mat)*n_mat) / rowSums(d_mat*n_mat) )
  C_mat <- matrix(rep(C_vec, num_spp), nrow = K)
  
  # time step 
  lam <- (1-d_mat) + exp(E_mat - C_mat) # finite rates of increase
  n_prime <- n_mat*lam # pop density after local growth 
  n_dispersal <- n_prime * dispersal_mat # dispersers (in density units) for each each location
  n_tp1 <- n_prime - n_dispersal + matrix(rep(colSums(n_dispersal)/K, K),nrow = K, byrow = TRUE) # pop density at time t + 1 is equal to population density after local growth, minus dispersal of local population, plus the contribution from global dispersal. 
  
  
  inv_pop_size <- sum(n_tp1[,inv_index])
  if(inv_pop_size >= inv_ceiling)
  {
    n_tp1[,inv_index] <- n_tp1[,inv_index] * (inv_ceiling / inv_pop_size)
  }
  
  return(list(n_mat = n_tp1, E_mat = E_mat, C_mat = C_mat, C_vec = C_vec))
}


# simulate and record E & C -----------------------------------------------

# simulation params

get_data <- function(n_mat_init, burn_in, record_time)
{
  t_counter <- 1
  # perform burn_in
  n_mat <- n_mat_init
  for(i in seq_len(burn_in)) {
    dat <- step(n_mat, t_counter)
    t_counter <- t_counter + 1
    n_mat <- dat[["n_mat"]]
  }
  
  
  # preallocate for each species
  spp_dat_list <- lapply(seq_len(num_spp), function(i) list(
    n = matrix(NA, nrow = K, ncol = record_time),
    E = matrix(NA, nrow = K, ncol = record_time),
    C = matrix(NA, nrow = K, ncol = record_time)))
  
  for(i in seq_len(record_time)) {
    dat <- step(n_mat, t_counter)
    t_counter <- t_counter + 1
    n_mat <- dat[["n_mat"]]
    for(j in seq_len(num_spp))
    {
      spp_dat_list[[j]]$n[,i] <- dat$n_mat[,j]
      spp_dat_list[[j]]$E[,i] <- dat$E_mat[,j]
      spp_dat_list[[j]]$C[,i] <- dat$C_vec
    }
  }
  
  return(spp_dat_list)
}

# get the equilibrium values of E and C -----------------------------------
# mean and sd lognormal parameters from which per-capita fecundities are drawn. For simplicity, we give these parameters to both species, but define each species' fecundities to be independent of the other

# 
# C_inv1_mean <- mean( get_data(n1_init = 0, n2_init = 1, burn_in, record_time)[["C_dat"]])  
# C_inv2_mean <- mean( get_data(n1_init = 1, n2_init = 0, burn_in, record_time)[["C_dat"]]) 
# 
# # define the equilbrium competition as the average competition experienced by the invader, which we then average across all species in the invader state
# C_star <- mean(c(C_inv1_mean, C_inv2_mean))
# 
# # use equilibrium C to solve for equilbrium E
# E1_star <- log(d1) + C_star 
# E2_star <- log(d2) + C_star 


# functions for calculating coexistence mechanisms ------------------------

lambda <- function(d, E, C){
  1 - d + exp(E-C)
}

lambdaTwiddle <- function(n, d, E, C){
  
  if(is.vector(n) && is.vector(E) && is.vector(C))
  {
    nu <- n/mean(n)
  } else
  {
    nu <- sapply(seq_len(ncol(n)), function(i) n[,i]/mean(n[,i]))
  }
  nu[which(!is.finite(nu), arr.ind = TRUE)] <- 1 # if spatial mean(n) = 0  at any time step, we will get Inf or NaNs. replace with 1.
  
  return(mean(nu * lambda(d, E, C)))
}

Exp <- function(dims = c("NA", "NA"), X){
  if(!("x" %in% dims) && !("t" %in% dims))
  {
    stop("must select an dimension to average over")
  } else if(is.vector(X))
  {
    return(mean(X))
  } else 
  {
    if("x" %in% dims && "t" %in% dims)
    {
      return(mean(as.vector(X)))
    } else if("x" %in% dims)
    {
      return(colMeans(X))
    } else if("t" %in% dims)
    {
      return(rowMeans(X))
    }else
    {
      stop("somthing is wrong; nothing happened")
    }
  }
}

Var <- function(dims = c("NA", "NA"), X){
  if(!("x" %in% dims) && !("t" %in% dims))
  {
    stop("must select an dimension to average over")
  } else if(is.vector(X))
  {
    return(var(X))
  } else 
  {
    if("x" %in% dims && "t" %in% dims)
    {
      return(var(as.vector(X)))
    } else if("x" %in% dims)
    {
      return(apply(X, MARGIN = 2, FUN = function(x) var(x)))
    } else if("t" %in% dims)
    {
      return(apply(X, MARGIN = 1, FUN = function(x) var(x)))
    }else
    {
      stop("somthing is wrong; nothing happened")
    }
  }
}

Cov <- function(dims = c("NA", "NA"), X, Y){
  if(!("x" %in% dims) && !("t" %in% dims))
  {
    stop("must select an dimension to average over")
  } else if(is.vector(X) && is.vector(Y))
  {
    return(cov(X,Y))
  } else 
  {
    if("x" %in% dims && "t" %in% dims)
    {
      return(cov(as.vector(X), as.vector(Y)))
    } else if("x" %in% dims)
    {
      return(sapply(seq_len(NCOL(X)), function(i) cov(X[,i], Y[,i])))
    } else if("t" %in% dims)
    {
      return(sapply(seq_len(NROW(X)), function(i) cov(X[i,], Y[i,])))
    } else
    {
      stop("somthing is wrong; nothing happened")
    }
  }
}

shuffle <- function(X)
{
  if(is.vector(X))
  {
    return(sample(X))
  } else
  {
    mat_dims <- dim(X)
    return(matrix(sample(X), nrow = mat_dims[1], ncol = mat_dims[2]))
  }
}

# functions for small-noise coexistence mechanisms and intermediaries -----------------

rPrimePre <- function(s) s$alphaSup1 * (Exp(c("x","t"), s$E) - s$E_star) + (1/2)*s$alphaSup2*Var(c("x","t"), s$E) - (1/2)*(s$alphaSup1^2)*Var("t", Exp("x", s$E))
rPrimePre_A <- function(s) s$alphaSup1 * (Exp(c("x","t"), s$E) - s$E_star)
rPrimePre_S <- function(s) (1/2)*(s$alphaSup2 * Var("x", Exp("t", s$E)))
rPrimePre_T <- function(s) (1/2)*(s$alphaSup2 - s$alphaSup1^2) * Var("t", Exp("x", s$E))
rPrimePre_R <- function(s){
  rPrimePre(s) - (rPrimePre_A(s) + rPrimePre_S(s) + rPrimePre_T(s))
}

deltaRhoPre <- function(s) s$betaSup1 * (Exp(c("x","t"), s$C) - s$C_star)

deltaNPre <- function(s) (1/2)*s$betaSup2*Var(c("x","t"), s$C) - (1/2)*(s$betaSup1^2)*Var("t", Exp("x", s$C))
deltaNPre_S <- function(s) (1/2)*(s$betaSup2 * Var("x", Exp("t", s$C)))
deltaNPre_T <- function(s) (1/2)*(s$betaSup2 - s$betaSup1^2) * Var("t", Exp("x", s$C))
deltaNPre_R <- function(s){
  deltaNPre(s) - (deltaNPre_S(s) + deltaNPre_T(s))
} 

s <- spp_dat_list[[1]]
s <- spp_dat_list[[2]]
s <- spp_dat_list[[3]]

deltaIPre <- function(s) s$zeta*Cov(c("x","t"), s$E, s$C) - s$alphaSup1*s$betaSup1*Cov("t", Exp("x", s$E), Exp("x", s$C))
deltaIPre_S <- function(s) s$zeta*Cov("x", Exp("t", s$E), Exp("t", s$C))
deltaIPre_T <- function(s) (s$zeta - s$alphaSup1*s$betaSup1)*Cov("t", Exp("x", s$E), Exp("x", s$C))
deltaIPre_R <- function(s) {
  deltaIPre(s) - (deltaIPre_S(s) + deltaIPre_T(s))
} 

deltaKappaPre <- function(s){
  nu <- sapply(seq_len(ncol(s$n)), function(i) s$n[,i]/mean(s$n[,i]))
  nu[which(!is.finite(nu), arr.ind = TRUE)] <- 1 # if spatial mean(n) = 0  at any time step, we will get Inf or NaNs. replace with 1.
  lambda_bit <- s$alphaSup1*(s$E - s$E_star) + s$betaSup1*(s$C-s$C_star)
  return(Exp("t", Cov("x", nu, lambda_bit)))
}

partition <- function(quantity_name, spp)
{
  val <-  quantity_name(spp[[inv_index]])
  
  res_list <- spp[-inv_index]
  num_residents <- length(res_list)
  for(j in seq_along(res_list))
  {
    val <- val - (1/num_residents) * res_list[[j]]$qij * quantity_name(res_list[[j]])
  }
  return(val)
}


# functions for exact coexistence mechanisms and intermediaries -----------------

curlyEBarPrime <- function(s)  Exp("t", log(Exp("x", lambda(s$d, s$E, s$C_star))))
curlyEBarPrime_A <- function(s) log(lambda(s$d, Exp(c("x","t"), s$E), s$C_star)) 
curlyEBarPrime_S <- function(s) log(Exp("x", lambda(s$d, Exp("t", s$E), s$C_star))) - curlyEBarPrime_A(s)
curlyEBarPrime_T <- function(s) Exp("t", log(lambda(s$d, Exp("x", s$E), s$C_star))) - curlyEBarPrime_A(s)
curlyEBarPrime_R <- function(s){
  curlyEBarPrime(s) - (curlyEBarPrime_A(s) + curlyEBarPrime_S(s) + curlyEBarPrime_T(s))
}

curlyCBarPrime <- function(s)  Exp("t", log(Exp("x", lambda(s$d, s$E_star, s$C))))
curlyCBarPrime_A <- function(s) log(lambda(s$d, s$E_star, Exp(c("x","t"), s$C))) 
curlyCBarPrime_S <- function(s) log(Exp("x",lambda(s$d, s$E_star, Exp("t", s$C)))) - curlyCBarPrime_A(s)
curlyCBarPrime_T <- function(s) Exp("t", log( lambda(s$d, s$E_star, Exp("x", s$C)))) - curlyCBarPrime_A(s)
curlyCBarPrime_R <- function(s){
  curlyCBarPrime(s) - (curlyCBarPrime_A(s) + curlyCBarPrime_S(s) + curlyCBarPrime_T(s))
}

curlyIBarPrime <- function(s) Exp("t", log(Exp("x", lambda(s$d, s$E, s$C)))) - Exp("t", log(Exp("x", lambda(s$d, shuffle(s$E), s$C))))
curlyIBarPrime_A <- function(s) log( lambda(s$d,  Exp(c("x","t"), s$E), Exp(c("x","t"),s$C))) - (curlyEBarPrime_A(s) + curlyCBarPrime_A(s))

curlyIBarPrime_S <- function(s) log(Exp("x", lambda(s$d, Exp("t", s$E), Exp("t", s$C)))) - log(Exp("x", lambda(s$d, shuffle(Exp("t", s$E)), Exp("t", s$C)))) 
curlyIBarPrime_T <- function(s)  Exp("t", log( lambda(s$d, Exp("x", s$E), Exp("x", s$C)))) - Exp("t", log( lambda(s$d, shuffle(Exp("x", s$E)), Exp("x",s$C))))
curlyIBarPrime_R <- function(s) curlyIBarPrime(s) - (curlyIBarPrime_A(s) + curlyIBarPrime_S(s) + curlyIBarPrime_T(s))


# calculate coexistence mechanisms, small-noise and exact -----------------

set.seed(20)
# parameters for three-species lottery model
K <- 500 # number of patches
num_spp <- 3
inv_index <- 1
inv_ceiling <- 0.001
n_mat <- matrix(c(0.0000001, 0.5, 0.5), nrow = K, ncol = num_spp, byrow = TRUE)
d_vec <- c(0.11, 0.1, 0.1)
dispersal_vec <- c(0.9,0.9,0.9)
d_mat <- sapply(seq_len(num_spp), function(i) rep(d_vec[i], K))
dispersal_mat <- sapply(seq_len(num_spp), function(i) rep(dispersal_vec[i], K))
record_time <- 1000
burn_in <- 200
total_time <- burn_in + record_time

env_sd <- 0.3
space_effect_mu <- 1
space_effect_sd <- env_sd
time_effect_mu <- 1 
time_effect_sd <- env_sd
st_interaction <- 1
extra_effect_mu <- 0 
extra_effect_sd <- 0
# make base_pre_E for each species
E_arr <- array(NA, dim = c(num_spp, K, total_time))
for(i in seq_len(num_spp))
{
  time_effect <- matrix(rep(rnorm(total_time, time_effect_mu, time_effect_sd), K), ncol = total_time, byrow = TRUE)
  space_effect <- matrix(rep(rnorm(K, space_effect_mu, space_effect_sd), total_time), nrow = K, byrow = FALSE)
  extra <- matrix(rnorm(length(time_effect), extra_effect_mu, extra_effect_sd), nrow = K)
  E_mat_1sp <- time_effect + space_effect + space_effect*time_effect*st_interaction + extra
  E_arr[i,,] <- E_mat_1sp
}
E_arr[which(E_arr < 0, arr.ind = TRUE)] <- 0 
spp_dat_list <- get_data(n_mat, burn_in, record_time)



for(i in seq_len(num_spp))
{
  E_det <- space_effect_mu + time_effect_mu + space_effect_mu*time_effect_mu*st_interaction + extra_effect_mu
  p <- list(
    E_star =  E_det,
    C_star = E_det - log(d_vec[i]),
    d = d_vec[i],
    qij = d_vec[inv_index]/d_vec[i],
    alphaSup1 = d_vec[i],
    betaSup1 = -d_vec[i],
    alphaSup2 = d_vec[i],
    betaSup2 = d_vec[i],
    zeta = -d_vec[i]
  )
  spp_dat_list[[i]] <- append(spp_dat_list[[i]], p)
}

str(spp_dat_list[[inv_index]]$n)

# small noise
rPrime <- partition(rPrimePre, spp_dat_list)
rPrime_A <- partition(rPrimePre_A, spp_dat_list)
rPrime_S <- partition(rPrimePre_S, spp_dat_list)
rPrime_T <- partition(rPrimePre_T, spp_dat_list)
rPrime_R <- partition(rPrimePre_R, spp_dat_list)

deltaRho <-  partition(deltaRhoPre, spp_dat_list)

deltaN_S <-  partition(deltaNPre_S, spp_dat_list)
deltaN_T <-  partition(deltaNPre_T, spp_dat_list)
deltaN_R <-  partition(deltaNPre_R, spp_dat_list)
deltaN <- deltaN_S + deltaN_T + deltaN_R

deltaI_A <- 0 
deltaI_S <-  partition(deltaIPre_S, spp_dat_list)
deltaI_T <-  partition(deltaIPre_T, spp_dat_list)
deltaI_R <-  partition(deltaIPre_R, spp_dat_list)
deltaI <- deltaI_S + deltaI_T + deltaI_R

deltaKappa <-  partition(deltaKappaPre, spp_dat_list)

inv_gr <- rPrime + deltaRho + deltaN + deltaI + deltaKappa


# exact
rPrimeExact <- partition(curlyEBarPrime, spp_dat_list)
rPrimeExact_A <- partition(curlyEBarPrime_A, spp_dat_list)
rPrimeExact_S <- partition(curlyEBarPrime_S, spp_dat_list)
rPrimeExact_T <- partition(curlyEBarPrime_T, spp_dat_list)
rPrimeExact_R <- partition(curlyEBarPrime_R, spp_dat_list)

deltaRhoExact <-  partition(curlyCBarPrime_A, spp_dat_list)

deltaNExact_S <-  partition(curlyCBarPrime_S, spp_dat_list)
deltaNExact_T <-  partition(curlyCBarPrime_T, spp_dat_list)
deltaNExact_R <-  partition(curlyCBarPrime_R, spp_dat_list)
deltaNExact <- deltaNExact_S + deltaNExact_T + deltaNExact_R

deltaIExact_A <-  partition(curlyIBarPrime_A, spp_dat_list)
deltaIExact_S <-  partition(curlyIBarPrime_S, spp_dat_list)
deltaIExact_T <-  partition(curlyIBarPrime_T, spp_dat_list)
deltaIExact_R <-  partition(curlyIBarPrime_R, spp_dat_list)
deltaIExact <-  partition(curlyIBarPrime, spp_dat_list)


si <- spp_dat_list[[inv_index]]
inv_gr_exact <- Exp("t",log(lambdaTwiddle(si$n, si$d, si$E, si$C)))

deltaKappaExact <-  inv_gr_exact - (rPrimeExact + deltaRhoExact + deltaNExact + deltaIExact)


# coef-plots of coexistence mechanisms: small nosie and exact -------------

# aggregate the values calculated in vectors. The zeros are to separate classes of coexistence mechanims in the plot
decomp_exact_values <- rev(c(inv_gr_exact, 0, rPrimeExact, rPrimeExact_A, rPrimeExact_S, rPrimeExact_T, rPrimeExact_R, 0, deltaRhoExact, 0, deltaNExact, deltaNExact_S, deltaNExact_T, deltaNExact_R, 0, deltaIExact, deltaIExact_A, deltaIExact_S, deltaIExact_T, deltaIExact_R, 0, deltaKappaExact))
decomp_values <- rev(c(inv_gr, 0, rPrime, rPrime_A, rPrime_S, rPrime_T, rPrime_R, 0, deltaRho, 0, deltaN, deltaN_S, deltaN_T, deltaN_R, 0, deltaI, deltaI_A, deltaI_S, deltaI_T, deltaI_R, 0, deltaKappa))

# define factors for each coexistence mechanims. these will be used to define colors. For example, frequency independent effects is a different factor than the space-time decomp components of the frequency independent effects.
col_groups <- as.factor(rev(c(1,0,2,3,3,3,3,0,4,0,5,6,6,6,0,7,8,8,8,8,0,9)))

# put it all in a data frame
inv_gr_decomp <- data.frame(names = seq_along(decomp_exact_values), values = decomp_values, exact_values = decomp_exact_values, col_groups = col_groups)

# define color pallete. One color for each factor level in col groups.
pal <- c("forestgreen", saturation("red3", 0.5), saturation("green3", 0.7), saturation("green3", 0.3), saturation("orange3", 0.9), saturation("purple3", 0.7), saturation("purple3", 0.4), saturation("blue2", 0.7), saturation("blue2", 0.4), "brown")

# define arguments for custom plot breaks and labelling
group_char <- as.character(rev(col_groups)) 
limits_inv <- 1:(length(group_char)) # ticks are associated with values 1 to the length of the data vector
breaks_inv <- limits_inv
breaks_inv[which(group_char == 0)] <- NA # ticks associated with the col_group factor = 0 (denoting no value) gets no tick break
# define_labels
labels_inv_small_noise <- rev(c(parse(text = TeX('$ \\approx E_{t} \\lbrack \\log ( \\widetilde{\\lambda_{i}} ) \\rbrack$')),
                                "",
                                parse(text = TeX("$\\Delta E_i$")),
                                parse(text = TeX("$\\Delta E_{i,A}$")),
                                parse(text = TeX("$\\Delta E_{i,S}$")),
                                parse(text = TeX("$\\Delta E_{i,T}$")),
                                parse(text = TeX("$\\Delta E_{i,R}$")),
                                "",
                                parse(text = TeX('$\\Delta \\rho_{i}$')),
                                "",
                                parse(text = TeX("$\\Delta N_{i}$")),
                                parse(text = TeX("$\\Delta N_{i,S}$")),
                                parse(text = TeX("$\\Delta N_{i,T}$")),
                                parse(text = TeX("$\\Delta N_{i,R}$")),
                                "",
                                parse(text = TeX("$\\Delta I_{i}$")),
                                parse(text = TeX("$\\Delta I_{i,A}$")),
                                parse(text = TeX("$\\Delta I_{i,S}$")),
                                parse(text = TeX("$\\Delta I_{i,T}$")),
                                parse(text = TeX("$\\Delta I_{i,R}$")),
                                "",
                                parse(text = TeX('$\\Delta \\kappa_i$'))))

#text = TeX('$E_{t} \\lbrack \\log ( \\widetilde{\\lambda_{i}} ) \\rbrack$')
#text = TeX('$r_i$')
labels_inv_exact <- rev(c(parse(text = TeX('$E_{t} \\lbrack \\log ( \\widetilde{\\lambda_{i}} ) \\rbrack$')),
                          "",
                          parse(text = TeX("$\\Delta E_i^{(e)}$")),
                          parse(text = TeX("$\\Delta E_{i,A}^{(e)}$")),
                          parse(text = TeX("$\\Delta E_{i,S}^{(e)}$")),
                          parse(text = TeX("$\\Delta E_{i,T}^{(e)}$")),
                          parse(text = TeX("$\\Delta E_{i,R}^{(e)}$")),
                          "",
                          parse(text = TeX('$\\Delta \\rho_{i}^{(e)}$')),
                          "",
                          parse(text = TeX("$\\Delta N_{i}^{(e)}$")),
                          parse(text = TeX("$\\Delta N_{i,S}^{(e)}$")),
                          parse(text = TeX("$\\Delta N_{i,T}^{(e)}$")),
                          parse(text = TeX("$\\Delta N_{i,R}^{(e)}$")),
                          "",
                          parse(text = TeX("$\\Delta I_{i}^{(e)}$")),
                          parse(text = TeX("$\\Delta I_{i,A}^{(e)}$")),
                          parse(text = TeX("$\\Delta I_{i,S}^{(e)}$")),
                          parse(text = TeX("$\\Delta I_{i,T}^{(e)}$")),
                          parse(text = TeX("$\\Delta I_{i,R}^{(e)}$")),
                          "",
                          parse(text = TeX('$\\Delta \\kappa_i^{(e)}$'))))

#font_import(paths = "C:/Windows/my_fonts/", recursive = FALSE)
loadfonts(device="win")       #Register fonts for Windows bitmap output
#loadfonts()
#fonts()      
# plot exact coexistence mechanisms
inv_gr_decomp %>% 
  ggplot(aes(names, exact_values, fill = col_groups)) +
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = pal) +
  coord_flip() + 
  scale_x_discrete( limits = limits_inv,
                    breaks = breaks_inv,
                    labels = labels_inv_exact) + 
  theme_classic() +
  xlab("Exact coexistence mechanism") + 
  ylab("Growth rate") + 
  theme( axis.text.x = element_text(size= 24, family="LM Roman 10"), axis.text.y = element_text(size= 18, family="LM Roman 10"),
         axis.title.x = element_text(size= 24, family="LM Roman 10"), axis.title.y = element_text(size= 24, family="LM Roman 10")) +
  geom_hline(yintercept = 0, lty = 1, lwd = 1)


# plot small-noise coexistence mechanisms
inv_gr_decomp %>% 
  ggplot(aes(names, values, fill = col_groups)) +
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = pal) +
  coord_flip() + 
  scale_x_discrete( limits = limits_inv,
                    breaks = breaks_inv,
                    labels = labels_inv_small_noise) + 
  xlab("small-noise coexistence mechanism") + 
  ylab("growth rate") + 
  theme( axis.text.x = element_text(size= 16), axis.text.y = element_text(size= 16),
         axis.title.x = element_text(size= 16), axis.title.y = element_text(size= 16)) +
  geom_hline(yintercept = 0, lty = 1, lwd = 1)


# make figure for paper
png("lottery_model.png", units="in", width=14, height=10, res=300)
inv_gr_decomp %>% 
  ggplot(aes(names, exact_values, fill = col_groups)) +
  geom_col(show.legend = FALSE) + 
  scale_fill_manual(values = pal) +
  coord_flip() + 
  scale_x_discrete( limits = limits_inv,
                    breaks = breaks_inv,
                    labels = labels_inv_exact) + 
  theme_classic() +
  xlab("Exact coexistence mechanism") + 
  ylab("Growth rate") + 
  theme( axis.text.x = element_text(size= 24, family="LM Roman 10"), axis.text.y = element_text(size= 18, family="LM Roman 10"),
         axis.title.x = element_text(size= 24, family="LM Roman 10"), axis.title.y = element_text(size= 24, family="LM Roman 10")) +
  geom_hline(yintercept = 0, lty = 1, lwd = 1)
dev.off()




# plot small-noise and exact storage effect against environmental noise -------------------

# check if they converge in the limit of small noise


set.seed(20)
# parameters for three-species lottery model
K <- 500 # number of patches
num_spp <- 3
inv_index <- 1
inv_ceiling <- 0.001
n_mat <- matrix(c(0.0000001, 0.5, 0.5), nrow = K, ncol = num_spp, byrow = TRUE)
d_vec <- c(0.11, 0.1, 0.1)
dispersal_vec <- c(0.9,0.9,0.9)
d_mat <- sapply(seq_len(num_spp), function(i) rep(d_vec[i], K))
dispersal_mat <- sapply(seq_len(num_spp), function(i) rep(dispersal_vec[i], K))
record_time <- 500
burn_in <- 20
total_time <- burn_in + record_time

# define a vector of standard deviations of the environment. This sd is the sd of space effects, time effects, and random effects. 
env_sd <- seq(from = 0, to = 1, length.out = 20)
deltaI_vec <- numeric(length(env_sd))
deltaIExact_vec <- numeric(length(env_sd))

for(k in seq_along(env_sd))
{
  
  space_effect_mu <- 1
  space_effect_sd <- env_sd[k]
  time_effect_mu <- 1
  time_effect_sd <- env_sd[k]
  st_interaction <- 0.1
  extra_effect_mu <- 0 
  extra_effect_sd <- 0
  
  # make the E ahead of time, for each species
  E_arr <- array(NA, dim = c(num_spp, K, total_time))
  for(i in seq_len(num_spp))
  {
    time_effect <- matrix(rep(rnorm(total_time, time_effect_mu, time_effect_sd), K), ncol = total_time, byrow = TRUE)
    space_effect <- matrix(rep(rnorm(K, space_effect_mu, space_effect_sd), total_time), nrow = K, byrow = FALSE)
    extra <- matrix(rnorm(length(time_effect), extra_effect_mu, extra_effect_sd), nrow = K)
    E_mat_1sp <- time_effect + space_effect + space_effect*time_effect*st_interaction + extra
    E_arr[i,,] <- E_mat_1sp
  }
  E_arr[which(E_arr < 0, arr.ind = TRUE)] <- 0 
  
  # get data
  spp_dat_list <- get_data(n_mat, burn_in, record_time)
  # associate data and parameters for each each species
  for(i in seq_len(num_spp))
  {
    E_det <- space_effect_mu + time_effect_mu + space_effect_mu*time_effect_mu*st_interaction + extra_effect_mu
    p <- list(
      E_star =  E_det,
      C_star = E_det - log(d_vec[i]),
      d = d_vec[i],
      qij = d_vec[inv_index]/d_vec[i],
      alphaSup1 = d_vec[i],
      betaSup1 = -d_vec[i],
      alphaSup2 = d_vec[i],
      betaSup2 = d_vec[i],
      zeta = -d_vec[i]
    )
    spp_dat_list[[i]] <- append(spp_dat_list[[i]], p)
  }
  
  # compute exact storage effect
  
  deltaIExact_vec[k] <- partition(curlyIBarPrime, spp_dat_list)
  deltaI_vec[k] <- partition(deltaIPre, spp_dat_list)
}

plot(x = env_sd, y =deltaIExact_vec, type = "b")
points(x = env_sd, y =deltaI_vec, type = "b")

data.frame(env_sd, deltaI_vec, deltaIExact_vec) %>%
  rename(Exact = deltaIExact_vec, `Small-noise` = deltaI_vec) %>%
  gather(key = "key", value = "value", -env_sd) %>%
  ggplot(aes(env_sd, value, color = key)) + 
  geom_point(size = 3) + 
  geom_line() +
  theme_bw() + 
  xlab("standard deviation of variation in the environment") + 
  ylab("storage effect") + 
  theme(axis.title = element_text(size = 16))

# plot the space component of the exact storage effect against spatial environmental noise -------------------

# check if they converge in the limit of small noise

set.seed(20)
# parameters for three-species lottery model
K <- 500 # number of patches
num_spp <- 3
inv_index <- 1
inv_ceiling <- 0.001
n_mat <- matrix(c(0.0000001, 0.5, 0.5), nrow = K, ncol = num_spp, byrow = TRUE)
d_vec <- c(0.11, 0.1, 0.1)
dispersal_vec <- c(0.9,0.9,0.9)
d_mat <- sapply(seq_len(num_spp), function(i) rep(d_vec[i], K))
dispersal_mat <- sapply(seq_len(num_spp), function(i) rep(dispersal_vec[i], K))
record_time <- 500
burn_in <- 20
total_time <- burn_in + record_time

# define a vector of standard deviations of the environment. This sd is the sd of space effects, time effects, and random effects. 
space_sd <- seq(from = 0, to = 1, length.out = 20)
deltaIExact_S_vec <- numeric(length(space_sd))

for(k in seq_along(env_sd))
{
  
  space_effect_mu <- 3
  space_effect_sd <- space_sd[k]
  time_effect_mu <- 3
  time_effect_sd <- 1
  st_interaction <- .1
  extra_effect_mu <- 0 
  extra_effect_sd <- 0
  
  # make the E ahead of time, for each species
  E_arr <- array(NA, dim = c(num_spp, K, total_time))
  for(i in seq_len(num_spp))
  {
    time_effect <- matrix(rep(rnorm(total_time, time_effect_mu, time_effect_sd), K), ncol = total_time, byrow = TRUE)
    space_effect <- matrix(rep(rnorm(K, space_effect_mu, space_effect_sd), total_time), nrow = K, byrow = FALSE)
    extra <- matrix(rnorm(length(time_effect), extra_effect_mu, extra_effect_sd), nrow = K)
    E_mat_1sp <- time_effect + space_effect + space_effect*time_effect*st_interaction + extra
    E_arr[i,,] <- E_mat_1sp
  }
  E_arr[which(E_arr < 0, arr.ind = TRUE)] <- 0 
  
  # get data
  spp_dat_list <- get_data(n_mat, burn_in, record_time)
  # associate data and parameters for each each species
  for(i in seq_len(num_spp))
  {
    E_det <- space_effect_mu + time_effect_mu + space_effect_mu*time_effect_mu*st_interaction + extra_effect_mu
    p <- list(
      E_star =  E_det,
      C_star = E_det - log(d_vec[i]),
      d = d_vec[i],
      qij = d_vec[inv_index]/d_vec[i],
      alphaSup1 = d_vec[i],
      betaSup1 = -d_vec[i],
      alphaSup2 = d_vec[i],
      betaSup2 = d_vec[i],
      zeta = -d_vec[i]
    )
    spp_dat_list[[i]] <- append(spp_dat_list[[i]], p)
  }
  
  # compute exact storage effect
  
  deltaIExact_S_vec[k] <- partition(curlyIBarPrime_S, spp_dat_list)
}


data.frame(space_sd, deltaIExact_S_vec) %>%
  ggplot(aes(space_sd, deltaIExact_S_vec)) + 
  geom_point(size = 3) + 
  geom_line() +
  theme_bw() + 
  xlab("standard deviation of spatial variation in the environment") + 
  ylab("space component of the exact storage effect") + 
  theme(axis.title = element_text(size = 16))


# plot and small-noise and exact space components of the storage effect against environmental noise -------------------

# check if they converge in the limit of small noise

set.seed(20)
# parameters for three-species lottery model
K <- 500 # number of patches
num_spp <- 3
inv_index <- 1
inv_ceiling <- 0.001
n_mat <- matrix(c(0.0000001, 0.5, 0.5), nrow = K, ncol = num_spp, byrow = TRUE)
d_vec <- c(0.11, 0.1, 0.1)
dispersal_vec <- c(0.9,0.9,0.9)
d_mat <- sapply(seq_len(num_spp), function(i) rep(d_vec[i], K))
dispersal_mat <- sapply(seq_len(num_spp), function(i) rep(dispersal_vec[i], K))
record_time <- 500
burn_in <- 100
total_time <- burn_in + record_time

# define a vector of standard deviations of the environment. This sd is the sd of space effects, time effects, and random effects. 
#space_sd <- seq(from = 0, to = 1, length.out = 20)
space_sd <- seq(from = 0, to = 1, length.out = 20)
deltaIExact_S_vec <- numeric(length(space_sd))
deltaI_S_vec <- numeric(length(space_sd))

for(k in seq_along(env_sd))
{
  
  space_effect_mu <- 3
  space_effect_sd <- space_sd[k]
  time_effect_mu <- 3
  time_effect_sd <- 1
  st_interaction <- .1
  extra_effect_mu <- 0 
  extra_effect_sd <- 0
  
  # make the E ahead of time, for each species
  E_arr <- array(NA, dim = c(num_spp, K, total_time))
  for(i in seq_len(num_spp))
  {
    time_effect <- matrix(rep(rnorm(total_time, time_effect_mu, time_effect_sd), K), ncol = total_time, byrow = TRUE)
    space_effect <- matrix(rep(rnorm(K, space_effect_mu, space_effect_sd), total_time), nrow = K, byrow = FALSE)
    extra <- matrix(rnorm(length(time_effect), extra_effect_mu, extra_effect_sd), nrow = K)
    E_mat_1sp <- time_effect + space_effect + space_effect*time_effect*st_interaction + extra
    E_arr[i,,] <- E_mat_1sp
  }
  E_arr[which(E_arr < 0, arr.ind = TRUE)] <- 0 
  
  # get data
  spp_dat_list <- get_data(n_mat, burn_in, record_time)
  # associate data and parameters for each each species
  for(i in seq_len(num_spp))
  {
    E_det <- space_effect_mu + time_effect_mu + space_effect_mu*time_effect_mu*st_interaction + extra_effect_mu
    p <- list(
      E_star =  E_det,
      C_star = E_det - log(d_vec[i]),
      d = d_vec[i],
      qij = d_vec[inv_index]/d_vec[i],
      alphaSup1 = d_vec[i],
      betaSup1 = -d_vec[i],
      alphaSup2 = d_vec[i],
      betaSup2 = d_vec[i],
      zeta = -d_vec[i]
    )
    spp_dat_list[[i]] <- append(spp_dat_list[[i]], p)
  }
  #curlyIBarPrime_S <- function(s) log(Exp("x", lambda(s$d, Exp("t", s$E), Exp("t", s$C)))) - log(Exp("x", lambda(s$d, shuffle(Exp("t", s$E)), Exp("t", s$C)))) 
  curlyIBarPrime_S <- function(s) log(Exp("x", lambda(s$d, Exp("t", s$E), Exp("t", s$C)))) - log(Exp("x", lambda(s$d, shuffle(Exp("t", s$E)), Exp("t", s$C)))) - curlyIBarPrime_A(s)
  
  # compute exact storage effect
  deltaIExact_S_vec[k] <- partition(curlyIBarPrime_S, spp_dat_list)
  # compute small-noise storage effect
  deltaI_S_vec[k] <-  partition(deltaIPre_S, spp_dat_list)
  
}


data.frame(space_sd, deltaIExact_S_vec, deltaI_S_vec) %>%
  gather(key = "key", value = "gr", deltaIExact_S_vec, deltaI_S_vec) %>%
  ggplot(aes(x = space_sd, y = gr, color = key)) + 
  geom_point(size = 3) + 
  geom_line() +
  theme_bw() + 
  xlab("standard deviation of spatial variation in the environment") + 
  ylab("space component of the exact storage effect") + 
  theme(axis.title = element_text(size = 16))

