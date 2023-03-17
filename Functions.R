# Define Matern covariance function
matern.cov<-function(d, nu = 2, kappa = 1,sig=1){
  if(nu == 0.5){
    return(sig^2*exp( - d * alpha))
  }
  gamma.d<- ifelse(d>0,sig^2 * (2^(1-nu)/ gamma(nu)) * ((d*kappa)^nu) * besselK(d*kappa, nu),sig^2)
  return(gamma.d)
}

# Function for simulating the data
GP.data <- function(coordinates,nu=1,kappa=1,sig=1,mu=50){
  coordinates = as.data.frame(coordinates)
  n <- nrow(coordinates)
  D <- dist(coordinates)
  sim_mat <- matern.cov(as.matrix(D),nu=nu,kappa=kappa,sig=sig)
  data <- mvrnorm(mu=rep(mu,n), Sigma=sim_mat)
  return(data)
}

GP.data(s.grid, nu=2, kappa=4, sig=50, mu=50)

# Function for testing the mesh
## including the data, locations, boundary,
## the parameters of mesh, and the parameters of SPDE
mesh.test <- function(data, coordinates, boundary, max_edge_1,
                      max_edge_2, min_angle_1, min_angle_2,
                      cutoff, offset_1, offset_2,
                      smoothness_alpha, range_1, range_2,
                      sigma_1, sigma_2){
  mesh0 <- inla.mesh.2d(loc = coordinates, boundary = boundary,
                      max.edge = c(max_edge_1,max_edge_2),
                      min.angle = c(min_angle_1,min_angle_2),
                      cutoff=cutoff,
                      offset=c(offset_1,offset_2))
  
  spde0 <- inla.spde2.pcmatern(
                      mesh = mesh0, alpha = smoothness_alpha,
                      prior.range = c(range_1,range_2),
                      prior.sigma = c(sigma_1,sigma_2))
  
  A <- inla.spde.make.A(mesh0, loc = coordinates)
  
  dat.stk <- inla.stack(
              data = list(resp = data),
              A = list(A, 1),
              effects = list(i = 1:spde0$n.spde,
                             beta0 = rep(1, length(data))))
  
  res <- inla(resp ~ 0 + beta0 + f(i, model = spde0),
            data = inla.stack.data(dat.stk),
            control.predictor = list(A = inla.stack.A(dat.stk)))
  
  sigma.hat = inla.emarginal(function(x) x, res$marginals.hyperpar[[2]])
  range.hat = inla.emarginal(function(x) x, res$marginals.hyperpar[[3]])
  
  hyperparameters = list(sigma_hat = sigma.hat, range_hat = range.hat)
  
  return(hyperparameters)
}

# Function for testing the mesh (2)
## including the data, locations, boundary,
## The priors of SPDE are compute automatically
mesh.test2 <- function(data, coordinates, boundary, max_edge_1,
                       max_edge_2, min_angle_1, min_angle_2,
                       cutoff, offset_1, offset_2){
  
  df <- as.data.frame(coordinates)
  df$dat = data
  colnames(df) <- c('lon','lat','dat')
  
  mesh0 <- inla.mesh.2d(loc = coordinates, boundary = boundary,
                        max.edge = c(max_edge_1,max_edge_2),
                        min.angle = c(min_angle_1,min_angle_2),
                        cutoff=cutoff,
                        offset=c(offset_1,offset_2))
  
  ## priors for SPDE
  ## sigma
  sigma_0 <- sd(df$dat)+0.1
  
  ## range
  vgm.emp <- variogram(dat~1, data = df,
                       locations = ~lon+lat)
  fitvariogram<-fit.variogram(vgm.emp, vgm(c( "Mat")), 
                              fit.kappa = TRUE)
  range_0 <- fitvariogram$range[2]
  
  spde0 <- inla.spde2.pcmatern(
    mesh = mesh0, alpha = 2,
    prior.range = c(range_0,0.85),
    prior.sigma = c(sigma_0,0.01))
  
  A <- inla.spde.make.A(mesh0, loc = coordinates)
  
  dat.stk <- inla.stack(
    data = list(resp = df$dat),
    A = list(A, 1),
    effects = list(i = 1:spde0$n.spde,
                   beta0 = rep(1, nrow(df))))
  
  res <- inla(resp ~ 0 + beta0 + f(i, model = spde0),
              data = inla.stack.data(dat.stk),
              control.predictor = list(A = inla.stack.A(dat.stk)))
  
  range.hat = inla.emarginal(function(x) x, res$marginals.hyperpar[[2]])
  sigma.hat = inla.emarginal(function(x) x, res$marginals.hyperpar[[3]])
  
  hyperparameters = list(range0 = range_0, sigma0 = sigma_0,
                         range_hat = range.hat, sigma_hat = sigma.hat)
  
  return(hyperparameters)
}





# Function for finding the best mesh

lost.func <- function(range_true, sigma_true, range_hat, sigma_hat){
  bias_range <- abs(range_hat - range_true)
  bias_sigma <- abs(sigma_hat - sigma_true)
  
}



















