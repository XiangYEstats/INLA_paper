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
# Function include constructing a mesh, fit an SPDE, and do the INLA to find the posteriors of the parameters
mesh.test2 <- function(data, coordinates, boundary, max_edge_1,
                       max_edge_2, min_angle_1, min_angle_2,
                       cutoff, offset_1, offset_2){
  
  ## data: the data
  ## coordinates: the longitude and latitude of the data point (recommended in kilometers)
  ## boundary: boundary of the mesh
  
  ## max_edge_1: maximum allowed triangle edge length in the inner domain
  ## max_edge_2: maximum allowed triangle edge length in the outer extension
  
  ## min_angle_1: the minimum internal angles of the triangles in the inner domain
  ## min_angle_2: the minimum internal angles of the triangles in the outer extension
  ### min_angle values up to 21 guarantee the convergence of the algorithm
  
  ## cutoff: to set the minimum allowed distance between points
  
  ## offset_1: set the automatic extension distance in direction 1
  ## offset_2: set the automatic extension distance in direction 2 
  
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
              control.predictor = list(A = inla.stack.A(dat.stk)),
              control.compute = list(dic = TRUE))
  
  ## expectation of the posterior standard error
  post.se <- inla.tmarginal(function(x) sqrt(1 / exp(x)),
                            res$internal.marginals.hyperpar[[1]])
  post.se.e <- inla.emarginal(function(x) x, post.se)
  
  precision = inla.emarginal(function(x) x, res$marginals.hyperpar[[1]])
  range.hat = inla.emarginal(function(x) x, res$marginals.hyperpar[[2]])
  sigma.hat = inla.emarginal(function(x) x, res$marginals.hyperpar[[3]])
  dic = res$dic$dic
  
  hyperparameters = list(range0 = range_0, sigma0 = sigma_0, range_hat = range.hat, sigma_hat = sigma.hat, dic = dic, post_expect_se = post.se.e, precision = precision)
  
  return(hyperparameters)
}





# Function for finding the best mesh

lost.func <- function(range_true, sigma_true, range_hat, sigma_hat){
  bias_range <- abs(range_hat - range_true)
  bias_sigma <- abs(sigma_hat - sigma_true)
  
}



















