require(tidyverse)
require(issb)
library(cowplot)

# Modify issb functions to return the drift path -----------------------------------

lnafun = function (t, state, model) 
{
      s = model$get_stoic()
      u = nrow(s)
      z0 = state[1:u]
      m0 = state[(u + 1):(2 * u)]
      V0 = matrix(state[(2 * u + 1):length(state)], ncol = u)
      Ft = model$get_jacobian(z0)
      h = model$get_haz(z0)
      z = s %*% ((exp(-z0) - 0.5*exp(-2*z0))*h)
      m = s %*% Ft %*% m0
      V = V0 %*% t(Ft) %*% t(s) + s %*% diag(exp(-2*z0) * h) %*% t(s) + s %*% 
            Ft %*% V0
      return(list(c(z, m, as.vector(V))))
}

lnastep = function(model, z, m, t, ddt)
{
      u = length(z)
      state = c(z, m, rep(0, u^2))	
      sol = deSolve:::ode(y=state, func=lnafun, times=seq(0, t, by = ddt), parms=model)
      
      #Remove first column and first row
      sol = sol[2:nrow(sol), -1]
      z = sol[nrow(sol), 1:u]
      m = sol[nrow(sol), (u+1):(2*u)]
      V = matrix(sol[nrow(sol), -(1:(2*u))], ncol=u)
      
      # get the path
      drift_path = sol[,1:u]
      
      return(list(z=z, m=m, V=V, drift_path = drift_path))
}

lna = function(model, maxtime, ddt, restart=FALSE) 
{
      
      z = model$get_initial()
      N = maxtime/ddt + 1
      ind_seq = seq_along(seq(ddt,1,by=ddt))
      
      # for N(0,1) draws
      draws = matrix(0, nrow = maxtime, ncol = length(z))
      
      # for the drift path
      pathmat = matrix(0, nrow = N, ncol = length(z))
      pathmat[1,] = z
      
      # for the centered sample
      xmat = matrix(0, nrow=maxtime+1, ncol=length(z))
      xmat[1,] = z
      
      # for the covariance matrix
      mean_vecs = matrix(0, ncol = length(z), nrow = maxtime)
      diff_mats = array(0, dim = c(length(z), length(z), maxtime))
      
      ##Initialise m and V
      m = z*0
      V = matrix(0, length(z), length(z))
      
      for (i in 2:nrow(xmat)){
            lsol = lnastep(model, z, m, 1, ddt)
            z = lsol[["z"]]
            e = eigen(lsol[["V"]])
            pathmat[1 + ind_seq + length(ind_seq) * (i-2),] = lsol[["drift_path"]]
            mean_vecs[i-1,]  = z
            diff_mats[,,i-1] = lsol[["V"]]
            
            draws[i-1,] = rnorm(length(z))
            xmat[i,] = z + lsol[["m"]] + 
                  ifelse(draws[i-1,]>0,
                         abs(e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% draws[i-1,]),
                         -abs(abs(e$vectors %*% diag(sqrt(pmax(e$values, 0))) %*% draws[i-1,])))

            # Reflecting barrier
            neg = xmat[i,] < 0
            xmat[i,][neg] = -xmat[i,][neg]
            
            over = xmat[i,] > model$get_pars()["popsize"]
            xmat[i,][over] = 2*model$get_pars() - xmat[i,][over]
            
            ##Reinitialise
            m = xmat[i, ] - lsol[["z"]] #set m
            if(restart){z = xmat[i,]; m = lsol[["m"]]*0}
      }
      pathmat = cbind(seq(0, maxtime, ddt), pathmat)
      xmat    = cbind(seq(0, nrow(xmat)-1), xmat)
      colnames(pathmat) = colnames(xmat) = c("Time", rownames(model$get_stoic()))
      
      return(list(pathmat = pathmat, 
                  xmat = xmat,
                  draws = draws,
                  mean_vecs = mean_vecs,
                  diff_mats = diff_mats))
}

# Set up SIR model --------------------------------------------------------
# Hazards
h = function(x, pars) {
      y = exp(x) - 1
      hazs = numeric(length(x))
      hazs[1] = pars[1] * (pars[3] - y[1]) * (pars[4] + y[1] - y[2])   # beta * S * I
      hazs[2] = pars[2] * (pars[4] + y[1] - y[2])        # mu * I
      return(hazs)
}

# Stoichiometry matrix
smat = diag(1,2)
rownames(smat) = c("N_SI", "N_IR")

# initial values
initial = c(0, 0)

# parameters
pars = c(beta = 1.5e-3, mu = 0.5, S0 = 1e3 - 7, I0 = 7, popsize = 1e3)

# jacobian
f = get_f = function(x, pars)
{
      exp_Z      = exp(x)
      exp_neg_Z  = exp(-x)
      exp_neg_2Z = exp(-2*x)
      expm1_Z    = exp(x)-1
      
      fmat = matrix(0, nrow=2, ncol=2)
      
      fmat[1,1] = (exp_neg_Z[1]-0.5*exp_neg_2Z[1])*(pars[1]*exp_Z[1]*(pars[3]-expm1_Z[1])-(pars[1]*(pars[4]+expm1_Z[1]-expm1_Z[2]))*exp_Z[1])-(exp_neg_Z[1]-0.5*(exp_neg_2Z[1]*2))*((pars[1]*(pars[4]+expm1_Z[1]-expm1_Z[2]))*(pars[3]-expm1_Z[1]))
      fmat[1,2] = -((exp_neg_Z[1]-0.5*exp_neg_2Z[1])*(pars[1]*exp_Z[2]*(pars[3]-expm1_Z[1])))
      fmat[2,1] = (exp_neg_Z[2]-0.5*exp_neg_2Z[2])*((pars[2])*exp_Z[1])
      fmat[2,2] = -((exp_neg_Z[2]-0.5*exp_neg_2Z[2])*((pars[2])*exp_Z[2])+(exp_neg_Z[2]-0.5*(exp_neg_2Z[2]*2))*((pars[2])*(pars[4]+expm1_Z[1]-expm1_Z[2])))
      return(fmat)
}

# create model
model = create_model(smat, h, initial, pars, f)

# simulate a path
set.seed(52787)
path = lna(model = model, maxtime = 10, ddt = 0.001, restart = TRUE)

# Make pretty plots -------------------------------------------------------

gaussians = data.frame(loc_x = rep(seq(1,10), each = length(seq(-4,4,by=0.01))),
                       loc_y = rep(seq(-4,4,by=0.01), each = length(seq(1,10))),
                       x = rep(seq(-4,4,by=0.01)), 
                       dens = dnorm(seq(-4,4,by=0.01)))

noncentered_path = data.frame(x = c(seq(1,10)),
                              loc_x = c(path$draws[,1]))

centered_path = data.frame(x = seq(1,10),
                           loc_x = path$xmat[-1,2],
                           interval = 1:10)

yseqs <- 
      rep(path$pathmat[path$pathmat[,1] %%1 == 0,2][-1], each = length(seq(-4,4,by=0.001))) + 
      rep(seq(-4,4,by=0.001), 10) * rep(sqrt(path$diff_mats[1,1,]), each = length(seq(-4,4,by=0.001)))

lna_dists = data.frame(loc_x = rep(seq(1,10), each = length(seq(-4,4,by=0.001))),
                       x = yseqs,
                       dens =
                             (dnorm(yseqs,
                                    rep(path$mean_vecs[,1], each = length(seq(-4,4,by=0.001))),
                                    rep(sqrt(path$diff_mats[1,1,]), each = length(seq(-4,4,by=0.001)))))*rep(seq(1,0.1,length=10)^3, each = length(seq(-4,4,by=0.001))))
                             # dnorm(seq(-4,4,by=0.001)))

lna_path = as.data.frame(path$pathmat)
lna_path$interval <- c(1,rep(1:10,each = 1000))

# plot
ggplot(lna_path, aes(x = Time, y = N_SI, group = as.factor(interval))) + 
      geom_line(linetype = 1) + 
      geom_line(data = lna_dists, aes(x = loc_x + dens, y = x, group = loc_x), alpha = 0.25) + 
      geom_point(data = centered_path, aes(x = x, y = loc_x), colour = "darkred", size = 3) + 
      theme_bw() + 
      labs(x = "Time", 
           y = bquote(paste(widetilde(bold(N))[SI](bold(theta)))), 
           title = bquote(paste("Centered parameterization: ", widetilde(bold(N)) ~"|"~ bold(theta) %~%" MVN"(bold(mu)(bold(theta)), bold(Sigma)(bold(theta)))))) + 
      scale_x_continuous(breaks = c(0,2,4,6,8,10), limits = c(0, 10.5)) + 
      theme(text = element_text(family = "Garamond", size = 40)) -> lna_N

ggplot(gaussians, aes(x = loc_x + dens, y = x, group = loc_x))+
      geom_line(alpha = 0.25) + 
      geom_point(data = noncentered_path, aes(x = x, y = loc_x), colour = "darkred", size = 3) + 
      theme_bw() + 
      labs(x = "Time", 
           y = bquote(paste(bold(Z)[SI], symbol("\136"), bold(theta))), 
           title = bquote(paste("Non-centered parameterization: ", doLNA(bold(Z)~","~theta) %=>% widetilde(bold(N))))) + 
      scale_x_continuous(limits = c(0,10.5),breaks = c(0,2,4,6,8,10), labels = c(0,2,4,6,8,10)) +
      geom_abline(slope = 0, intercept = 0, linetype = 2) + 
      theme(text = element_text(family = "Garamond", size = 40)) -> lna_Z

cowplot::plot_grid(lna_N, lna_Z, ncol = 1, rel_heights = c(2,1), align = "hv") -> lna_diagram

save_plot("lna_sampling_diagram.svg", lna_diagram, base_height = 12, base_width = 20)
