#simulate data

library(TMB)

# set.seed(123) 

n_locations <- 3 #Number of survey locations
n_individuals <- 10000     # Number of time points
capture <- matrix(NA , n_individuals , n_locations) #Capture history matrix
z <- matrix(NA , n_individuals , n_locations) #Capture history matrix

time <- sample(0:5,n_individuals,replace = TRUE) # Time points

# names(capture) <- c("Wen","MCN","DNSTRM") 
logit_phi <- qlogis(rep(0.6, n_locations)) # Logit-transformed capture probabilities
logit_p <- qlogis(rep(0.1, n_locations)) # Logit-transformed detection probabilities
lam <- 0     # Log lambda for time-density model
  


# Loop over individuals
for (i in 1:nrow(capture)) {
    n_captures <- ncol(capture);  # Number of capture occasions for the individual
    # Loop over survey locations
    eta_phi = logit_phi[1] + lam * time[i]
    phi = plogis(eta_phi)
    
    z[i,1] <- rbinom(1,1,phi)
    capture[i,1] <- rbinom(1,z[i,1],p)
    
    for (j in 2:n_locations) {
      # Loop over capture occasions
      # Probability of capture at location j and time k

      # Probability of detection at location j
      p = plogis(logit_p)[1];
      # Calculate the log-likelihood contribution for this capture occasion
      z[i,j] <- rbinom(1,z[i,j-1],phi)
      capture[i,j] <- rbinom(1,z[i,j],p)
    }
}

#Get first caption "location"
f <- IPMbook::getFirst(as.matrix(capture))

#Remove capture where f is equal to last location
time <- na.omit(time[f < ncol(capture)])
capture <- capture[ f < ncol(capture),]
capture <- na.omit(as.data.frame(capture))
f <- na.omit(f[f<ncol(capture)])

#Get hidden process
z <- IPMbook::zKnown(capture)
z[is.na(z)] <- 0

try(dyn.unload("code/wenatchee.dll"))

compile("code/wenatchee.cpp")
dyn.load(dynlib("code/wenatchee"))
# # 
data <- list(c_it = as.matrix(capture), 
             z_it = as.matrix(z), 
             f_i = as.vector(f), 
             t_i = as.vector(time))

parameters <-list(f_phi = rep(1,1), f_p = rep(1,1), lam = 1)

obj <- MakeADFun(data, 
                 parameters,
                 map = list(lam = as.factor(NA)),
                 DLL = "wenatchee")
obj$hessian <- TRUE
opt<-nlminb(obj$par,obj$fn,obj$gr)

print(round(c(logit_phi[1],logit_p[1],lam),2))
print(round(opt$par,2))

# paste(round(plogis(logit_phi[1]),2),round(plogis(logit_p[1]),2))# opt$hessian ## <-- FD hessian from optim
# obj$he()    ## <-- Analytical hessian
# sdreport(obj)
