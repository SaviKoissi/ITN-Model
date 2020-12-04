## Load deSolve package
library(deSolve)

#Parameters 
#Human hosts
Lambda_h<-matrix(c(0.005369485,0.000001,0.000001,0.000001,0.005369485,0.000001,0.000001,0.000001), 2, 4,byrow=T, dimnames = list(c("WU","Sl"),c("A1", "A2","A3","A4"))) #Recruitment rate of humans
#Lambda_h<-matrix(c(0.005369485,0.005369485,0.005369485,0.005369485,0.005369485,0.005369485,0.005369485,0.005369485), 2, 4,byrow=T, dimnames = list(c("WU","Sl"),c("A1", "A2","A3","A4"))) #Recruitment rate of humans
#mu_h<-matrix(c(1/63.78,1/63.78,1/63.78,1/63.78,1/63.78,1/63.78,1/63.78,1/63.78)/365, 2, 4)#Death rate of human
mu_h<-matrix(c(1/63.78,1/63.78,1/63.78,1/63.78,1/63.48,1/63.48,1/63.48,1/63.48)/365, 2, 4)#Death rate of human
gamma_h<-matrix(c(0.0023,0.0023,0.0023,0.0023,0.0023,0.0023,0.0023,0.0023), 2, 4,byrow=T)#Recovery rate
sigma_h<-matrix(c(1/91.3125,1/91.3125,1/91.3125,1/91.3125,1/91.3125,1/91.3125,1/91.3125,1/91.3125), 2, 4,byrow=T)#Proportion of getting immune 
beta_h<-matrix(c(0.42,0.42,0.42,0.42,0.42,0.42,0.42,0.42),2, 4,byrow=T)#Transmission rate from infectious human to mosquito
delta_h<-matrix(c(0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001),2, 4,byrow=T)#Disease induced-death
m<-matrix(c(1,1,1,1,1,1,1,1)*1e-10,nrow = 2, ncol= 4,byrow=T)#Between patches migration
psi<-matrix(c(0.6,0.6,0.6,0.6,0.6,0.3,0.3,0.3), 2, 4,byrow=T)#Proportion of ITN use 

#Mosquito vectors
mu_v<-matrix(c(0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05),2, 4,byrow=T)#Death rate of mosquito
a<-matrix(c(0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5),2, 4,byrow=T)#Biting rate
Lambda_v<-matrix(c(0.3,0.3,0.3,0.3,0.3,0.3,0.3,0.3)/365, 2, 4,byrow=T)#Recruitment/birth rate of mosquitoes
beta_v<-matrix(c(0.42,0.42,0.42,0.42,0.42,0.42,0.42,0.42), 2, 4,byrow=T)#Transmission rate from infected mosquito to human


params<-list(Lambda_h=Lambda_h, Lambda_v=Lambda_v, beta_v=beta_v, beta_h=beta_h, sigma_h=sigma_h, gamma_h=gamma_h,delta_h=delta_h,mu_h=mu_h, mu_v=mu_v, psi=psi, a=a, m=m)
#Variables
#Initial conditions at DFE
Sh<-matrix(c(1250,1250,1250,1250,1250,1250,1250,1250),2,4,byrow=T,dimnames = list(c("WU","Sl"),c("A1", "A2","A3","A4")))#Initial value of susceptible human
Rh<-matrix(c(100, 100, 100, 100,100, 100, 100, 100),2,4,byrow=T)#Initial value of recovered human
Ih<-matrix(c(300,300,300,300,600,600,600,600),2,4,byrow=T)#Initial value of infected human
Sv<-matrix(c(1000,1000,1000,1000,1500,1500,1500,1500),2,4,byrow=T)
Iv<-matrix(c(100,100,100,100,160,160,160,160),2,4,byrow=T)

# Constructing time vector
t_start<-0 #Start
t_end<-(1825)-1 # End 
t_inc<-.2
t_range<-seq(t_start,t_end+t_inc, t_inc)
# Create a function Sir_si for the differential equations
sir_si <- function(time, state, params) {
  nh<-matrix(2,4)
  Sh<-matrix(2,4, data=state[2:8]) 
  Ih<-matrix(2,4, data=state[9:16])
  Rh<-matrix(2,4, data=state[17:24])
  Sv<-matrix(2,4, data=state[25:32])  
  Iv<-matrix(2,4, data=state[33:40]) 
  nh<- Sh+ Ih + Rh
  with(as.list(c(state, params)), { 
    #for(j in 1:2){#Loop for the 2 patches
      #for (i in 1: 4){#Loop for age groups
          dSh <-(Lambda_h) - (1-psi)*beta_v*a*Iv*Sh/(nh) - (sigma_h*Rh)-(mu_h*Sh)+(m*Sh)
          dIh <- (((1-psi)*beta_v*a*Sh*Iv)/(nh))  - ((gamma_h + mu_h + delta_h)*Ih)
          dRh <- (gamma_h*Ih) - ((mu_h+sigma_h)*Rh) + (m*Rh)
          dSv <- (Lambda_v) - (mu_v*Sv) - (((1-psi)*beta_h*a*Ih*Sv)/(nh)) - (psi*a*Sv)
          dIv <- (((1-psi)*beta_h*a*Ih*Sv)/(nh)) - (psi*a*Iv) - (mu_v*Iv)
        #}
      #}
    list(c(dSh, dIh, dRh, dSv, dIv))
    
  })
}

init<- c(Sh=Sh, Ih=Ih, Rh=Rh, Sv=Sv, Iv=Iv)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = t_range, func = sir_si, parms = params)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
#out$time <- NULL
## Show data
head(out, 5)
View(out)

## Plot
 # A general concern about the plot they are all straight lines plus infected number are negative look like something is wrong
matplot(x = out[,1], y = out[,2:9], type = "l",
        xlab = "Time", ylab = "Number of people", main = "Susceptible by age group",
        lwd = 1, lty = 1, bty = "l", col = 1:8)
legend("topright", c("Sh11", "Sh12", "Sh13", "Sh14", "Sh21", "Sh22", "Sh23", "Sh24"), pch = 1, col = 1:8, bty = "n")

matplot(x = out[,1], y = out[,10:17], type = "l",
        xlab = "Time", ylab = "Number of people", main = "Infected by age group",
        lwd = 1, lty = 1, bty = "l", col = 1:8)
legend("topright", c("Ih11", "Ih12", "Ih13", "Ih14", "Ih21", "Ih22", "Ih23", "Ih24"), pch = 1, col = 1:8, bty = "n")

matplot(x = out[,1], y = out[,18:25], type = "l",
        xlab = "Time", ylab = "Number of people", main = "Recovered by age group",
        lwd = 1, lty = 1, bty = "l", col = 1:8)
legend("topright", c("Rh11", "Rh12", "Rh13", "Rh14", "Rh21", "Rh22", "Rh23", "Rh24"), pch = 1, col = 1:8, bty = "n")

#matplot(x = t_range, y = out[,25:32], type = "l",
 #       xlab = "Time", ylab = "Number of mosquitoes", main = "Susceptible mosquitoes ",
#      lwd = 1, lty = 1, bty = "l", col = 1:8)

#matplot(x = t_range, y = out[,33:40], type = "l",
#        xlab = "Time", ylab = "Number of mosquitoes", main = "Infected mosquitoes ",
#        lwd = 1, lty = 1, bty = "l", col = 1:8)


 #Compute the reproductive rate per age group within patches
#NGM method
r0<-function(parms){
  wf<-Lambda_v*beta_v*a*(1-psi)
  th<-(mu_v+a*psi)^2
  pg<-beta_h*a*(1-psi)*(mu_h-m)
  kt<-Lambda_h*(gamma_h+mu_h+sigma_h)
  r0<-sqrt((wf*pg)/(th*kt))
     return(r0)
 }

r0(params)
# Second scenario decade rate 
x1<-c(0.6,0.6,0.6,0.6,0.6,0.6,0.6,0.6)
x2<-c(1.62,1.62,1.62,1.62,1,1,1,1)*365
b<-function(x1,x2){x1*exp(-x2*t_range)} # x1 represent the value of psi and x2 is either 1.62 or 1
psi<-matrix(c(b(x1[1],x2[1]),b(x1[2],x2[2]), b(x1[3],x2[3]),b(x1[4],x2[4]), b(x1[5],x2[5]), b(x1[6],x2[6]), b(x1[7],x2[7]), b(x1[8],x2[8]) ), 2, 4,byrow=T)#Proportion of ITN use 
params<-list(Lambda_h=Lambda_h, Lambda_v=Lambda_v, beta_v=beta_v, beta_h=beta_h, sigma_h=sigma_h, gamma_h=gamma_h,delta_h=delta_h,mu_h=mu_h, mu_v=mu_v, psi=psi, a=a, m=m)
init<- c(Sh=Sh, Ih=Ih, Rh=Rh, Sv=Sv, Iv=Iv)

## Solve using ode (General Solver for Ordinary Differential Equations)
out <- ode(y = init, times = t_range, func = sir_si, parms = params)
## change to data frame
out <- as.data.frame(out)
## Delete time variable
#out$time <- NULL
## Show data
head(out, 5)
View(out)
# Third scenario 
d<-function(x1,x2){x1*x2}
x3<-c(0.8,0.8,0.8,0.8,0.5,0.5,0.5,0.5)
psi<-matrix(c(d(x1[1],x3[1]),d(x1[2],x3[2]), d(x1[3],x3[3]),d(x1[4],x3[4]), d(x1[5],x3[5]), d(x1[6],x3[6]), d(x1[7],x3[7]), d(x1[8],x3[8]) ), 2, 4,byrow=T)#Proportion of ITN use 



matplot(x = out[,1], y = out[,2:9], type = "l",
        xlab = "Time", ylab = "Number of people", main = "Susceptible by age group",
        lwd = 1, lty = 1, bty = "l", col = 1:8)
legend("topright", c("Sh11", "Sh12", "Sh13", "Sh14", "Sh21", "Sh22", "Sh23", "Sh24"), pch = 1, col = 1:8, bty = "n")

matplot(x = out[,1], y = out[,10:17], type = "l",
        xlab = "Time", ylab = "Number of people", main = "Infected by age group",
        lwd = 1, lty = 1, bty = "l", col = 1:8)
legend("topright", c("Ih11", "Ih12", "Ih13", "Ih14", "Ih21", "Ih22", "Ih23", "Ih24"), pch = 1, col = 1:8, bty = "n")

matplot(x = out[,1], y = out[,18:25], type = "l",
        xlab = "Time", ylab = "Number of people", main = "Recovered by age group",
        lwd = 1, lty = 1, bty = "l", col = 1:8)
legend("topright", c("Rh11", "Rh12", "Rh13", "Rh14", "Rh21", "Rh22", "Rh23", "Rh24"), pch = 1, col = 1:8, bty = "n")
# Sensitivity analysis
#require(ODEsensitivity)
#require(fast)#Fourier Amplitude Sensitivity Analysis
# Latin hypercube sampling 
require(lhs)

