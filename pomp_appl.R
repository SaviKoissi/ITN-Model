## Iterative filtering method for the detection of the parameter 

# Data importation 
require(tidyverse)
setwd("C:/Users/ZEF/Desktop/MS3/Parameters/ITN-Model")
data<-read.csv("C:/Users/ZEF/Downloads/data_spatial.csv", h=T, na.strings="NA")
#head(data)
dt<-data %>%
  tibble()%>%
  filter (Region == "Greater Accra")%>%
  unite(MonthYear, Month, Year, sep = " ")%>%
  mutate(MonthYear= parse_date(MonthYear, format = "%B %Y")) %>%
  mutate(#Age = stringr::str_replace(Age, "days.", ""),
         Age = stringr::str_replace(Age, "....28days.", "1"),
         Age = stringr::str_replace(Age, "...1.11mths.", "1"),
         Age = stringr::str_replace(Age, "...1.4.", "1"),
         Age = stringr::str_replace(Age, "...5.9.", "1"),
         Age = stringr::str_replace(Age, "...10.14.", "1"),
         Age = stringr::str_replace(Age, "...15.17.", "2"),
         Age = stringr::str_replace(Age, "...18.19.", "2"),
         Age = stringr::str_replace(Age, "...20.34.", "2"),
         Age = stringr::str_replace(Age, "...35.49.", "2"),
         Age = stringr::str_replace(Age, "...50.59.", "2"),
         Age = stringr::str_replace(Age, "...60.69.", "2"),
         Age = stringr::str_replace(Age, "...70..", "2"),
         Age = stringr::str_replace(Age, "..28days.", "1"),
         Age = stringr::str_replace(Age, ".1.11mths.", "1"),
         Age = stringr::str_replace(Age, ".1.4.", "1"),
         Age = stringr::str_replace(Age, ".5.9.", "1"),
         Age = stringr::str_replace(Age, ".10.14.", "1"),
         Age = stringr::str_replace(Age, ".15.17.", "2"),
         Age = stringr::str_replace(Age, ".18.19.", "2"),
         Age = stringr::str_replace(Age, ".20.34.", "2"),
         Age = stringr::str_replace(Age, ".35.49.", "2"),
         Age = stringr::str_replace(Age, ".50.59.", "2"),
         Age = stringr::str_replace(Age, ".60.69.", "2"),
         Age = stringr::str_replace(Age, ".70...", "2"),
         Age = stringr::str_replace(Age, ".70..", "2")
                  ) %>%
  group_by(Region, MonthYear, Age) %>%
  mutate(Infected = sum(value, na.rm = TRUE)) %>%
  select(-c(1, 3, 5:7)) %>% 
  
  dt %>% 
    ggplot(aes(x = MonthYear, y= Infected, group = Age))+
    geom_line()
 
 # Demographic data 

demo_data<-read.csv("C:/Users/ZEF/Desktop/MS3/WorlPop/Popu_15-18_gha.csv")

demo<-demo_data %>%
  tibble() %>%
  select(-c(1:7,9:11)) %>%
  filter(RGN_NM2012 == "Greater Accra") %>%
  group_by(Year, Age) %>%
  summarise(summa = sum(Estimate, na.rm = TRUE)) %>%
  mutate(birth = ifelse(Age == 0 , summa, NA)) %>%
  mutate(pop = sum(summa)) %>%
  select(-c(2,3)) %>%
  drop_na() %>%
  mutate(frac_birth = birth /pop)

#Adjust the birthrate using the spline function. In this case, we are running to a problem cause the size of data== 4 < 30
#covar<-demo %>%
 # summarise(
  #  time = seq(from = min(Year), to = max(Year), by = 1/12),
   # pop = predict(smooth.spline(x = Year, y = pop), x = time)$y, 
    #birthrate = predict(smooth.spline(x = Year + 0.5, y = birth), x = time - 4)$y
  #)

  
  
# Partially Observed Markov Process for the determination of the parameters of the model

dt %>%
  filter(Age == 1) ->dat_inf
  #select(MonthYear, Infected)%>%
  
dat_inf %>%
  ggplot(aes(x=MonthYear, y=Infected)) +
  geom_line(col = "blue")+
  theme_bw()

#double dN_SIv = rbinom (Sv, 1 - exp(( beta_h * a * Ih * Sv) / nh) * dt));
#Sv -= dN_SIv;
#Iv += dN_SIv;
#Sv = nearbyint (eta2 * nh);
require(pomp)
sir_si_step <- Csnippet( "
double dN_SI = rbinom (Sh, 1- exp (- beta_v * a  * Ih * Sh / nh * dt)); 
double dN_IR = rbinom (Ih, 1- exp(-(gamma_h + mu_h + delta_h) * dt )); 
Sh -= dN_SI; 
Ih += dN_SI-dN_IR; 
Hh += dN_SI;
                         ")

sir_si_init<- Csnippet("
Sh = nearbyint (eta * nh);
Ih = 1; 
Hh = 0; 
                     ")

dmal<- Csnippet("
double ll = dbinom(Infected, Hh, rho, give_log);
lik =  (!isfinite(ll) ? -1000 : ll );
                ")

rmal<- Csnippet("
Infected = rbinom (Hh, rho);
                ")

dat_inf %>%
  unique()%>%
  mutate(Weeks= as.numeric(difftime(MonthYear,as.Date("2015-01-01","%Y-%m-%d") ,units = "weeks")))%>%
  pomp(
    times = "Weeks", t0 = 0, 
    rprocess = euler (sir_si_step, delta.t = 1/7), 
    rinit = sir_si_init, 
    rmeasure = rmal, 
    dmeasure = dmal, 
    accumvars = "Hh", 
    partrans = parameter_trans(
      log = c("beta_v", "beta_h", "a", "gamma_h", "mu_h", "delta_h"),
      logit = c("rho", "eta")
    ), 
    statenames = c("Sh", "Ih", "Hh"), 
    paramnames = c("beta_v", "beta_h", "a", "gamma_h", "rho", "eta", "nh","mu_h", "delta_h")
    
  )-> malSIRSI

params <- c(beta_v = 0.024, beta_h = 0.024, a = 0.35, gamma_h = 0.00274, rho = 0.5, eta = 0.1,
            nh = 5298358, mu_h = 1/(63*365), delta_h = 0.0003454)

malSIRSI %>%
  simulate(params= params, nsim = 10, format = "data.frame") -> y

y %>%
  ggplot(aes(x= Weeks, y = Infected, group = .id, color = factor(.id)))+
  geom_line()+
  theme_bw()+
  scale_color_brewer(type="qual",palette=3)+
  guides(color=FALSE)

malSIRSI %>%
  pfilter(Np=1000,params=params) -> pf

plot(pf)

fixed_params <- c(N=5298358, mu_h=1/(63*365))

# Parameter estimation

require(foreach)
require(doParallel)
registerDoParallel()
require(doRNG)
registerDoRNG(625904618)

foreach(i=1:10,.combine=c) %dopar% {
  require(pomp)
  require(dplyr)
  malSIRSI %>% pfilter(params=params,Np=10000)
   }-> pf


pf %>% 
  logLik() %>% 
  logmeanexp(se=TRUE) -> L_pf
L_pf


pf[[1]] %>% 
  coef() %>% 
  bind_rows() %>%
  bind_cols(loglik=L_pf[1],loglik.se=L_pf[2])%>%
  saveRDS("mal_params.rds")


foreach(i=1:20,.combine=c) %dopar% {
  require(pomp)
  require(dplyr)
  malSIRSI %>%
    mif2(
      params=params,
      Np=2000, Nmif=50,
      cooling.fraction.50=0.5,
      rw.sd=rw.sd(beta_v=0.1, beta_h=1, gamma_h = 0.1, a= 0.01, mu_h= 0.1, delta_h = 0.1, rho = 0.02,
                  eta=ivp(0.02))
    )
  } -> mifs_local


mifs_local %>%
  traces() %>%
  melt() %>%
  ggplot(aes(x=iteration,y=value,group=L1,color=factor(L1)))+
  geom_line()+
  guides(color=FALSE)+
  facet_wrap(~variable,scales="free_y")

foreach(mf=mifs_local,.combine=rbind) %dopar% {
  require(pomp)
  require(dplyr)
  evals <- replicate(10, logLik(pfilter(mf,Np=20000)))
  ll <- logmeanexp(evals,se=TRUE)
  mf %>% coef() %>% bind_rows() %>%
  bind_cols(loglik=ll[1],loglik.se=ll[2])
  } -> results

pairs(~loglik+beta_v+eta+rho,data=results,pch=16)

readRDS(file = "mal_params.rds")%>%
  bind_rows(results)%>%
  arrange(-loglik)%>%
  saveRDS("mal_params.rds")

# Searching the MLE









with(as.list(c(state, params)), {
    dSh <- (Lambda_h) - (1 - psi) * beta_v * a * Iv * Sh / nh - (sigma_h * Rh) - (mu_h * Sh) + (m * Sh)
    dIh <- (((1 - psi) * beta_v * a * Sh * Iv) / nh) - ((gamma_h + mu_h + delta_h) * Ih)
    dRh <- (gamma_h * Ih) - ((mu_h + sigma_h) * Rh) + (m * Rh)
    dSv <- (Lambda_v) - (mu_v * Sv) - (((1 - psi) * beta_h * a * Ih * Sv) / (nh)) - (psi * a * Sv)
    dIv <- (((1 - psi) * beta_h * a * Ih * Sv) / (nh)) - (psi * a * Iv) - (mu_v * Iv)
    list(c(dSh, dIh, dRh, dSv, dIv))
  })
                        