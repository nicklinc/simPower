notice = "
# The simPower_b() and simPower_w() functions (Nicklin & Vitta, in preparation) 
# are powered by Norouzian et al.'s (2018; 2019) \"bayesL2\" a suite of R functions
# for Bayesian estimation in L2 research.
# Copyright (C) 2017-present  Reza Norouzian, rnorouzian@gmail.com
# For the full set of packages and more information:
# source(\"https://raw.githubusercontent.com/rnorouzian/i/master/i.r\")
#
# NOTE: Requires the RVAideMemoire package"


message(notice)

library(RVAideMemoire)
# power.t() (Norouzian et al., 2018; 2019)
power.t <- function(d = .1, sig.level = .05, power = .8, base.rate = 1, paired = FALSE, two.tailed = TRUE)
{
  
  pwr <- Vectorize(function(d, sig.level, power, base.rate, paired, two.tailed)
  {
    
    d <- abs(d)
    sig.level <- if(two.tailed) sig.level/2 else sig.level
    k <- base.rate / (1 + base.rate)
    options(warn = -1)
    
    f <- if(two.tailed){ function(x){
      
      power - (pt(qt(sig.level, df = x), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2))) + pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE))
    }
      
    } else {
      
      function(x){
        
        power - pt(qt(sig.level, df = x, lower.tail = FALSE), df = x, ncp = d*sqrt(if(paired) x + 1 else k*(x + 2)), lower.tail = FALSE)
        
      }
    }
    
    df <- ceiling(uniroot(f, c(1e-8, 1e7), extendInt = "downX")[[1]])
    
    n1 <- df + 1
    n2 <- if(paired) NA else round(base.rate*n1)
    
    return(c(n1 = n1, n2 = n2))
  })
  
  data.frame(t(pwr(d = d, sig.level = sig.level, power = power, base.rate = base.rate, paired = paired, two.tailed = two.tailed)))
}

# BF.t() (Norouzian et al., 2018; 2019)
BF.t <- function(t, n1, n2 = NA, scale = sqrt(2)/2, log.BF = FALSE){
  f <- Vectorize(function(t, n1, n2, scale, log.BF){
    options(warn = -1)  
    t = abs(t)
    N = ifelse(is.na(n2), n1, (n1 * n2)/(n1 + n2))
    df = ifelse(is.na(n2), n1 - 1, (n1 + n2) - 2)
    d = t / sqrt(N)
    H1 = integrate(function(delta) dcauchy(delta, 0, scale) * dt(t, df, delta*sqrt(N)), -Inf, Inf)[[1]]
    H0 = dt(t, df)
    BF10 = ifelse(log.BF, log(H1/H0), H1/H0)
    p.value = 2*(1-pt(t, df))
    return(c(t = t, BF10 = BF10, H0 = H0, H1 = H1, p.value = p.value, d = d,  n1 = n1, n2 = n2))
  })
  data.frame(t(f(t = t, n1 = n1, n2 = n2, scale = scale, log.BF = log.BF)))
} 

# simPower_b (Nicklin & Vitta, in preparation)
simPower_b <- function(d, n.eq = TRUE, csv = TRUE){
  n1 = vector()
  n2 = vector()
  n_tot = vector()
  t_stat = vector()
  BF10_med = vector()
  BF10_wide = vector()
  BF10_v.wide = vector()
  p=seq(0.99,0.01,-.01)
  i=1
if(n.eq == TRUE){for(i in i:99){
  n1 <- append(n1, as.numeric(power.t(d=d, power = p[i])[1]))
  n2 <- append(n2, as.numeric(power.t(d=d, power = p[i])[2]))
  n_tot <- append(n_tot, as.numeric(power.t(d=d, power = p[i])[1]) + as.numeric(power.t(d=d, power = p[i])[2]))
  t_stat <- append(t_stat, as.numeric(d*sqrt((as.numeric(power.t(d=d, power = p[i])[1])-1)+(as.numeric(power.t(d=d, power = p[i])[2])-1))/2))
  BF10_med <- append(BF10_med, BF.t(t_stat[i],n1[i],n2[i],.5)[1,2])
  BF10_wide <- append(BF10_wide, BF.t(t_stat[i],n1[i],n2[i],.707)[1,2])
  BF10_v.wide <- append(BF10_v.wide, BF.t(t_stat[i],n1[i],n2[i],1)[1,2])
  print(i)
  i=i+1
}
  sim_data <- data.frame(p, n1, n2, n_tot, t_stat, BF10_med, BF10_wide, BF10_v.wide)
  print(sim_data)
  if (csv == FALSE){print("NOTE: No .csv called")}  else if (csv == TRUE){write.csv(sim_data, "simPower_data(btw_n_equal).csv")} 
  plot(p, BF10_v.wide)
  spearman.ci(p, BF10_v.wide, 1000, .95)
} else if(n.eq == FALSE){for(i in i:99){
  n1 <- append(n1, floor(as.numeric(power.t(d=d, power = p[i])[2]) - (as.numeric(power.t(d=d, power = p[i])[2])*.05)))
  n2 <- append(n2, ceiling(as.numeric(power.t(d=d, power = p[i])[2]) + (as.numeric(power.t(d=d, power = p[i])[2])*.05)))
  n_tot <- append(n_tot, as.numeric(power.t(d=d, power = p[i])[1]) + as.numeric(power.t(d=d, power = p[i])[2]))
  t_stat <- append(t_stat, as.numeric(d*sqrt((as.numeric(power.t(d=d, power = p[i])[1])-1)+(as.numeric(power.t(d=d, power = p[i])[2])-1))/2))
  BF10_med <- append(BF10_med, BF.t(t_stat[i],n1[i],n2[i],.5)[1,2])
  BF10_wide <- append(BF10_wide, BF.t(t_stat[i],n1[i],n2[i],.707)[1,2])
  BF10_v.wide <- append(BF10_v.wide, BF.t(t_stat[i],n1[i],n2[i],1)[1,2])
  print(i)
  i=i+1
}
  sim_data <- data.frame(p, n1, n2, n_tot, t_stat, BF10_med, BF10_wide, BF10_v.wide)
  print(sim_data)
  if (csv == FALSE){print("NOTE: No .csv called")}  else if (csv == TRUE){write.csv(sim_data, "simPower_data(btw_n_dif).csv")} 
  plot(p, BF10_v.wide)
  spearman.ci(p, BF10_v.wide, 1000, .95)
  } 
}

# simPower_w() (Nicklin & Vitta, in preparation)
simPower_w <- function(d, csv = TRUE){
  n1 = vector()
  n2 = vector()
  n_tot = vector()
  t_stat = vector()
  BF10_med = vector()
  BF10_wide = vector()
  BF10_v.wide = vector()
  p=seq(0.99,0.01,-.01)
  i=1
  for(i in i:99){
    n1 <- append(n1, as.numeric(power.t(d=d, power = p[i], paired = TRUE)[1]))
    n2 <- append(n2, NA)
    n_tot <- append(n_tot, as.numeric(power.t(d=d, power = p[i], paired = TRUE)[1]))
    t_stat <- append(t_stat, as.numeric(d*sqrt(as.numeric(power.t(d=d, power = p[i], paired = TRUE)[1]))))
    w_t <- as.numeric(d*sqrt(as.numeric(power.t(d=d, power = p[i], paired = TRUE)[1])))
    BF10_med <- append(BF10_med, BF.t(w_t, as.numeric(power.t(d=d, power = p[i], paired = TRUE)[1]), n2=NA, .5)[1,2])
    BF10_wide <- append(BF10_wide, BF.t(w_t, as.numeric(power.t(d=d, power = p[i], paired = TRUE)[1]), n2=NA, .707)[1,2])
    BF10_v.wide <- append(BF10_v.wide, BF.t(w_t, as.numeric(power.t(d=d, power = p[i], paired = TRUE)[1]), n2=NA, 1)[1,2])
    print(i)
    i=i+1
  }
  sim <- data.frame(p, n1, n2, n_tot, t_stat, BF10_med, BF10_wide, BF10_v.wide)
  print(sim)
  if (csv == FALSE){print("NOTE: No .csv called")}  else if (csv == TRUE){write.csv(sim, "simPower_data(within).csv")} 
  plot(p, BF10_v.wide)
  spearman.ci(p, BF10_v.wide, 1000, .95)
}

