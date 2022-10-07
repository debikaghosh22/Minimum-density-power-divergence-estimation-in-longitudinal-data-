install.packages("JM")
library(JM)
library(dplyr)
library(tidyverse)
library(data.table)
library(Rlab)
library(HDInterval)
library(mvtnorm)
library(matrixStats)
data("pbc")
data <- pbcseq
data <- data %>%  select(id,hepato,age,trt,ascites,spiders,bili,albumin,platelet,protime)

u <- data[rowSums(is.na(data)) > 0, ]
dat <- drop_na(data)
## transforming some variables:
dat$bili <- log(dat$bili)
dat$albumin <- log(dat$albumin)
dat$platelet <- log(dat$platelet)
dat$protime <- log(dat$protime)

colnames(dat)[1] <- "index"
num <- c()
for(i in 1:312){
  if((colCounts(as.matrix(dat),value = i)[1])==1) {
    num <- append(num,i)
  }
}
new_dat <- subset(dat, !(index %in% num))


## Get the DPD  estimates: 



dpd_pairwise <- function(x,fixed = c(rep(FALSE,11))){
  params <- fixed
  dpd_1 <- function(p){
    params[!fixed] <- p
    alpha <- params[1]
    beta_0 <- params[2]
    beta_1 <- params[3]
    beta_2 <- params[4]
    beta_3 <- params[5]
    beta_4 <- params[6]
    beta_5 <- params[7]
    beta_6 <- params[8]
    beta_7 <- params[9]
    beta_8 <- params[10]
    rho <- params[11]
    add_pi <- function(d){
      k <- beta_0+(d[3]*beta_1)+(d[4]*beta_2)+(d[5]*beta_3)+(d[6]*beta_4)+(d[7]*beta_5)+(d[8]*beta_6)+(d[9]*beta_7)+(d[10]*beta_8)
      k1 <- exp(k)/(1+exp(k))
      d <- cbind(d,k1)
    }
    dat_split <- split(x , f  = x$index)
    result <- lapply(dat_split, add_pi)
    
    result <- rbindlist(result)
    result <- as.data.frame(result)
    colnames(result) <- c('index','y','x1','x2','x3','x4','x5','x6','x7','x8','exp_prob')
    result_split <- split(result, f = result$index)
    
    expression <- function(d){
      bin <- as.data.frame(combn(d$y , 2))
      pr <- as.data.frame(combn(d$exp_prob , 2))
      
      ## Evaluation of the probabilities:
      f_jk <- function(u,v){
        
        
        dummy_func <- function(x,y){
          pbern(x, prob = y)
        }
        dummy_func_1 <- function(x,y){
          pbern(x-1, prob = y)
        }
        
        k <- mapply(dummy_func,u,v)
        k_1 <- mapply(dummy_func_1,u,v)
        inv1 <- inverseCDF(as.matrix(k), pnorm)
        inv2 <- inverseCDF(as.matrix(k_1), pnorm)
        mean <- rep(0,2)
        lower <- inv2
        upper <- inv1
        corr <- diag(2)
        corr[lower.tri(corr)] <- rho
        corr[upper.tri(corr)] <- rho
        prob <- pmvnorm(lower = lower, upper = upper, mean = mean, corr = corr)
        prob <- (1+(1/alpha))*(prob^alpha)
        ## First expression:
        all_possib <- as.data.frame(expand.grid(c(0,1),c(0,1)))
        all_possib <- as.data.frame(t(all_possib))
        val <- c()
        for(i in 1:ncol(all_possib)){
          u <- all_possib[i]
          u <- u[[1]]
          k <- mapply(dummy_func,u,v)
          k_1 <- mapply(dummy_func_1,u,v)
          inv1_1 <- inverseCDF(as.matrix(k), pnorm)
          inv2_1 <- inverseCDF(as.matrix(k_1), pnorm)
          mean1 <- rep(0,2)
          lower1 <- inv2_1
          upper1 <- inv1_1
          corr1 <- diag(2)
          corr1[lower.tri(corr1)] <- rho
          corr1[upper.tri(corr1)] <- rho
          prob1 <- pmvnorm(lower = lower1, upper = upper1, mean = mean1, corr = corr1)
          prob1 <- prob1^(1+alpha)
          
          val <- append(val,prob1)
        }
        
        val_s <- sum(val)
        
        return(val_s - prob)
      }
      
      final_res <- mapply(f_jk, bin, pr)
      
      final_value <- sum(final_res)
    }
    
    u <- sapply(result_split,expression)
    return(sum(u))
  }
}

val <- dpd_pairwise(new_dat,c(0.1,rep(FALSE,10)))
### The following IVs have lead to convergence = 0:
optim(c(beta_0 = 5,beta_1=0.01,beta_2= -0.02,beta_3= 0.4,beta_4= 0.5,beta_5 = 0.5,beta_6 = -2,beta_7 = -0.3,beta_8 = -0.8,rho = 0.4),val, control = list(maxit=3000))









## for variance:
one <- rep(1,1834)
new_dat <- add_column(new_dat,one, .after =2) ## only for variance
colnames(new_dat) <- c("index","y","x0","x1","x2","x3",'x4','x5','x6','x7','x8')

#VARIANCE:

## Omega matrix:


omega <- function(dat,alpha,beta_0,beta_1,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7,beta_8,rho){
  dat_split <- split(dat, f= dat$index)
  
  xi <- function(y1,y2,b1,b2,xit1,xit2){
    z1 <- pbern(y1,prob = b1/(1+b1))
    z1 <- inverseCDF(z1,pnorm) # final
    z2 <- pbern(y2,prob = b2/(1+b2))
    z2 <- inverseCDF(z2,pnorm) # final
    z1_bar <- pbern(y1-1,prob = b1/(1+b1))
    z1_bar <- inverseCDF(z1_bar,pnorm) # final
    z2_bar <- pbern(y2-1,prob = b2/(1+b2))
    z2_bar <- inverseCDF(z2_bar,pnorm) # final
    
    t1 <- (z1 - rho*z2)/sqrt(1-rho^2)
    t2 <- (z2 - rho*z1)/sqrt(1-rho^2)
    t3 <- (z1_bar - rho*z2)/sqrt(1-rho^2)
    t4 <- (z2 - rho*z1_bar)/sqrt(1-rho^2)
    t5 <- (z2_bar - rho*z1)/sqrt(1-rho^2)
    t6 <- (z1 - rho*z2_bar)/sqrt(1-rho^2)
    t7 <- (z1_bar - rho*z2_bar)/sqrt(1-rho^2)
    t8 <- (z2_bar - rho*z1_bar)/sqrt(1-rho^2)
    
    k1 <- ( pnorm(t2,0,1)*(-b1/(1+b1)^2)*(1-y1) - pnorm(t4,0,1)*(-b1/(1+b1)^2)*y1  - pnorm(t5,0,1)*(-b1/(1+b1)^2)*(1-y1) + pnorm(t8,0,1)*(-b1/(1+b1)^2)*y1 ) 
    k2 <- ( pnorm(t1,0,1)*(-b2/(1+b2)^2)*(1-y2) - pnorm(t3,0,1)*(-b2/(1+b2)^2)*(1-y2) - pnorm(t6,0,1)*(-b2/(1+b2)^2)*y2  + pnorm(t7,0,1)*(-b2/(1+b2)^2)*y2 )
    
    k <- (k1 * xit1 + k2* xit2)
    
    ## for rho:
    m1 <- (z1^2 + z2^2 - 2*z1*z2*rho )/(2*(1-rho^2))
    m2 <- exp(-m1)
    m3 <- (z1_bar^2 + z2^2 - 2*z1_bar*z2*rho )/(2*(1-rho^2))
    m4 <- exp(-m3)
    m5 <- (z1^2 + z2_bar^2 - 2*z1*z2_bar*rho )/(2*(1-rho^2))
    m6 <- exp(-m5)
    m7 <- (z1_bar^2 + z2_bar^2 - 2*z1_bar*z2_bar*rho )/(2*(1-rho^2))
    m8 <- exp(-m7)
    
    
    
    r <- 1/(2*pi*sqrt(1-rho^2))
    p <- ( m2-m4-m6+m8)
    exp_rho <- r*p
    res1 <- c(k,exp_rho)
    
    f <- pbinorm(z1,z2,0,0,1,1,rho)-pbinorm(z1_bar,z2,0,0,1,1,rho)-pbinorm(z1,z2_bar,0,0,1,1,rho)+pbinorm(z1_bar,z2_bar,0,0,1,1,rho)
    final <- (1/f)*f^(1+alpha)*res1
    return(final)
  }
  
  k_mat <- function(y1,y2,b1,b2,xit1,xit2){
    
    z1 <- pbern(y1,prob = b1/(1+b1))
    z1 <- inverseCDF(z1,pnorm) # final
    z2 <- pbern(y2,prob = b2/(1+b2))
    z2 <- inverseCDF(z2,pnorm) # final
    z1_bar <- pbern(y1-1,prob = b1/(1+b1))
    z1_bar <- inverseCDF(z1_bar,pnorm) # final
    z2_bar <- pbern(y2-1,prob = b2/(1+b2))
    z2_bar <- inverseCDF(z2_bar,pnorm) # final
    
    t1 <- (z1 - rho*z2)/sqrt(1-rho^2)
    t2 <- (z2 - rho*z1)/sqrt(1-rho^2)
    t3 <- (z1_bar - rho*z2)/sqrt(1-rho^2)
    t4 <- (z2 - rho*z1_bar)/sqrt(1-rho^2)
    t5 <- (z2_bar - rho*z1)/sqrt(1-rho^2)
    t6 <- (z1 - rho*z2_bar)/sqrt(1-rho^2)
    t7 <- (z1_bar - rho*z2_bar)/sqrt(1-rho^2)
    t8 <- (z2_bar - rho*z1_bar)/sqrt(1-rho^2)
    
    k1 <- ( pnorm(t2,0,1)*(-b1/(1+b1)^2)*(1-y1) - pnorm(t4,0,1)*(-b1/(1+b1)^2)*y1  - pnorm(t5,0,1)*(-b1/(1+b1)^2)*(1-y1) + pnorm(t8,0,1)*(-b1/(1+b1)^2)*y1 ) 
    k2 <- ( pnorm(t1,0,1)*(-b2/(1+b2)^2)*(1-y2) - pnorm(t3,0,1)*(-b2/(1+b2)^2)*(1-y2) - pnorm(t6,0,1)*(-b2/(1+b2)^2)*y2  + pnorm(t7,0,1)*(-b2/(1+b2)^2)*y2 )
    
    k <- (k1 * xit1 + k2* xit2)
    
    ## for rho:
    m1 <- (z1^2 + z2^2 - 2*z1*z2*rho )/(2*(1-rho^2))
    m2 <- exp(-m1)
    m3 <- (z1_bar^2 + z2^2 - 2*z1_bar*z2*rho )/(2*(1-rho^2))
    m4 <- exp(-m3)
    m5 <- (z1^2 + z2_bar^2 - 2*z1*z2_bar*rho )/(2*(1-rho^2))
    m6 <- exp(-m5)
    m7 <- (z1_bar^2 + z2_bar^2 - 2*z1_bar*z2_bar*rho )/(2*(1-rho^2))
    m8 <- exp(-m7)
    
    
    
    r <- 1/(2*pi*sqrt(1-rho^2))
    p <- ( m2-m4-m6+m8)
    exp_rho <- r*p
    res1 <- c(k,exp_rho)
    
    f <- pbinorm(z1,z2,0,0,1,1,rho)-pbinorm(z1_bar,z2,0,0,1,1,rho)-pbinorm(z1,z2_bar,0,0,1,1,rho)+pbinorm(z1_bar,z2_bar,0,0,1,1,rho)
    
    
    u1 <- (1/f)*res1
    u2 <- t((1/f)*res1)
    u3 <- u1 %*% u2
    u4 <- u3 * f^(2*alpha+1)
    return(u4)
  }
  omega_1 <- function(y1,y2,b1,b2,xit1,xit2){
    t1 <- xi(0,0,b1,b2,xit1,xit2)+xi(0,1,b1,b2,xit1,xit2)+xi(1,0,b1,b2,xit1,xit2)+xi(1,1,b1,b2,xit1,xit2)
    t2 <- t1 %*% t(t1)
    t3 <- k_mat(0,0,b1,b2,xit1,xit2)+k_mat(0,1,b1,b2,xit1,xit2)+k_mat(1,0,b1,b2,xit1,xit2)+k_mat(1,1,b1,b2,xit1,xit2)
    return(t3-t2)
  }
  final_var <- function(d){
    a <- data.frame(combn(d$y,2))
    y1 <- as.list(a[1,])
    y2 <- as.list(a[2,])   ## responses
    
    w <- exp(beta_0*d[3]+beta_1*d[4]+beta_2*d[5]+beta_3*d[6]+beta_4*d[7]+beta_5*d[8]+beta_6*d[9]+beta_7*d[10]+beta_8*d[11])
    beta <- data.frame(combn(w$x0,2))
    c1 <- as.list(beta[1,])
    c2 <- as.list(beta[2,])  ## values of b
    
    d2 <- data.frame(combn(d$x0,2))
    d3 <- data.frame(combn(d$x1,2))
    d4 <- data.frame(combn(d$x2,2))
    d5 <- data.frame(combn(d$x3,2))
    d6 <- data.frame(combn(d$x4,2))
    d7 <- data.frame(combn(d$x5,2))
    d8 <- data.frame(combn(d$x6,2))
    d9 <- data.frame(combn(d$x7,2))
    d10 <- data.frame(combn(d$x8,2))
    
    make_list_of_vec <- function(a,b,d,e,f,g,h,j,k){   # Will change with covariates
      m <- list()
      for (i in 1:ncol(a)){
        vec <- c(a[1,i],b[1,i],d[1,i],e[1,i],f[1,i],g[1,i],h[1,i],j[1,i],k[1,i])
        m[[i]] <- vec
      }
      m
    }
    make_list_of_vec_1 <- function(a,b,d,e,f,g,h,j,k){   # Will change with covariates
      m <- list()
      for (i in 1:ncol(a)){
        vec <- c(a[2,i],b[2,i],d[2,i],e[2,i],f[2,i],g[2,i],h[2,i],j[2,i],k[2,i])
        m[[i]] <- vec
      }
      m
    }
    
    xit1 <- make_list_of_vec(d2,d3,d4,d5,d6,d7,d8,d9,d10)
    xit2 <- make_list_of_vec_1(d2,d3,d4,d5,d6,d7,d8,d9,d10) ## values of xit
    
    var_3 <- mapply(omega_1,y1,y2,c1,c2,xit1,xit2, SIMPLIFY = FALSE)
    var_4 <- Reduce('+',var_3)
    return(var_4)
  }
  final <- lapply(dat_split,final_var)
}

## PSI matrix:

psi <- function(dat,alpha,beta_0,beta_1,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7,beta_8,rho){
  dat_split <- split(dat, f= dat$index)
  J_mat <- function(y1,y2,b1,b2,xit1,xit2){
    
    z1 <- pbern(y1,prob = b1/(1+b1))
    z1 <- inverseCDF(z1,pnorm) # final
    z2 <- pbern(y2,prob = b2/(1+b2))
    z2 <- inverseCDF(z2,pnorm) # final
    z1_bar <- pbern(y1-1,prob = b1/(1+b1))
    z1_bar <- inverseCDF(z1_bar,pnorm) # final
    z2_bar <- pbern(y2-1,prob = b2/(1+b2))
    z2_bar <- inverseCDF(z2_bar,pnorm) # final
    
    t1 <- (z1 - rho*z2)/sqrt(1-rho^2)
    t2 <- (z2 - rho*z1)/sqrt(1-rho^2)
    t3 <- (z1_bar - rho*z2)/sqrt(1-rho^2)
    t4 <- (z2 - rho*z1_bar)/sqrt(1-rho^2)
    t5 <- (z2_bar - rho*z1)/sqrt(1-rho^2)
    t6 <- (z1 - rho*z2_bar)/sqrt(1-rho^2)
    t7 <- (z1_bar - rho*z2_bar)/sqrt(1-rho^2)
    t8 <- (z2_bar - rho*z1_bar)/sqrt(1-rho^2)
    
    k1 <- ( pnorm(t2,0,1)*(-b1/(1+b1)^2)*(1-y1) - pnorm(t4,0,1)*(-b1/(1+b1)^2)*y1  - pnorm(t5,0,1)*(-b1/(1+b1)^2)*(1-y1) + pnorm(t8,0,1)*(-b1/(1+b1)^2)*y1 ) 
    k2 <- ( pnorm(t1,0,1)*(-b2/(1+b2)^2)*(1-y2) - pnorm(t3,0,1)*(-b2/(1+b2)^2)*(1-y2) - pnorm(t6,0,1)*(-b2/(1+b2)^2)*y2  + pnorm(t7,0,1)*(-b2/(1+b2)^2)*y2 )
    
    k <- (k1 * xit1 + k2* xit2)
    
    ## for rho:
    m1 <- (z1^2 + z2^2 - 2*z1*z2*rho )/(2*(1-rho^2))
    m2 <- exp(-m1)
    m3 <- (z1_bar^2 + z2^2 - 2*z1_bar*z2*rho )/(2*(1-rho^2))
    m4 <- exp(-m3)
    m5 <- (z1^2 + z2_bar^2 - 2*z1*z2_bar*rho )/(2*(1-rho^2))
    m6 <- exp(-m5)
    m7 <- (z1_bar^2 + z2_bar^2 - 2*z1_bar*z2_bar*rho )/(2*(1-rho^2))
    m8 <- exp(-m7)
    
    
    
    r <- 1/(2*pi*sqrt(1-rho^2))
    p <- ( m2-m4-m6+m8)
    exp_rho <- r*p
    res1 <- c(k,exp_rho)
    
    f <- pbinorm(z1,z2,0,0,1,1,rho)-pbinorm(z1_bar,z2,0,0,1,1,rho)-pbinorm(z1,z2_bar,0,0,1,1,rho)+pbinorm(z1_bar,z2_bar,0,0,1,1,rho)
    
    u1 <- (1/f)*res1
    u2 <- t((1/f)*res1)
    u3 <- u1 %*% u2
    u4 <- u3 * f^(alpha+1)
    return(u4) 
  }
  psi_1 <- function(y1,y2,b1,b2,xit1,xit2){
    
    t1 <- J_mat(0,0,b1,b2,xit1,xit2)+J_mat(0,1,b1,b2,xit1,xit2)+J_mat(1,0,b1,b2,xit1,xit2)+J_mat(1,1,b1,b2,xit1,xit2)
  }
  
  final_var_1 <- function(d){
    a <- data.frame(combn(d$y,2))
    y1 <- as.list(a[1,])
    y2 <- as.list(a[2,])   ## responses
    
    w <- exp(beta_0*d[3]+beta_1*d[4]+beta_2*d[5]+beta_3*d[6]+beta_4*d[7]+beta_5*d[8]+beta_6*d[9]+beta_7*d[10]+beta_8*d[11])
    beta <- data.frame(combn(w$x0,2))
    c1 <- as.list(beta[1,])
    c2 <- as.list(beta[2,])  ## values of b
    
    d2 <- data.frame(combn(d$x0,2))
    d3 <- data.frame(combn(d$x1,2))
    d4 <- data.frame(combn(d$x2,2))
    d5 <- data.frame(combn(d$x3,2))
    d6 <- data.frame(combn(d$x4,2))
    d7 <- data.frame(combn(d$x5,2))
    d8 <- data.frame(combn(d$x6,2))
    d9 <- data.frame(combn(d$x7,2))
    d10 <- data.frame(combn(d$x8,2))
    
    make_list_of_vec <- function(a,b,d,e,f,g,h,j,k){   # Will change with covariates
      m <- list()
      for (i in 1:ncol(a)){
        vec <- c(a[1,i],b[1,i],d[1,i],e[1,i],f[1,i],g[1,i],h[1,i],j[1,i],k[1,i])
        m[[i]] <- vec
      }
      m
    }
    make_list_of_vec_1 <- function(a,b,d,e,f,g,h,j,k){   # Will change with covariates
      m <- list()
      for (i in 1:ncol(a)){
        vec <- c(a[2,i],b[2,i],d[2,i],e[2,i],f[2,i],g[2,i],h[2,i],j[2,i],k[2,i])
        m[[i]] <- vec
      }
      m
    }
    
    xit1 <- make_list_of_vec(d2,d3,d4,d5,d6,d7,d8,d9,d10)
    xit2 <- make_list_of_vec_1(d2,d3,d4,d5,d6,d7,d8,d9,d10) ## values of xit
    
    var_3 <- mapply(psi_1,y1,y2,c1,c2,xit1,xit2, SIMPLIFY = FALSE)
    var_4 <- Reduce('+',var_3)
    return(var_4)
  }
  final <- lapply(dat_split,final_var_1)
  return(final)
}









dummy_dat <- data.frame(c(1,1,1,2,2,2,2),c(1,0,0,1,0,1,0),rep(1,7),c(0.2,0.1,0.3,2.3,2,1,0.3),c(0.1,0.2,0.5,0.5,3.2,0.1,1.2))
colnames(dummy_dat) <- c("index","y","x0","x1","x2")

s <- omega(dummy_dat,0,1,0.02,0.07,0.05)
r <- psi(dummy_dat,0,1,0.02,0.07,0.05)
solve(Reduce('+',r)) %*% Reduce('+',s)%*%solve(Reduce('+',r))


variance_3 <- function(data,alpha,beta_0,beta_1,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7,beta_8,rho){
  t1 <- psi(data,alpha,beta_0,beta_1,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7,beta_8,rho) 
  psi_k <- Reduce('+',t1)
  t2 <- omega(data,alpha,beta_0,beta_1,beta_2,beta_3,beta_4,beta_5,beta_6,beta_7,beta_8,rho)
  omega_k <- Reduce('+',t2)
  cov_matrix <- solve(psi_k) %*% omega_k %*% solve(psi_k)
  cov_matrix
}

variance_3(new_dat, 0.1, 6.092,0.012,-0.041,0.553,0.481,0.494,-2.756,-0.554,-0.376,0.431)
g <- psi(new_dat, 0, 5.292,0.018,-0.023,0.418,0.503,0.547,-2.183,-0.380,-0.860,0.422)
h <- omega(new_dat, 0, 5.292,0.018,-0.023,0.418,0.503,0.547,-2.183,-0.380,-0.860,0.422)




