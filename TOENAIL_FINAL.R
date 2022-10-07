install.packages("mice")
library(mice)
library(data.table)
library(Rlab)
library(HDInterval)
library(mvtnorm)
library(matrixStats)
library(VGAM)
library(tibble)






data("toenail")
data <- toenail
data <- subset(data, select = -visit)
colnames(data) <- c("index","y","Ti","tij")
data$Tixtij <- data$Ti*data$tij

num <- c()
for(i in 1:383){
  if((colCounts(as.matrix(data),value = i)[1])==1) {
    num <- append(num,i)
  }
}
data <- subset(data, !(index %in% num))





## Objective function for MLE:

mle_pairwise <- function(x,fixed = c(rep(FALSE,5))){
  params <- fixed
  dpd_1 <- function(p){
    params[!fixed] <- p
    beta_0 <- params[1]
    beta_1 <- params[2]
    beta_2 <- params[3]
    beta_3 <- params[4]
    rho <- params[5]
    add_pi <- function(d){
      k <- beta_0+(d[3]*beta_1)+(d[4]*beta_2)+(d[5]*beta_3)
      k1 <- exp(k)/(1+exp(k))
      d <- cbind(d,k1)
    }
    dat_split <- split(x , f  = x$index)
    result <- lapply(dat_split, add_pi)
    
    result <- rbindlist(result)
    result <- as.data.frame(result)
    colnames(result) <- c('index','y','x1','x2','x3','exp_prob')
    result_split <- split(result, f = result$index)
    
    ml_expression <- function(d){
      bin <- as.data.frame(combn(d$y , 2))
      pr <- as.data.frame(combn(d$exp_prob , 2))
      t <- as.data.frame(combn(d$x2 , 2))
      
      ## Evaluation of the probabilities:
      ml_f_jk <- function(u,v,w){
        
        
        dummy_func <- function(a,b){
          pbern(a, prob = b)
        }
        dummy_func_1 <- function(a,b){
          pbern(a-1, prob = b)
        }
        
        k <- mapply(dummy_func,u,v)
        k_1 <- mapply(dummy_func_1,u,v)
        inv1 <- inverseCDF(as.matrix(k), pnorm)
        inv2 <- inverseCDF(as.matrix(k_1), pnorm)
        mean <- rep(0,2)
        lower <- inv2
        upper <- inv1
        corr <- diag(2)
        corr[lower.tri(corr)] <- rho^(w[2] - w[1])
        corr[upper.tri(corr)] <- rho^(w[2] - w[1])
        probab <- pmvnorm(lower = lower, upper = upper, mean = mean, corr = corr)
        l_prob <- log(abs(probab+0.000001))
        return(l_prob)
      }
      
      final_res <- mapply(ml_f_jk, bin, pr,t)
      
      final_value <- sum(final_res)
    }
    
    u <- sapply(result_split,ml_expression)
    return( -sum(u) )
  }
}

val1 <- mle_pairwise(data,c(rep(FALSE,5)))
### The following IVs have lead to convergence = 0:
optim(c(beta_0 = -0.1,beta_1= 1,beta_2= -0.01,beta_3= -0.001,rho = 0.002),val1, method = "L-BFGS-B", lower = c(-Inf,0.0239,-Inf,-Inf,-1), upper = c(Inf,Inf,Inf,Inf,1), control = list(maxit=500) )

optim(c(beta_0 = -0.1,beta_1 = 0.0027,beta_2= - 0.01,beta_3= -0.001,rho = 0.002), val1 , control = list(maxit=1000))







## DPD objective function: 

dpd_pairwise <- function(x,fixed = c(rep(FALSE,6))){
  params <- fixed
  dpd_1 <- function(p){
    params[!fixed] <- p
    alpha <- params[1]
    beta_0 <- params[2]
    beta_1 <- params[3]
    beta_2 <- params[4]
    beta_3 <- params[5]
    rho <- params[6]
    add_pi <- function(d){
      k <- beta_0+(d[3]*beta_1)+(d[4]*beta_2)+(d[5]*beta_3)
      k1 <- exp(k)/(1+exp(k))
      d <- cbind(d,k1)
    }
    dat_split <- split(x , f  = x$index)
    result <- lapply(dat_split, add_pi)
    
    result <- rbindlist(result)
    result <- as.data.frame(result)
    colnames(result) <- c('index','y','x1','x2','x3','exp_prob')
    result_split <- split(result, f = result$index)
    
    expression <- function(d){
      bin <- as.data.frame(combn(d$y , 2))
      pr <- as.data.frame(combn(d$exp_prob , 2))
      t <- as.data.frame(combn(d$x2 , 2))
      
      ## Evaluation of the probabilities:
      f_jk <- function(u,v,w){
        
        
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
        corr[lower.tri(corr)] <- rho^(w[2] - w[1])
        corr[upper.tri(corr)] <- rho^(w[2] - w[1])
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
          corr1[lower.tri(corr1)] <- rho^(w[2] - w[1])
          corr1[upper.tri(corr1)] <- rho^(w[2] - w[1])
          prob1 <- pmvnorm(lower = lower1, upper = upper1, mean = mean1, corr = corr1)
          prob1 <- prob1^(1+alpha)
          
          val <- append(val,prob1)
        }
        
        val_s <- sum(val)
        
        return(val_s - prob)
      }
      
      final_res <- mapply(f_jk, bin, pr,t)
      
      final_value <- sum(final_res)
    }
    
    u <- sapply(result_split,expression)
    return(sum(u))
  }
}

val <- dpd_pairwise(data,c(1,rep(FALSE,5)))
### The following IVs have lead to convergence = 0:
optim(c(beta_0 = -0.4,beta_1=0.1,beta_2= -0.15,beta_3= -0.07,rho = 0.01),val, method = "L-BFGS-B", lower = c(rep(-Inf, 4),0), upper = c(rep(Inf,4),1) )
optim(c(beta_0 = -0.4,beta_1=0.01,beta_2= -0.15,beta_3= -0.05,rho = 0.02), val )

## for variance:
one <- rep(1,1903)
data_1 <- add_column(data,one, .after =2) ## only for variance
colnames(data_1) <- c("index","y","x0","x1","x2","x3")

#VARIANCE:

## Omega matrix:


omega <- function(dat,alpha,beta_0,beta_1,beta_2,beta_3,rho){
  dat_split <- split(dat, f= dat$index)
  
  xi <- function(y1,y2,b1,b2,xit1,xit2,tim1,tim2){
    rho <- rho^(tim2-tim1) ## modification of rho
    
    
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
  
  k_mat <- function(y1,y2,b1,b2,xit1,xit2,tim1,tim2){
    rho <- rho^(tim2-tim1) ## modification of rho
    
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
  omega_1 <- function(y1,y2,b1,b2,xit1,xit2,tim1,tim2){
    t1 <- xi(0,0,b1,b2,xit1,xit2,tim1,tim2)+xi(0,1,b1,b2,xit1,xit2,tim1,tim2)+xi(1,0,b1,b2,xit1,xit2,tim1,tim2)+xi(1,1,b1,b2,xit1,xit2,tim1,tim2)
    t2 <- t1 %*% t(t1)
    t3 <- k_mat(0,0,b1,b2,xit1,xit2,tim1,tim2)+k_mat(0,1,b1,b2,xit1,xit2,tim1,tim2)+k_mat(1,0,b1,b2,xit1,xit2,tim1,tim2)+k_mat(1,1,b1,b2,xit1,xit2,tim1,tim2)
    return(t3-t2)
  }
  final_var <- function(d){
    a <- data.frame(combn(d$y,2))
    y1 <- as.list(a[1,])
    y2 <- as.list(a[2,])   ## responses
    
    q <- data.frame(combn(d$x2,2))
    tim1 <- as.list(q[1,])
    tim2 <- as.list(q[2,])  ## time points
    
    w <- exp(beta_0*d[3]+beta_1*d[4]+beta_2*d[5]+beta_3*d[6])
    beta <- data.frame(combn(w$x0,2))
    c1 <- as.list(beta[1,])
    c2 <- as.list(beta[2,])  ## values of b
    
    d2 <- data.frame(combn(d$x0,2))
    d3 <- data.frame(combn(d$x1,2))
    d4 <- data.frame(combn(d$x2,2))
    d5 <- data.frame(combn(d$x3,2))
    
    
    make_list_of_vec <- function(a,b,d,e){   # Will change with covariates
      m <- list()
      for (i in 1:ncol(a)){
        vec <- c(a[1,i],b[1,i],d[1,i],e[1,i])
        m[[i]] <- vec
      }
      m
    }
    make_list_of_vec_1 <- function(a,b,d,e){   # Will change with covariates
      m <- list()
      for (i in 1:ncol(a)){
        vec <- c(a[2,i],b[2,i],d[2,i],e[2,i])
        m[[i]] <- vec
      }
      m
    }
    
    xit1 <- make_list_of_vec(d2,d3,d4,d5)
    xit2 <- make_list_of_vec_1(d2,d3,d4,d5) ## values of xit
    
    var_3 <- mapply(omega_1,y1,y2,c1,c2,xit1,xit2,tim1,tim2, SIMPLIFY = FALSE)
    var_4 <- Reduce('+',var_3)
    return(var_4)
  }
  final <- lapply(dat_split,final_var)
}

## PSI matrix:

psi <- function(dat,alpha,beta_0,beta_1,beta_2,beta_3,rho){
  dat_split <- split(dat, f= dat$index)
  J_mat <- function(y1,y2,b1,b2,xit1,xit2,tim1,tim2){
    rho <- rho^(tim2-tim1) ## modification of rho
    
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
  psi_1 <- function(y1,y2,b1,b2,xit1,xit2,tim1,tim2){
    
    t1 <- J_mat(0,0,b1,b2,xit1,xit2,tim1,tim2)+J_mat(0,1,b1,b2,xit1,xit2,tim1,tim2)+J_mat(1,0,b1,b2,xit1,xit2,tim1,tim2)+J_mat(1,1,b1,b2,xit1,xit2,tim1,tim2)
  }
  
  final_var_1 <- function(d){
    a <- data.frame(combn(d$y,2))
    y1 <- as.list(a[1,])
    y2 <- as.list(a[2,])   ## responses
    
    q <- data.frame(combn(d$x2,2))
    tim1 <- as.list(q[1,])
    tim2 <- as.list(q[2,])  ## time points
    
    
    w <- exp(beta_0*d[3]+beta_1*d[4]+beta_2*d[5]+beta_3*d[6])
    beta <- data.frame(combn(w$x0,2))
    c1 <- as.list(beta[1,])
    c2 <- as.list(beta[2,])  ## values of b
    
    d2 <- data.frame(combn(d$x0,2))
    d3 <- data.frame(combn(d$x1,2))
    d4 <- data.frame(combn(d$x2,2))
    d5 <- data.frame(combn(d$x3,2))
    
    
    make_list_of_vec <- function(a,b,d,e){   # Will change with covariates
      m <- list()
      for (i in 1:ncol(a)){
        vec <- c(a[1,i],b[1,i],d[1,i],e[1,i])
        m[[i]] <- vec
      }
      m
    }
    make_list_of_vec_1 <- function(a,b,d,e){   # Will change with covariates
      m <- list()
      for (i in 1:ncol(a)){
        vec <- c(a[2,i],b[2,i],d[2,i],e[2,i])
        m[[i]] <- vec
      }
      m
    }
    
    xit1 <- make_list_of_vec(d2,d3,d4,d5)
    xit2 <- make_list_of_vec_1(d2,d3,d4,d5) ## values of xit
    
    var_3 <- mapply(psi_1,y1,y2,c1,c2,xit1,xit2,tim1,tim2, SIMPLIFY = FALSE)
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


variance_3 <- function(data,alpha,beta_0,beta_1,beta_2,beta_3,rho){
  t1 <- psi(data,alpha,beta_0,beta_1,beta_2,beta_3,rho) 
  psi_k <- Reduce('+',t1)
  t2 <- omega(data,alpha,beta_0,beta_1,beta_2,beta_3,rho)
  omega_k <- Reduce('+',t2)
  cov_matrix <- solve(psi_k) %*% omega_k %*% solve(psi_k)
  cov_matrix
}

variance_3(data_1, 1, -0.476,0.029,-0.211,-0.073,0.946)
variance_3(data_1, 0, -0.553,0.024,-0.179,-0.064,0.937)
variance_3(data_1, 0.1, -0.535,0.010,-0.183,-0.065,0.937)
variance_3(data_1, 0.3, -0.525,0.026,-0.188,-0.071,0.938)
variance_3(data_1, 0.5, -0.510,0.032,-0.195,-0.075,0.939)
variance_3(data_1, 0.7, -0.492,0.035,-0.204,-0.077,0.942)



g <- psi(data_1, 0, -0.6018, 0.0239,-0.1669,-0.0583,0.7636)
h <- omega(data_1, 0, -0.6018, 0.0239,-0.1669,-0.0583,0.7636)

