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

## Final form to obtain the variance:

variance_3 <- function(data,alpha,beta_0,beta_1,beta_2,beta_3,rho){
  t1 <- psi(data,alpha,beta_0,beta_1,beta_2,beta_3,rho) 
  psi_k <- Reduce('+',t1)
  t2 <- omega(data,alpha,beta_0,beta_1,beta_2,beta_3,rho)
  omega_k <- Reduce('+',t2)
  cov_matrix <- solve(psi_k) %*% omega_k %*% solve(psi_k)
  cov_matrix
}


