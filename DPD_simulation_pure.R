library(foreach)
library(doParallel)
library(data.table)
library(Rlab)
library(HDInterval)
library(mvtnorm)
library(matrixStats)
# DATA
datalist <- list()
for (i in 1:600){
  set.seed(i)
  x1 <- rnorm(600,1.5,2)
  x2 <- rnorm(600,0,1)
  data_2 <- data.frame(x1,x2)
  lin_pred <- 1+(0.02 * data_2[1]) + (0.05*data_2[2]) 
  prob <- 1/(1+exp(-lin_pred))
  index <- rep(1:100, each = 6)
  data <- data.frame(index,prob)
  data_split <- split(data,f=data$index)
  create_y <- function(p){
    
    y <- c()
    y_1 <- rbinom(1,1,p[1,2])
    y <- append(y,y_1)
    for(i in 2:6){
      pr <- p[i,2] + 0.4*(y[i-1]-p[i-1,2])
      u <- rbinom(1,1,pr)
      y <- append(y,u)
    }
    return(y)
  }
  
  res <- lapply(data_split,create_y)
  y <- data.frame(x=unlist(res))
  final_data <- data.frame(index,y,x1,x2)
  datalist[[i]] <- final_data
}
## Objective function:

dpd_pairwise <- function(x,fixed = c(rep(FALSE,5))){
  params <- fixed
  dpd_1 <- function(p){
    params[!fixed] <- p
    alpha <- params[1]
    beta_0 <- params[2]
    beta_1 <- params[3]
    beta_2 <- params[4]
    rho <- params[5]
    add_pi <- function(d){
      k <- beta_0+(d[3]*beta_1)+(d[4]*beta_2)
      k1 <- exp(k)/(1+exp(k))
      d <- cbind(d,k1)
    }
    dat_split <- split(x , f  = x$index)
    result <- lapply(dat_split, add_pi)
    
    result <- rbindlist(result)
    result <- as.data.frame(result)
    colnames(result) <- c('index','y','x1','x2','exp_prob')
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
#Simulation: (pure data)

cl = makeCluster(39)
registerDoParallel(cl)



mat_1 <- foreach(i = 1:600, .packages = c('data.table','Rlab','HDInterval','mvtnorm'), .errorhandling = 'remove',.combine = rbind) %dopar% {
  beta <- rbind(1,0.02,0.05,0.4)
  val <- dpd_pairwise(datalist[[i]], c(1,FALSE,FALSE,FALSE,FALSE))
  b_s <- as.vector(optim(c(beta_0 =0.9, beta_1 =0.005 ,beta_2 = 0.009,rho=0.1),val)$par)
  conv <- optim(c(beta_0 =0.9, beta_1 =0.005 ,beta_2 = 0.009,rho = 0.1),val)$convergence
  k <- (b_s -beta)
  k <- append(k,conv)
}
stopCluster(cl)


m1 <- head(mat_0.1,501)
m1 <- subset(m1,m1[,5] ==0) ## for alpha = 0.1
m2 <- head(mat_0.3,501)
m2 <- subset(m2,m2[,5] ==0)
m3 <- head(mat_0.5,503)
m3 <- subset(m3,m3[,5] ==0)


m4 <- head(mat_0.7,501)
m4 <- subset(m4,m4[,5] ==0)

m5 <- head(mat_1,505)
m5 <- subset(m5,m5[,5] ==0)

m9 <- head(mat_0.7_cont,503)
m9 <- subset(m9,m9[,5] ==0)


m10 <- head(mat_1_cont,505)
m10 <- subset(m10,m10[,5] ==0)
