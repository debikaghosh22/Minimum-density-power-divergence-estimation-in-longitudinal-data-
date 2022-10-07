## MLE on contaminated data:


datalist_cont <- list()
for (j in 1:600){
  set.seed(j)
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
    prob <- c()
    y_1 <- rbinom(1,1,p[1,2])
    y <- append(y,y_1)
    prob <- append(prob,p[1,2])
    for(i in 2:6){
      pr <- p[i,2] + 0.4*(y[i-1]-p[i-1,2])
      u <- rbinom(1,1,pr)
      prob <- append(prob,pr)
      y <- append(y,u)
    }
    df <- data.frame(y,prob)
    return(df)
  }
  
  res <- lapply(data_split,create_y)
  res <- rbindlist(res)
  res <- as.data.frame(res)
  ##y <- data.frame(x=unlist(res)) 
  o <- which(res$y == 0)[order(res$prob[res$y == 0])]
  o <- o[1:60]
  
  for (i in o) {
    res[i,1] <- ifelse(res[i,1]==0 , 1,res[i,1])
  }
  y1 <- res[,1]
  
  
  final_data <- data.frame(index,y1,x1,x2)
  colnames(final_data) <- c("index","y","x1","x2")
  datalist_cont[[j]] <- final_data
}


## Objective function:
mle_pairwise <- function(x,fixed = c(rep(FALSE,4))){
  params <- fixed
  dpd_1 <- function(p){
    params[!fixed] <- p
    beta_0 <- params[1]
    beta_1 <- params[2]
    beta_2 <- params[3]
    rho <- params[4]
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
    
    ml_expression <- function(d){
      bin <- as.data.frame(combn(d$y , 2))
      pr <- as.data.frame(combn(d$exp_prob , 2))
      
      ## Evaluation of the probabilities:
      ml_f_jk <- function(u,v){
        
        
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
        corr[lower.tri(corr)] <- rho
        corr[upper.tri(corr)] <- rho
        probab <- pmvnorm(lower = lower, upper = upper, mean = mean, corr = corr)
        l_prob <- log(abs(probab))
        return(l_prob)
      }
      
      final_res <- mapply(ml_f_jk, bin, pr)
      
      final_value <- sum(final_res)
    }
    
    u <- sapply(result_split,ml_expression)
    return(-sum(u))
  }
}
#Simulation: (contaminated data)

cl = makeCluster(39)
registerDoParallel(cl)


m_cont <- foreach(i = 1:600, .packages = c('data.table','Rlab','HDInterval','mvtnorm'), .errorhandling = 'remove',.combine = rbind) %dopar% {
  beta <- rbind(1,0.02,0.05,0.4)
  val <- mle_pairwise(datalist_cont[[i]], c(FALSE,FALSE,FALSE,FALSE))
  b_s <- as.vector(optim(c(beta_0 =0.9, beta_1 =0.005 ,beta_2 = 0.009,rho=0.1),val)$par)
  conv <- optim(c(beta_0 =0.9, beta_1 =0.005 ,beta_2 = 0.009,rho = 0.1),val)$convergence
  k <- (b_s -beta)
  k <- append(k,conv)
}
stopCluster(cl)


