set.seed(12)
# =====================================
# Dataset generating function
# =====================================
rdat <- function(n=100, beta =c(-0.7, -0.7, 0, 0, 0, -0.7, 0, 0, 0),
                 rho=0.5, censor.control=1)
{
  p <- length(beta)
  # GENERATE X
  mu <- rep(0, p)
  S <- matrix(1, p, p)
  for (i in 1:p){
    for (j in 1:p){
      S[i, j] <- rho^(abs(i-j))
    }
  }
  X <- mvrnorm(n = n, mu=mu, Sigma=S, tol=1e-6, empirical=F)
  # GENERATE SURVIVAL TIME AND CENSORING TIME
  rate <- exp(X%*%beta)
  T0 <- rexp(n, rate)
  C0 <- rexp(n, rate)
  time <- pmin(T0, (C0*censor.control))
  status <- sign(T0 <= (C0*censor.control))
  # OUTPUT THE DATA
  dat <- data.frame(cbind(id=1:n, time, status, X, T0, C0))
  colnames(dat) <- c("id", "time", "status", paste("x", 1:p, sep=""),
                     "true.time", "true.censor")
  return(list(dat=dat, beta.true=beta, S=S))
}

# =====================================
# Group Data Generation
# =====================================

genData <- function(nSim=100, n=100, beta=c(-0.7, -0.7, 0, 0, 0, -0.7, 0, 0, 0),rho=0.5,censor.control=1)
{
  gennedData <- list()
  for (i in 1:nSim) {
    gennedData[[i]] <- list(rdat(n,beta,rho,censor.control))
  }
  return(gennedData)
}

# =====================================
# Fit Lasso using glmnet
# =====================================
fitLasso <- function(z, tsurv)
{
  lasso_cv <- cv.glmnet(x = z, y = tsurv, family = "cox", nfolds = 10, alpha = 1)
  lasso_l1se <- lasso_cv$lambda.1se
  lasso_fit <- glmnet:::glmnet(x = z, y = tsurv, family = "cox", lambda = lasso_l1se, alpha = 1)
  return(as.matrix(lasso_fit$beta, ncol = p, nrow = 1))
}

# =====================================
# Fit Ridge using glmnet
# =====================================
fitRidge <- function(z, tsurv)
{
  ridge_cv <- cv.glmnet(x = z, y = tsurv, family = "cox", nfolds = 10, alpha = 0)
  ridge_l1se <- ridge_cv$lambda.1se
  ridge_fit <- glmnet:::glmnet(x = z, y = tsurv, family = "cox", lambda = ridge_l1se, alpha = 0)
  return(as.matrix(ridge_fit$beta,ncol=p,nrow=1))
}

# =====================================
# Fit Elastic Net using glmnet
# =====================================
fitEnet <- function(z, tsurv, a = 0.5)
{
  enet_cv <- cv.glmnet(x = z, y = tsurv, family = "cox", nfolds = 10, alpha = a)
  enet_l1se <- enet_cv$lambda.1se
  enet_fit <- glmnet:::glmnet(x = z, y = tsurv, family = "cox", lambda = enet_l1se, alpha = a)
  return(as.matrix(enet_fit$beta, ncol = p, nrow = 1))
}

# =====================================
# Fit Relaxed LASSO using glmnet
# =====================================
fitRlasso <- function(z, tsurv)
{
  rlasso_cv <- cv.glmnet(x = z, y = tsurv, family = "cox", nfolds = 10, alpha = 1, relax = TRUE)
  l1se <- rlasso_cv$lambda.1se
  rlasso_fit <- glmnet:::glmnet(x = z, y = tsurv, family = "cox", relax = TRUE, lambda = l1se, alpha = 1)
  return(as.matrix(rlasso_fit$beta, ncol = p, nrow = 1))
}

# =====================================
# Fit Adaptive LASSO using glmnet
# =====================================
fitAlasso <- function(z, tsurv)
{
  # STAGE 1: initial ridge regression
  ridge0_cv <- cv.glmnet(x = z, y = tsurv, family = "cox", nfolds = 10, alpha = 0)
  ridge0_l1se <- ridge0_cv$lambda.1se
  ridge0_fit <- glmnet:::glmnet(x = z, y = tsurv, family = "cox", lambda = ridge0_l1se, alpha = 0)
  ridge0_beta <- as.matrix(ridge0_fit$beta, ncol = p, nrow = 1)
  penalty <- 1/abs(ridge0_beta)
  
  # STAGE 2: adaptively weighted LASSO regression
  alasso_cv <- cv.glmnet(x = z, y = tsurv, family = "cox", nfolds = 10, alpha = 1)
  alasso_l1se <- alasso_cv$lambda.1se
  alasso_fit <- glmnet:::glmnet(x = z, y = tsurv, family = "cox", lambda = alasso_l1se, alpha = 1, penalty.factor = penalty)
  
  return(as.matrix(alasso_fit$beta, ncol = p, nrow = 1))
}

# =====================================
# Fit Stepwise Regression using My.stepwise
# =====================================
stepwiseRegression <- function(Time=NULL, T1=NULL, T2=NULL, Status=NULL, variable.list, in.variable="NULL", data, sle=0.15, sls=0.15)
{
  
  univar.pvalue <- NULL
  temp.model <- NULL
  
  if (is.null(T2)) {
    initial.model <- coxph(as.formula(paste("Surv(", Time,", ", Status,") ~ ", paste(in.variable, collapse="+"), sep="")), data=data, method="efron")
  } else if (is.null(Time)){
    initial.model <- coxph(as.formula(paste("Surv(", T1,", ", T2,", ", Status,") ~ ", paste(in.variable, collapse="+"), sep="")), data=data, method="efron")
  }
  
  if (is.null(initial.model$coefficients)) {
    
    for (i in 1:length(variable.list))
    {
      if (is.null(T2))
      {
        uni.model <- coxph(as.formula(paste("Surv(", Time,", ", Status,") ~ ", variable.list[i], sep="")), data=data, method="efron")
      }
      if (is.null(Time))
      {
        uni.model <- coxph(as.formula(paste("Surv(", T1,", ", T2,", ", Status,") ~ ", variable.list[i], sep="")), data=data, method="efron")
      }
      univar.pvalue[i] <- summary(uni.model)$coefficients[5]
    }
    
    variable.list1 <- variable.list[univar.pvalue<=0.9 & !is.na(univar.pvalue)]
    univar.pvalue1 <- univar.pvalue[univar.pvalue<=0.9 & !is.na(univar.pvalue)]
    uni.x <- variable.list1[which.min(univar.pvalue1)]
    if (length(uni.x) > 0) {
      if (is.null(T2))
      {
        formula <- as.formula(paste("Surv(", Time,", ", Status,") ~ ", uni.x, sep=""))
        temp.model <- coxph(formula, data=data, method="efron")
      }
      if (is.null(Time))
      {
        formula <- as.formula(paste("Surv(", T1,", ", T2,", ", Status,") ~ ", uni.x, sep=""))
        temp.model <- coxph(formula, data=data, method="efron")
      }
    }
    
  } else if (!is.null(initial.model$coefficients)) {
    temp.model <- initial.model
  }
  
  i <- 0
  break.rule <- TRUE
  while (break.rule)
  {
    i <- i + 1
    if (i == 1)
    {
      variable.list2 <- setdiff(variable.list, all.vars(temp.model$formula))
    } else
    {
      variable.list2 <- setdiff(variable.list, c(all.vars(temp.model$formula), out.x))
      out.x <- NULL
    }
    
    if (length(variable.list2) != 0)
    {
      anova.pvalue <- NULL
      mv.pvalue <- NULL
      for (k in 1:length(variable.list2))
      {
        model <- update(temp.model, as.formula(paste(". ~ . + ", variable.list2[k], sep="")))
        if (length(model$coefficients) > 1)
        {
          if (sum(is.na(model$coefficients)) != 0)
          {
            anova.pvalue[k] <- 1
            mv.pvalue[k] <- 1
          } else {
            anova.pvalue[k] <- anova(temp.model, model)[2,"P(>|Chi|)"]
            mv.pvalue[k] <- summary(model)$coefficients[nrow(summary(model)$coefficients),"Pr(>|z|)"]
          }
        }
      }
      
      variable.list2.1 <- variable.list2[mv.pvalue<=0.9 & !is.na(mv.pvalue)]
      anova.pvalue2 <- anova.pvalue[mv.pvalue<=0.9 & !is.na(mv.pvalue)]
      mv.pvalue2 <- mv.pvalue[mv.pvalue<=0.9 & !is.na(mv.pvalue)]
      enter.x <- variable.list2.1[anova.pvalue2==min(anova.pvalue2, na.rm=TRUE) & anova.pvalue2 <= sle]
      wald.p <- mv.pvalue2[anova.pvalue2==min(anova.pvalue2, na.rm=TRUE) & anova.pvalue2 <= sle]
      if (length(setdiff(enter.x, NA)) != 0)
      {
        if (length(enter.x) > 1)
        {
          enter.x <- enter.x[which.min(wald.p)]
        }
        temp.model <- update(temp.model, as.formula(paste(". ~ . + ", enter.x, sep="")))
      }
    } else {enter.x <- NULL}
    
    if (i != 1 || length(enter.x) != 0) {
      variable.list3 <- setdiff(rownames(summary(temp.model)$coefficients), c(enter.x, in.variable))
      if (length(variable.list3) != 0)
      {
        anova.pvalue <- NULL
        for (k in 1:length(variable.list3))
        {
          model <- update(temp.model, as.formula(paste(". ~ . - ", variable.list3[k], sep="")))
          anova.pvalue[k] <- anova(model, temp.model)[2,"P(>|Chi|)"]
        }
        
        out.x <- variable.list3[anova.pvalue==max(anova.pvalue, na.rm=TRUE) & anova.pvalue > sls]
        out.x <- setdiff(out.x, NA)
        if (length(out.x) != 0)
        {
          if (length(out.x) > 1)
          {
            out.x.1 <- out.x
            for (j in 1:length(out.x)) {
              out.x[j] <- out.x.1[(length(out.x)-j+1)]
            }
            
            wald.p <- rep(NA, length(out.x))
            for (j in 1:length(out.x)) {
              wald.p[j] <- summary(temp.model)$coefficients[,"Pr(>|z|)"][rownames(summary(temp.model)$coefficients)==out.x[j]]
            }
            out.x <- out.x[which.max(wald.p)]
          }
          temp.model <- update(temp.model, as.formula(paste(". ~ . - ", out.x, sep="")))
        }
      } else {out.x <- NULL}
    }
    
    if ((length(enter.x) + length(out.x)) == 0)
    {
      return(temp.model)
      break.rule <- FALSE
    }
    
    enter.x <- NULL
  }
}

fitStepwise <- function(dataset, sle = 0.15, sls = 0.15) 
{
  temp_coef <- coef(stepwiseRegression(Time = "time", Status = "status", variable.list = c("x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9"),
                                data = dataset, sle = sle, sls = sls))
  return(c(temp_coef['x1'], temp_coef['x2'], temp_coef['x3'], temp_coef['x4'], temp_coef['x5'], temp_coef['x6'], temp_coef['x7'], temp_coef['x8'], temp_coef['x9']))
}

# =====================================
# Fit Models
# =====================================

fitModels <- function(dataset, nSim, n, beta, p)
{
  lasso_coef <- data.frame(array(0, dim = c(nSim, p)))
  ridge_coef <- data.frame(array(0, dim = c(nSim, p)))
  enet_coef <- data.frame(array(0, dim = c(nSim, p)))
  rlasso_coef <- data.frame(array(0, dim = c(nSim, p)))
  alasso_coef <- data.frame(array(0, dim = c(nSim, p)))
  #step_coef <- data.frame(array(0, dim = c(nSim, p)))
  scad_coef <- data.frame(array(0, dim = c(nSim, p)))
  mcp_coef <- data.frame(array(0, dim = c(nSim, p)))
  for (i in 1:nSim){
    
    # break up the data in frame i into covariates, event time and censoring indicator
    z <- as.matrix(dataset[[i]][[1]]$dat[,4:12], ncol=p, nrow=n)
    t <- dataset[[i]][[1]]$dat[,2]
    delta <- dataset[[i]][[1]]$dat[,3]
    
    # create survival object for glmnet
    tsurv <- survival:::Surv(t, delta)
    
    # Estimate coefficients for glmnet functions
    lasso_coef[i,] <- fitLasso(z = z, tsurv = tsurv)
    ridge_coef[i,] <- fitRidge(z = z, tsurv = tsurv)
    enet_coef[i,] <- fitEnet(z = z, tsurv = tsurv, a = 0.5)
    rlasso_coef[i,] <- fitRlasso(z = z, tsurv = tsurv)
    alasso_coef[i,] <- fitAlasso(z = z, tsurv = tsurv)
    
    # Estimate coefficients for the simplified stepwise model
    # step_coef[i,] <- fitStepwise(data = dataset[[i]][[1]]$dat, sle = 0.25, sls = 0.25)
    
    # Estimate coefficients for the nonconvex penalty functions
    # Same cross validation lambda is used in both methods
    NCV_cv <- cv.ncvsurv(X = z, y = tsurv, nfolds = 10)
    NCV_lmin <- NCV_cv$lambda.min
    
    # Fit the models using the cv-lambda
    scad_coef[i,] <- as.matrix(coef(ncvsurv(X = z, y = tsurv, penalty = "SCAD"), lambda = NCV_lmin), ncol = p, nrow = 1)
    mcp_coef[i,] <- as.matrix(coef(ncvsurv(X = z, y = tsurv, penalty = "MCP"), lambda = NCV_lmin), ncol = p, nrow = 1)
    print(i)
  }
  #step_coef[is.na(step_coef)] <- 0
  return(list("lasso" = lasso_coef, "ridge" = ridge_coef, "enet" = enet_coef, "rlasso" = rlasso_coef, "alasso" = alasso_coef, "scad" = scad_coef, "mcp" = mcp_coef))
}

# =====================================
# Calculate Squared Error 
# =====================================

calcMSE <- function(df, beta_true) {
  df <- t(t(df)-beta)^2
  MSE <- rowSums(df)
  avMSE <- mean(MSE)
  meMSE <- median(MSE)
  return(list("avMSE" = avMSE, "meMSE" = meMSE))
}

# =====================================
# Calculate Success Rate
# =====================================

calcSR <- function(df, beta_true) {
  trials <- NROW(df)
  for(i in 1:length(beta_true)) {
    if(NROW(df) == 0) {
      break
    } else if(beta_true[i] == 0) {
      df <- df[df[,i]==0,]
    } else {
      df <- df[df[,i]!=0,]
    }
  }
  return(NROW(df)/trials)
}

# =====================================
# Calculate Inclusion and Exclusion Errors
# =====================================
calcIncE <- function(df, beta_true) {
  inc <- 0
  incC <- 0
  incI <- 0
  exc <- 0
  excC <- 0
  excI <- 0
  for(i in 1:length(beta_true)) {
    if(beta_true[i] == 0) {
      excC <- excC + NROW(df[df[,i]==0,])
      incI <- incI + NROW(df[df[,i]!=0,])
      exc <- exc + 1
    } else {
      incC <- incC + NROW(df[df[,i]!=0,])
      excI <- excI + NROW(df[df[,i]==0,])
      inc <- inc + 1
    }
  }
  s <- NROW(df)*length(beta_true)
  return(list("iI" = incI/s, "cI" = incC/s, "iE" = excI/s, "cE" = excC/s, "c" = (incC + excC)/s, "i" = (excI + incI)/s))
}

simFit <- function(nSim, n, beta, p, censor.control, rho) {

}

# =====================================
# Setting Simulation Variables
# =====================================

simulating <- FALSE
if(simulating) {
  output <- data.frame("Method"=NA, "Success Rate"=NA, "Average Model Size"=NA, "Mean MSE"=NA, "Median MSE"=NA, "Correct Inclusion"=NA, "Incorrect Inclusion"=NA, "Correct Exclusion"=NA, "Incorrect Exclusion"=NA)
  nSim <- 100 #number of simulations to perform
  n <- 150
  beta <- c(-0.4, -0.3, 0, 0, 0, -0.2, -0.2, 0, 0)
  p <- length(beta)
  censor.control <- 1
  rho <- 0.2
  simData <- genData(nSim = nSim, n = n, beta = beta, rho = rho, censor.control = censor.control)
  beta_est <- fitModels(dataset = simData, nSim = nSim, n = n, beta = beta, p = p)
  
  for (i in 1:length(beta_est)) {
    
    methodname <- switch(i, "lasso", "ridge", "enet", "rlasso", "alasso", "scad", "mcp")
    MSE <- calcMSE(beta_est[[i]], beta)
    incRates <- calcIncE(beta_est[[i]], beta)
    ssR <- calcSR(beta_est[[i]], beta)
    output[NROW(output) + 1,] <- c(methodname, ssR, mean(rowSums(beta_est[[i]]!=0)), MSE["avMSE"], MSE["meMSE"], incRates["cI"], incRates["iI"], incRates["cE"], incRates["iE"])
    
  }
  
  output[,-1] <- round(output[,-1],2)
  output <- output[-1,]
  write.table(output, paste("Documents/Projects/MSCi\ project/Results/Scheme-2.n-", n), sep = "\t", row.names = FALSE)
}
