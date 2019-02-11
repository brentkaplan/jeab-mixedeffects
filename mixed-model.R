#### Mixed-Effects Model on data from Study One of Ackerlund Brandt et al. (2015) ####

#### Authors: William Brady DeHart & Brent A. Kaplan ####

## Check directory
pltdir <- "plots/"
if (!dir.exists(pltdir)) dir.create(pltdir)

## Load libraries
## For any packages not downloaded, please use install.packages("package_name")
library(tidyverse) # Nice package for data manipulation/preparation. Includes ggplot for graphs
library(broom) # Data manipulation
library(ggpubr)
library(lme4) # Needed for R2 calculation
library(psych)
library(glmmTMB) # Package for zero-inflated mixed-model
library(DHARMa) # glmmTMB model diagnostics
library(sjstats) # Random Effects ICC
library(vcdExtra) # Zero-inflation test


#### Functions ####
## https://stackoverflow.com/questions/16249570/uppercase-the-first-letter-in-data-frame
capFirst <- function(s) {
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

#### Read in data ####
Final <- read.csv("data/extracted-data.csv", header = T, sep = ",")
Final <- Final %>%
  select(child, condition, session, frequency)

# Re-order session so that control is the comparison group
Final$child <- capFirst(Final$child)
Final$child <- as.factor(Final$child)
Final$condition <- capFirst(Final$condition)
Final$condition <- factor(Final$condition, 
                          levels = c("Control", "Experimenter", "Child"))
Final$Condition <- Final$condition


#### Plot of data by subject ####
png(paste0(pltdir, "Single Subject.png"), width = 20, height = 12, units = "in", res = 300)
ggplot(Final, aes(x = session, y = frequency, group = Condition)) + 
  geom_line(size = 1) +
  geom_point(aes(shape = Condition, fill = Condition), size = 5, stroke = 1.5) +
  scale_shape_manual(values = c(22, 21, 21)) + 
  scale_fill_manual(values = c("white", "white", "black")) +
  scale_x_continuous(limits = c(1,14), breaks = c(seq(1,14,2)), 
                     labels = c( "1", "3", "5", "7", "9", "11", "13")) +
  scale_y_continuous(limits = c(-1,17), 
                     breaks = c(0, 5, 10, 15), 
                     labels = c("0", "5", "10", "15")) +
  theme_bw(base_size = 28) +
  ylab("Selections") +
  xlab("Session") +
  ggtitle("Single-Subject Results") +
  geom_text(data=Final,x=14.5,y=16.5,aes(label=child),
            vjust = "inward", hjust = "inward", size = 7)+
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black", size = 1), 
        plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(size = 1), 
        axis.ticks = element_line(size = 1), 
        strip.background = element_blank(), 
        strip.text = element_blank()) +
  facet_wrap(~child, ncol = 3)
dev.off()


####Data Exploration####

###Test for zero-inflation###
zero.test(Final$frequency)

###Overdispersion###
var(Final$frequency) / mean(Final$frequency)


#### Mixed-Model ####

### Model Fit
## Note that the following code may result in model estimates that are within
## a range of ~.05 depending on R version/operating system
model.1 <- glmmTMB(frequency ~ condition*I(log(session)) + (condition + I(log(session))| child), 
                   data = Final, ziformula = ~ condition*I(log(session)),  
                   family = nbinom2(), REML = FALSE)
summary(model.1)
fixedeff <- glmmTMB::fixef(model.1)
exp(fixedeff$cond)

###Custom Function for Model Fit R2

##' extract the 'conditional-model' term from a glmmTMB object;
##' otherwise, return x unchanged
collapse_cond <- function(x) {
  if (is.list(x) && "cond" %in% names(x)) x[["cond"]] else x
  
}

##' Cleaned-up/adapted version of Jon Lefcheck's code from SEMfit;
##' also incorporates some stuff from MuMIn::rsquaredGLMM.
##' Computes Nakagawa/Schielzeth/Johnson analogue of R^2 for
##' GLMMs. Should work for [g]lmer(.nb), glmmTMB models ...
##'
##' @param model a fitted model
##' @return a list composed of elements "family", "link", "marginal", "conditional"
my_rsq <- function(model) {
  
  ## get basics from model (as generally as possible)
  vals <- list(
    beta=fixef(model),
    X=getME(model,"X"),
    vc=VarCorr(model),
    re=ranef(model))
  
  ## glmmTMB-safety
  if (is(model,"glmmTMB")) {
    vals <- lapply(vals,collapse_cond)
    nullEnv <- function(x) {
      environment(x) <- NULL
      return(x)
    }
    if (!identical(nullEnv(model$modelInfo$allForm$ziformula),nullEnv(~0)))
      warning("R2 ignores effects of zero-inflation")
    dform <- nullEnv(model$modelInfo$allForm$dispformula)
    if (!identical(dform,nullEnv(~1)) &&
        (!identical(dform,nullEnv(~0))))
      warning("R2 ignores effects of dispersion model")
  }
  
  ## Test for non-zero random effects
  if (any(sapply(vals$vc, function(x) any(diag(x)==0)))) {
    ## FIXME: test more generally for singularity, via theta?
    stop("Some variance components equal zero. Respecify random structure!")
  }
  
  ## set family/link info
  ret <- list()
  if (is(model,"glmmTMB") || is(model,"glmerMod")) {
    ret$family <- family(model)$family
    ret$link <- family(model)$link
  } else {
    ret$family <- "gaussian"; ret$link <- "identity"
  }
  
  ## Get variance of fixed effects: multiply coefs by design matrix
  varF <- with(vals,var(as.vector(beta %*% t(X))))
  
  ## Are random slopes present as fixed effects? Warn.
  random.slopes <- if("list" %in% class(vals$re)) {
    ## multiple RE
    unique(c(sapply(vals$re,colnames)))
  } else {
    colnames(vals$re)
  }
  if (!all(random.slopes %in% names(vals$beta))) 
    warning("Random slopes not present as fixed effects. This artificially inflates the conditional R2. Respecify fixed structure!")
  
  ## Separate observation variance from variance of random effects
  nr <- sapply(vals$re, nrow)
  not.obs.terms <- names(nr[nr != nobs(model)])
  obs.terms <- names(nr[nr==nobs(model)])
  
  ## Compute variance associated with a random-effects term
  ## (Johnson 2014)
  getVarRand <- function(terms) {
    sum(
      sapply(vals$vc[terms],
             function(Sigma) {
               Z <- vals$X[, rownames(Sigma), drop = FALSE]
               Z.m <- Z %*% Sigma
               return(sum(diag(crossprod(Z.m, Z))) / nobs(model))
             } )
    )
  }
  
  ## Variance of random effects 
  varRand <- getVarRand(not.obs.terms)
  
  if (is(model,"lmerMod") ||
      (ret$family=="gaussian" && ret$link=="identity")) {
    ## Get residual variance
    varDist <- sigma(model)^2
    varDisp <- 0
  } else {
    varDisp <- if (length(obs.terms)==0) 0 else getVarRand(obs.terms)
    
    badlink <- function(link,family) {
      warning(sprintf("Model link '%s' is not yet supported for the %s distribution",link,family))
      return(NA)
    }
    
    if(ret$family == "binomial") {
      varDist <- switch(ret$link,
                        logit=pi^2/3,
                        probit=1,
                        badlink(ret$link,ret$family))
    } else if (ret$family == "poisson" ||
               grepl("nbinom",ret$family) ||
               grepl("Negative Binomial", ret$family)) {
      ## Generate null model (intercept and random effects only, no fixed effects)
      
      ## https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q4/023013.html
      ## FIXME: deparse is a *little* dangerous
      rterms <- paste0("(",sapply(findbars(formula(model)),deparse),")")
      nullform <- reformulate(rterms,response=".")
      null.model <- update(model,nullform)
      
      ## from MuMIn::rsquaredGLMM
      
      ## Get the fixed effects of the null model
      null.fixef <- unname(collapse_cond(fixef(null.model)))
      
      ## in general want log(1+var(x)/mu^2)
      logVarDist <- function(null.fixef) {
        mu <- exp(null.fixef)
        if (mu < 6)
          warning(sprintf("mu of %0.1f is too close to zero, estimate may be unreliable \n",mu))
        vv <- switch(ret$family,
                     poisson=mu,
                     nbinom1=,
                     nbinom2=family(model)$variance(mu,sigma(model)),
                     if (is(model,"merMod"))
                       mu*(1+mu/getME(model,"glmer.nb.theta"))
                     else mu*(1+mu/model$theta))
        cvsquared <- vv/mu^2
        return(log1p(cvsquared))
      }
      
      varDist <- switch(ret$link,
                        log=logVarDist(null.fixef),
                        sqrt=0.25,
                        badlink(ret$link,ret$family))
    }
  }
  ## Calculate R2 values
  ret$Marginal = varF / (varF + varRand + varDisp + varDist)
  ret$Conditional = (varF + varRand) / (varF + varRand + varDisp + varDist)
  return(ret)
}

my_rsq(model.1)

###Random Effects ICC
sjstats::icc(model.1)

###Simulate Model Fit###
res <- simulateResiduals(model.1, refit = F, n = 1000)

###Test to see if model accounts for zero inflation
testZeroInflation(res)

###Test to see if model accounts for overdispersion
testDispersion(res)

####Diagnostic Plots####

###Plots###

##Test of Linearity##
lindf <- data.frame(Final$frequency, res$fittedPredictedResponse)
Plot1 <- ggplot(lindf, aes(x = Final.frequency, y = res.fittedPredictedResponse)) +
  geom_point(size = 4) +
  stat_smooth(method = "loess", size = 1.5, color = "black", lty = "dashed") +
  scale_x_continuous(limits = c(0,15), breaks = c(0, 3, 6, 9, 12, 15), 
                     labels = c("0", "3", "6", "9", "12", "15")) +
  scale_y_continuous(limits = c(-2,15), breaks = c(0, 5, 10, 15), 
                     labels = c("0", "5", "10", "15")) +
  xlab("Frequency") +
  ylab("Fitted") +
  ggtitle("Test of Linearity") +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(size = 1), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(hjust = 0.5))


##Plot Residuals##
modelResid <- res$scaledResiduals
predicted <- res$fittedPredictedResponse
predicted <- rank(predicted, ties.method = "average")
predicted <- predicted/max(predicted)
residfull <- data.frame(modelResid,predicted)
Plot2 <- ggplot(residfull, aes(x = predicted, y = modelResid)) +
  geom_point(size = 4, alpha = .7) +
  theme_bw(base_size = 24) +
  ggtitle("Ranked Residuals vs Predicted Plot") +
  xlab("Predicted Value") +
  ylab("Residuals") +
  geom_hline(yintercept = .25, color = "black", size = 1.5) +
  geom_hline(yintercept = .50, color = "black", size = 1.5) +
  geom_hline(yintercept = .75, color = "black", size = 1.5) +
  geom_quantile(color = "black", linetype = "dashed", size = 2) +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(size = 1), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(hjust = 0.5))


##QQplot##
n <- length(res$scaledResiduals)
m <- (1:n)/(n + 1)
qqdataframe <- data.frame(m, res$scaledResiduals)
sx <- sort(qqdataframe$m)
sy <- sort(qqdataframe$res.scaledResiduals)
lenx <- length(sx)
leny <- length(sy)
if (leny < lenx) 
  sx <- approx(1L:lenx, sx, n = leny)$y
if (leny > lenx) 
  sy <- approx(1L:leny, sy, n = lenx)$y
Plot3 <- ggplot(qqdataframe, aes(x = sx, y = sy)) +
  geom_point(size = 4, shape = 1, stroke = 1.5) +
  xlab("Expected") +
  ylab("Observed") +
  xlim(0,1) +
  ylim(0,1) +
  geom_abline(intercept = 0, size = 1.5, linetype = "dashed") +
  ggtitle("Q-Q Plot Residuals") + 
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(size = 1), 
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(hjust = 0.5))


##Histogram of Residuals##
histdf <- data.frame(res$fittedResiduals)
Plot4 <- ggplot(histdf, aes(x = res$fittedResiduals)) + 
  geom_histogram(color = "black", fill = "lightgray", bins = 12) +
  scale_x_continuous(limits = c(-6,6), breaks = c( -6, -4, -2, 0, 2, 4, 6)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 150)) +
  theme_bw(base_size = 24) +
  ggtitle("Histogram of Residuals") +
  xlab("Residuals") +
  ylab("Count") +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        panel.border = element_blank(), 
        axis.line = element_line(size = 1),
        axis.ticks = element_line(size = 1), 
        plot.title = element_text(hjust = 0.5))

png(paste0(pltdir, "Figure 5.png"), width = 20, height = 12, units = "in", res = 300)
ggarrange(Plot1, Plot2, Plot3, Plot4)
dev.off()

#### Prediction Graphs ####

### Extract predicted values
Final2 <- Final %>% 
  mutate(pred_dist = fitted(model.1)) 
png(paste0(pltdir, "Figure 6.png"), width = 20, height = 12, units = "in", res = 300)
ggplot(Final2, aes(x=session, y=pred_dist, group=Condition)) + 
  geom_smooth(aes(group = Condition, linetype = Condition), method='lm', 
              fullrange = TRUE, color = "black", size = 1.5, se=FALSE) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  scale_x_continuous(limits = c(1,14), breaks = c(seq(1,14,2)), 
                     labels = c( "1", "3", "5", "7", "9", "11", "13")) +
  scale_y_continuous(limits = c(-2,18), breaks = c(0, 5, 10, 15), 
                     labels = c("0", "5", "10", "15")) +
  geom_text(data=Final,x=.75,y=18,aes(label=child),
            vjust = "inward", hjust = "inward", size = 7)+
  xlab("Session") +
  ylab("Predicted Behavior") +
  ggtitle("Individual Participant Model Predictions") +
  theme_bw(base_size = 24) +
  theme(axis.text = element_text(color = "black"), 
        panel.grid = element_blank(), 
        axis.line = element_line(color = "black", size = 1), 
        plot.title = element_text(hjust = 0.5), 
        panel.border = element_rect(size = 1), 
        axis.ticks = element_line(size = 1), 
        strip.background = element_blank(), 
        strip.text = element_blank(), 
        legend.key.width = unit(3, "line")) +
  facet_wrap(~child, ncol = 3)
dev.off()



