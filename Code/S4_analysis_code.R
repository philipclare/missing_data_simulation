
args <- commandArgs(trailingOnly = TRUE)
# args <- c(1,10,"C:/Users/z3312911/Cloudstor/PhD/Katana/missingsim/")
start <- as.numeric(args[1])
stop <- as.numeric(args[2])
n <- stop-start+1
nimp <- 40

packages <- paste0("/home/z3312911/RPackages/")
# packages <- paste0("C:/Users/z3312911/Cloudstor/R Library")
.libPaths(packages)

load(paste0(args[3],"seeds.RData"))
set.seed(eval(seeds[as.numeric(stop/n)]))

filepath1 <- paste0(args[3],"Data/S1/")
filepath2 <- paste0(args[3],"Data/S2/")
filepath3 <- paste0(args[3],"Data/S3/")

## Check required packages are installed, and if not, install them
libs <- c("nnls","SuperLearner","Matrix","foreach","glmnet","ranger","lme4","geepack","parallel",
          "doParallel","mvtnorm","survival","TH.data","MASS","splines","boot","haven","ggplot2",
          "multcomp","doBy","gam","future","stats","data.table","optimx","mitml","Amelia","mice","norm")
missing <- !libs %in% installed.packages()
if (any(missing)) {
  install.packages(libs[missing],repos="https://cloud.r-project.org")
}

library("nnls")
library("SuperLearner")
library("Matrix")
library("foreach")
library("glmnet")
library("ranger")
library("lme4")
library("geepack")
library("parallel")
library("doParallel")
library("boot")
library("haven")
library("ggplot2")
library("multcomp")
library("doBy")
library("gam")
library("future")
library("stats")
library("data.table")
library("optimx")
library("plyr")
library("ltmle")
library("mitml")
library("mice")
library("norm")
library("Amelia")
library("clubSandwich")

## Define custom RF wrapper with slightly fewer trees (based on previous work, this is all that is required for this data)
create.Learner("SL.ranger", params = list(num.trees = 250))

## Program to truncate inverse probability weights
## Used to reduce variability in weights, which can cause issues
checkrange <- function(v) {
  v <- ifelse(v<0.001,0.001,v)
  v <- ifelse(v>0.999,ifelse(is.na(v),v,0.999),v)
}

## Complete-case analysis syntax
simcomp <- function(data) {
  data<-data[complete.cases(data),]
  
  # Define analytic models for outcome, propensities, etc
  outmodel <- "y ~ a + la + l + ll + oa + ob + oc + obs"
  outmodel2 <- "y ~ a + la + ll + oa + ob + oc + obs"
  
  propa0model <- "la ~ ll + oa + ob + oc + obs"
  propa1model <- "a ~ la + l + ll + oa + ob + oc + obs"
  
  # qform/gform is essentially the same, but in the form required by the package ltmle
  qform <- c(l="Q.kplus1 ~ la + ll + obs + oa + ob + oc",
             y="Q.kplus1 ~ a + la + l + la + oa + ob + oc + obs")
  gform1 <- c(la="la ~ ll + oa + ob + oc + obs",
              a="a ~ la + l + ll + oa + ob + oc + obs")
  
  create.Learner("SL.ranger", params = list(num.trees = 250)) # same as in top-level code - replicated because sometimes goes wrong in parallel
  
  # Define libraries to be used by SuperLearner
  # Top set is set to be used in final analysis. Second set is used for testing.
  SLlib1 <- c("SL.glm","SL.glm.interaction","SL.gam","SL.ranger_1")
  SLlib2 <- list(g=c("SL.glm","SL.glm.interaction","SL.gam","SL.ranger_1"),
                 Q=c("SL.glm","SL.glm.interaction","SL.gam"))
  
  SLlib1 <- c("SL.glm","SL.glm.interaction")
  SLlib2 <- list(g=c("SL.glm","SL.glm.interaction"),
                 Q=c("SL.glm","SL.glm.interaction"))
  
  # GLM IPTW-based analysis
  GLexp0 <- glm(formula=propa0model,data=data,family=binomial)
  GLexp1 <- glm(formula=propa1model,data=data,family=binomial)
  GLexpp <- data.table(cbind(id=data$id,obs=data$obs,a=data$a,la=data$la,
                             propa=ifelse(data$a==1,checkrange(predict(GLexp1,type="response")),checkrange(1-predict(GLexp1,type="response"))),
                             propla=ifelse(data$la==1,checkrange(predict(GLexp0,type="response")),checkrange(1-predict(GLexp0,type="response")))))
  GLexpp$p <- GLexpp$propla*GLexpp$propa
  GLexpp$GLwt <- 1/GLexpp$p
  
  # GLM IPT and DR-IPT analysis
  GLiptw <- glm(y~a+la,data=merge(data,GLexpp[,c("id","obs","GLwt")]),weight=GLwt,family=gaussian)
  GLdriptw <- glm(outmodel2,data=merge(data,GLexpp[,c("id","obs","GLwt")]),weight=GLwt,family=gaussian) 
  
  # SuperLearner IPTW-based analysis
  SLexp0 <- SuperLearner(Y=as.vector(data[,]$la),X=data[,c("obs","oa","ob","oc","ll")],id=data[,1],SL.library=SLlib1,family=binomial)
  SLexp1 <- SuperLearner(Y=as.vector(data[,]$a),X=data[,c("obs","oa","ob","oc","ll","la","l")],id=data[,1],SL.library=SLlib1,family=binomial)
  SLexpp <- data.table(cbind(id=data$id,obs=data$obs,a=data$a,la=data$la,
                             propa=ifelse(data$a==1,checkrange(predict(SLexp1)$pred),checkrange(1-predict(SLexp1)$pred)),
                             propla=ifelse(data$la==1,checkrange(predict(SLexp0)$pred),checkrange(1-predict(SLexp0)$pred))))
  SLexpp$p <- SLexpp$propla*SLexpp$propa
  SLexpp$SLwt <- 1/SLexpp$p
  
  # SL IPT and DR-IPT analysis
  SLiptw <- glm(y~a+la,data=merge(data,SLexpp[,c("id","obs","SLwt")]),weight=SLwt,family=gaussian)
  SLdriptw <- glm(outmodel2,data=merge(data,SLexpp[,c("id","obs","SLwt")]),weight=SLwt,family=gaussian)
  
  # GLM TMLE
  GLtmle <- ltmle(data[,c("obs","oa","ob","oc","ll","la","l","a","y")],
                  id=data[,1],
                  Anodes=c("la","a"),
                  Lnodes=c("l"),
                  Ynodes="y",
                  Qform=qform,
                  gform=gform1,
                  abar=list(c(1,1),c(0,0)),
                  estimate.time = FALSE)
  
  # SuperLearner TMLE    
  SLtmle <- ltmle(data[,c("obs","oa","ob","oc","ll","la","l","a","y")],
                  id=data[,1],
                  Anodes=c("la","a"),
                  Lnodes=c("l"),
                  Ynodes="y",
                  Qform=qform,
                  gform=gform1,
                  abar=list(c(1,1),c(0,0)),
                  SL.library = SLlib2,
                  estimate.time = FALSE)
  
  # Naive analysis
  GLM <- glm(outmodel,data=data,family=gaussian)
  RI <- lmer(paste0(outmodel,"+(1|id)"),data=data)
  GEE <- geeglm(formula(outmodel),data=data,id=data$id,waves=data$obs,family=gaussian)
  
  c(coef(summary(GLM))[2,1]+coef(summary(GLM))[3,1],sqrt(vcovCR(GLM, cluster=data$id, type = "CR3")[2,2] + vcovCR(GLM, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(GLM, cluster=data$id, type = "CR3")[2,3]),
    coef(summary(RI))[2,1]+coef(summary(RI))[3,1],sqrt(vcov(RI)[2,2] + vcov(RI)[3,3] + 2*vcov(RI)[2,3]),
    coef(summary(GEE))[2,1]+coef(summary(GEE))[3,1],sqrt(summary(GEE)$cov.scaled[2,2] + summary(GEE)$cov.scaled[3,3] + 2*summary(GEE)$cov.scaled[2,3]),
    coef(summary(GLiptw))[2,1]+coef(summary(GLiptw))[3,1],sqrt(vcovCR(GLiptw, cluster=data$id, type = "CR3")[2,2] + vcovCR(GLiptw, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(GLiptw, cluster=data$id, type = "CR3")[2,3]),
    coef(summary(SLiptw))[2,1]+coef(summary(SLiptw))[3,1],sqrt(vcovCR(SLiptw, cluster=data$id, type = "CR3")[2,2] + vcovCR(SLiptw, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(SLiptw, cluster=data$id, type = "CR3")[2,3]),
    coef(summary(GLdriptw))[2,1]+coef(summary(GLdriptw))[3,1],sqrt(vcovCR(GLdriptw, cluster=data$id, type = "CR3")[2,2] + vcovCR(GLdriptw, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(GLdriptw, cluster=data$id, type = "CR3")[2,3]),
    coef(summary(SLdriptw))[2,1]+coef(summary(SLdriptw))[3,1],sqrt(vcovCR(SLdriptw, cluster=data$id, type = "CR3")[2,2] + vcovCR(SLdriptw, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(SLdriptw, cluster=data$id, type = "CR3")[2,3]),
    summary(GLtmle)$effect.measures$ATE$estimate,summary(GLtmle)$effect.measures$ATE$std.dev,
    summary(SLtmle)$effect.measures$ATE$estimate,summary(SLtmle)$effect.measures$ATE$std.dev)
}

## Missing Data analysis syntax
simmiss <- function(data,nimpute,scenario) {
  
  if (scenario==1) {
    # Define analytic models for outcome, propensities, etc
    outmodel <- "y ~ a + la + l + ll + oa + ob + oc + obs"
    outmodel2 <- "y ~ a + la + ll + oa + ob + oc + obs"
    
    propa0model <- "la ~ ll + oa + ob + oc + obs"
    propa1model <- "a ~ la + l + ll + oa + ob + oc + obs"
    
    missmodel <- "c ~ a + la + l + ll + oa + ob + oc + obs"
    
    # qform/gform is essentially the same, but in the form required by the package ltmle
    qform <- c(l="Q.kplus1 ~ la + ll + obs + oa + ob + oc",
               y="Q.kplus1 ~ a + la + l + la + oa + ob + oc + obs")
    gform1 <- c(la="la ~ ll + oa + ob + oc + obs",
                a="a ~ la + l + ll + oa + ob + oc + obs")
    gform2 <- c(la="la ~ ll + oa + ob + oc + obs",
                a="a ~ la + l + ll + oa + ob + oc + obs",
                c="c ~ a + la + l + ll + oa + ob + oc + obs")
    
    # Define imputation formulae for jomo impute
    impfmlc <- y + a ~ la + l + ll + oa + ob + oc + obs
    impfml <- y + a ~ l + ll + oa + ob + oc + obs + (1|id)
    
    # Reshape data into wide format and drop all lagged variables except those at baseline
    datawide <- data
    datawide$id <- as.numeric(levels(datawide$id))[datawide$id]
    datawide <- reshape(datawide[,c(-10)],
                        v.names=c("ll","la","l","a","y"),
                        idvar="id",
                        timevar="obs",
                        direction="wide")[,c(-10,-11,-15,-16,-20,-21)]
    
    # Standard-FCS imputation via MICE
    # Run imputation and save as list of data frames
    impdatafcs <- complete(mice(datawide[,-1],
                                m=nimpute,
                                maxit=10),
                           action="all")
    
    # Reshape back into long form for analysis, included re-creating lagged variables
    impdatafcs <- lapply(impdatafcs, function (x) {
      x <- cbind(id=seq.int(nrow(x)),x)
      x$ll.2 <- x$l.1
      x$ll.3 <- x$l.2
      x$ll.4 <- x$l.3
      x$la.2 <- x$a.1
      x$la.3 <- x$a.2
      x$la.4 <- x$a.3
      x <- data.frame(x[,c(1:9,19,22,10:12,20,23,13:15,21,24,16:18)])
      x <- reshape(x,varying=list(c(5,10,15,20),
                                  c(6,11,16,21),
                                  c(7,12,17,22),
                                  c(8,13,18,23),
                                  c(9,14,19,24)),
                   v.names=c("ll","la","l","a","y"),
                   times=c(1,2,3,4),
                   sep=".",
                   idvar="id",
                   timevar="obs",
                   direction="long")
      x <- data.frame(x)
    })
    
    # JM imputation via MVN
    # Run imputation and save as list of data frames
    s <- prelim.norm(as.matrix(datawide[,-1]))
    thetahat <- em.norm(s)
    rngseed(269012)
    theta <- da.norm(s, thetahat, steps=100, showits=TRUE)
    ximp <- lapply(vector("list", length=nimpute),function (x) {data.frame(imp.norm(s,theta,as.matrix(datawide[,-1])))})
    # Reshape back into long form for analysis, included re-creating lagged variables
    impdatamvn <- lapply(ximp, function (x) {
      x <- cbind(id=seq.int(nrow(x)),x)
      x$ll.2 <- x$l.1
      x$ll.3 <- x$l.2
      x$ll.4 <- x$l.3
      x$la.2 <- x$a.1
      x$la.3 <- x$a.2
      x$la.4 <- x$a.3
      x <- data.frame(x[,c(1:9,19,22,10:12,20,23,13:15,21,24,16:18)])
      x <- reshape(x,varying=list(c(5,10,15,20),
                                  c(6,11,16,21),
                                  c(7,12,17,22),
                                  c(8,13,18,23),
                                  c(9,14,19,24)),
                   v.names=c("ll","la","l","a","y"),
                   times=c(1,2,3,4),
                   sep=".",
                   idvar="id",
                   timevar="obs",
                   direction="long")
      x <- data.frame(x)
    })
    
    # Random-intercept imputation via jomo
    # Transform binary variables into factors (required for JOMO to run)
    data$ll <- factor(data$ll)
    data$la <- factor(data$la)
    data$l <- factor(data$l)
    data$a <- factor(data$a)
    # Run imputation and save as list of data frames
    impdatajomo <- mitmlComplete(jomoImpute(data=data.frame(data[,!(names(data) %in% "c")]),formula=impfml,m=nimpute,n.burn=100,n.iter=10),print="all")
    data$ll <- as.numeric(levels(data$ll))[data$ll]
    data$la <- as.numeric(levels(data$la))[data$la]
    data$l <- as.numeric(levels(data$l))[data$l]
    data$a <- as.numeric(levels(data$a))[data$a]
    
  } else {
    # Define analytic models for outcome, propensities, etc
    outmodel <- "y ~ a + la + l + ll + oa + ob + oc + obs"
    outmodel2 <- "y ~ a + la + ll + oa + ob + oc + obs"
    
    propa0model <- "la ~ ll + oa + ob + oc + obs"
    propa1model <- "a ~ la + l + ll + oa + ob + oc + obs"
    
    missmodel <- "c ~ a + la + w + lw + oa + ob + oc + obs"
    
    # qform/gform is essentially the same, but in the form required by the package ltmle
    qform <- c(l="Q.kplus1 ~ la + ll + obs + oa + ob + oc",
               y="Q.kplus1 ~ a + la + l + la + oa + ob + oc + obs")
    gform1 <- c(la="la ~ ll + oa + ob + oc + obs",
                a="a ~ la + l + ll + oa + ob + oc + obs")
    gform2 <- c(la="la ~ ll + oa + ob + oc + obs",
                a="a ~ la + l + ll + oa + ob + oc + obs",
                c="c ~ a + la + w + lw + oa + ob + oc + obs")
    
    # Define imputation formulae for jomo impute
    impfmlc <- y + a ~ la + l + ll + w + lw + oa + ob + oc + obs
    impfml <- y + a ~ l + ll + w + lw + oa + ob + oc + obs + (1|id)
    
    # Reshape data into wide format and drop all lagged variables except those at baseline
    datawide <- data
    datawide$id <- as.numeric(levels(datawide$id))[datawide$id]
    datawide <- reshape(datawide[,c(-12)],
                        v.names=c("ll","lw","la","l","w","a","y"),
                        idvar="id",
                        timevar="obs",
                        direction="wide")[,c(-12,-13,-14,-19,-20,-21,-26,-27,-28)]
    
    # Standard-FCS imputation via MICE
    # Run imputation and save as list of data frames
    impdatafcs <- complete(mice(datawide[,-1],
                                m=nimpute,
                                maxit=10),
                           action="all")
    
    # Reshape back into long form for analysis, included re-creating lagged variables
    impdatafcs <- lapply(impdatafcs, function (x) {
      x <- cbind(id=seq.int(nrow(x)),x)
      x$ll.2 <- x$l.1
      x$ll.3 <- x$l.2
      x$ll.4 <- x$l.3
      x$lw.2 <- x$w.1
      x$lw.3 <- x$w.2
      x$lw.4 <- x$w.3
      x$la.2 <- x$a.1
      x$la.3 <- x$a.2
      x$la.4 <- x$a.3
      x <- data.frame(x[,c(1:11,24,27,30,12:15,25,28,31,16:19,26,29,32,20:23)])
      x <- reshape(x,varying=list(c(5,12,19,26),
                                  c(6,13,20,27),
                                  c(7,14,21,28),
                                  c(8,15,22,29),
                                  c(9,16,23,30),
                                  c(10,17,24,31),
                                  c(11,18,25,32)),
                   v.names=c("ll","lw","la","l","w","a","y"),
                   times=c(1,2,3,4),
                   sep=".",
                   idvar="id",
                   timevar="obs",
                   direction="long")
      x <- data.frame(x)
    })
    
    # Random-intercept imputation via MVN
    # Run imputation and save as list of data frames
    s <- prelim.norm(as.matrix(datawide[,-1]))
    thetahat <- em.norm(s)
    rngseed(269012)
    theta <- da.norm(s, thetahat, steps=100, showits=TRUE)
    ximp <- lapply(vector("list", length=nimpute),function (x) {data.frame(imp.norm(s,theta,as.matrix(datawide[,-1])))})
    # Reshape back into long form for analysis, included re-creating lagged variables
    impdatamvn <- lapply(ximp, function (x) {
      x <- cbind(id=seq.int(nrow(x)),x)
      x$ll.2 <- x$l.1
      x$ll.3 <- x$l.2
      x$ll.4 <- x$l.3
      x$lw.2 <- x$w.1
      x$lw.3 <- x$w.2
      x$lw.4 <- x$w.3
      x$la.2 <- x$a.1
      x$la.3 <- x$a.2
      x$la.4 <- x$a.3
      x <- data.frame(x[,c(1:11,24,27,30,12:15,25,28,31,16:19,26,29,32,20:23)])
      x <- reshape(x,varying=list(c(5,12,19,26),
                                  c(6,13,20,27),
                                  c(7,14,21,28),
                                  c(8,15,22,29),
                                  c(9,16,23,30),
                                  c(10,17,24,31),
                                  c(11,18,25,32)),
                   v.names=c("ll","lw","la","l","w","a","y"),
                   times=c(1,2,3,4),
                   sep=".",
                   idvar="id",
                   timevar="obs",
                   direction="long")
      x <- data.frame(x)
    })
    
    # Random-intercept imputation via jomo
    # Transform binary variables into factors (required for JOMO to run)
    data$ll <- factor(data$ll)
    data$la <- factor(data$la)
    data$lw <- factor(data$la)
    data$l <- factor(data$l)
    data$a <- factor(data$a)
    data$w <- factor(data$a)
    # Run imputation and save as list of data frames
    impdatajomo <- mitmlComplete(jomoImpute(data=data.frame(data[,!(names(data) %in% "c")]),formula=impfml,m=nimpute,n.burn=100,n.iter=10),print="all")
    data$ll <- as.numeric(levels(data$ll))[data$ll]
    data$la <- as.numeric(levels(data$la))[data$la]
    data$l <- as.numeric(levels(data$l))[data$l]
    data$a <- as.numeric(levels(data$a))[data$a]
  }
  
  create.Learner("SL.ranger", params = list(num.trees = 250)) # same as in top-level code - replicated because sometimes goes wrong in parallel
  
  # Define libraries to be used by SuperLearner
  # Top set is set to be used in final analysis. Second set is used for testing.
  SLlib1 <- c("SL.glm","SL.glm.interaction","SL.gam","SL.ranger_1")
  SLlib2 <- list(g=c("SL.glm","SL.glm.interaction","SL.gam","SL.ranger_1"),
                 Q=c("SL.glm","SL.glm.interaction","SL.gam"))
  
  SLlib1 <- c("SL.glm","SL.glm.interaction")
  SLlib2 <- list(g=c("SL.glm","SL.glm.interaction"),
                 Q=c("SL.glm","SL.glm.interaction"))
  
  # GLM IPTW-based analysis
  GLexp0 <- glm(formula=propa0model,data=data,family=binomial)
  GLexp1 <- glm(formula=propa1model,data=data,family=binomial)
  GLexpc <- glm(formula=missmodel,data=data,family=binomial)
  GLexpp <- data.table(cbind(id=data$id,obs=data$obs,a=data$a,la=data$la,
                             propa=ifelse(data$a==1,checkrange(predict(GLexp1,type="response")),checkrange(1-predict(GLexp1,type="response"))),
                             propla=ifelse(data$la==1,checkrange(predict(GLexp0,type="response")),checkrange(1-predict(GLexp0,type="response"))),
                             propc=checkrange(1-predict(GLexpc,type="response"))))
  GLexpp$p <- GLexpp$propla*GLexpp$propa
  GLexpp$pc <- GLexpp$propla*GLexpp$propa*GLexpp$propc
  GLexpp$GLwt <- 1/GLexpp$p
  GLexpp$GLwtc <- 1/GLexpp$pc
  GLexpp$GLwc <- 1/GLexpp$propc
  
  # IPC-weighted GLM IPT and DR-IPT analysis
  GLiptcw <- glm(y~a+la,data=merge(data,GLexpp[,c("id","obs","GLwtc")]),weight=GLwtc,family=gaussian)
  GLdriptcw <- glm(outmodel2,data=merge(data,GLexpp[,c("id","obs","GLwtc")]),weight=GLwtc,family=gaussian) 
  
  # SuperLearner IPTW-based analysis
  SLexp0 <- SuperLearner(Y=as.vector(data[,]$la),X=data[,c("obs","oa","ob","oc","ll")],id=data[,1],SL.library=SLlib1,family=binomial)
  SLexp1 <- SuperLearner(Y=as.vector(data[,]$a),X=data[,c("obs","oa","ob","oc","ll","la","l")],id=data[,1],SL.library=SLlib1,family=binomial)
  SLexpc <- SuperLearner(Y=as.vector(data$c),X=data[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a") else c("obs","oa","ob","oc","lw","la","w","a")],id=data[,1],SL.library=SLlib1,family=binomial)
  SLexpp <- data.table(cbind(id=data$id,obs=data$obs,a=data$a,la=data$la,
                             propa=ifelse(data$a==1,checkrange(predict(SLexp1)$pred),checkrange(1-predict(SLexp1)$pred)),
                             propla=ifelse(data$la==1,checkrange(predict(SLexp0)$pred),checkrange(1-predict(SLexp0)$pred)),
                             propc=as.vector(checkrange(1-predict(SLexpc)$pred))))
  SLexpp$p <- SLexpp$propla*SLexpp$propa
  SLexpp$pc <- SLexpp$propla*SLexpp$propa*SLexpp$propc
  SLexpp$SLwt <- 1/SLexpp$p
  SLexpp$SLwtc <- 1/SLexpp$pc
  SLexpp$SLwc <- 1/SLexpp$propc
  
  # IPC-weighted SL IPT and DR-IPT analysis
  SLiptcw <- glm(y~a+la,data=merge(data,SLexpp[,c("id","obs","SLwtc")]),weight=SLwtc,family=gaussian)
  SLdriptcw <- glm(outmodel2,data=merge(data,SLexpp[,c("id","obs","SLwtc")]),weight=SLwtc,family=gaussian)
  
  data$c <- BinaryToCensoring(is.censored=data$c)
  # IPC-weighted GLM TMLE
  GLtmlec <- ltmle(data[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a","c","y") else c("obs","oa","ob","oc","ll","lw","la","l","w","a","c","y")],
                   id=data[,"id"],
                   Anodes=c("la","a"),
                   Lnodes=if (scenario==1) c("l") else c("l","w"),
                   Ynodes="y",
                   Cnodes="c",
                   Qform=qform,
                   gform=gform2,
                   abar=list(c(1,1),c(0,0)),
                   estimate.time = FALSE)
  
  # IPC-weighted SuperLearner TMLE
  SLtmlec <- ltmle(data[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a","c","y") else c("obs","oa","ob","oc","ll","lw","la","l","w","a","c","y")],
                   id=data[,"id"],
                   Anodes=c("la","a"),
                   Lnodes=if (scenario==1) c("l") else c("l","w"),
                   Ynodes="y",
                   Cnodes="c",
                   Qform=qform,
                   gform=gform2,
                   abar=list(c(1,1),c(0,0)),
                   SL.library = SLlib2,
                   estimate.time = FALSE)
  
  # GLM IPC-weighted naive analysis
  GLGLMc <- glm(outmodel,data=merge(data,GLexpp[,c("id","obs","GLwc")]),weight=GLwc,family=gaussian)
  GLRIc <- lmer(paste0(outmodel,"+(1|id)"),data=merge(data,GLexpp[,c("id","obs","GLwc")]),weight=GLwc)
  GLGEEc <- geeglm(formula(outmodel),data=merge(data,GLexpp[,c("id","obs","GLwc")]),id=data$id,waves=data$obs,weight=GLwc,family=gaussian)
  
  # SL IPC-weighted naive analysis
  SLGLMc <- glm(outmodel,data=merge(data,SLexpp[,c("id","obs","SLwc")]),weight=SLwc,family=gaussian)
  SLRIc <- lmer(paste0(outmodel,"+(1|id)"),data=merge(data,SLexpp[,c("id","obs","SLwc")]),weight=SLwc)
  SLGEEc <- geeglm(formula(outmodel),data=merge(data,SLexpp[,c("id","obs","SLwc")]),id=data$id,waves=data$obs,weight=SLwc,family=gaussian)
  
  # FCS Multiple Imputation Analyses
  # Imputed naive analysis
  # GLM
  GLMimpfcs <- lapply(impdatafcs, function (x) {glm(formula=formula(outmodel),data=x,family=gaussian)})
  GLMimpcofcs <- matrix(unlist(lapply(GLMimpfcs, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLMimpsefcs <- matrix(unlist(lapply(GLMimpfcs, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLMfcs <- matrix(unlist(mi.meld(q=GLMimpcofcs, se=GLMimpsefcs)),nrow=1,ncol=2)
  
  # Random intercept
  RIimpfcs <- lapply(impdatafcs, function (x) {lmer(paste0(outmodel,"+(1|id)"),data=x)})
  RIimpcofcs <- matrix(unlist(lapply(RIimpfcs, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  RIimpsefcs <- matrix(unlist(lapply(RIimpfcs, function (x) {sqrt(vcov(x)[2,2] + vcov(x)[3,3] + 2*vcov(x)[2,3])})),nrow=nimpute,ncol=1)
  RIfcs <- matrix(unlist(mi.meld(q=RIimpcofcs, se=RIimpsefcs)),nrow=1,ncol=2)
  
  # GEE
  GEEimpfcs <- lapply(impdatafcs, function (x) {geeglm(formula(outmodel),data=x,id=x$id,waves=x$obs,family=gaussian)})
  GEEimpcofcs <- matrix(unlist(lapply(GEEimpfcs, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GEEimpsefcs <- matrix(unlist(lapply(GEEimpfcs, function (x) {sqrt(summary(x)$cov.scaled[2,2] + summary(x)$cov.scaled[3,3] + 2*summary(x)$cov.scaled[2,3])})),nrow=nimpute,ncol=1)
  GEEfcs <- matrix(unlist(mi.meld(q=GEEimpcofcs, se=GEEimpsefcs)),nrow=1,ncol=2)
  
  # GLM-based IPTW and DR-IPTW
  GLiptwimpfcs <- lapply(impdatafcs, function (x) {
    GLexp0 <- glm(formula=propa0model,data=x,family=binomial)
    GLexp1 <- glm(formula=propa1model,data=x,family=binomial)
    GLexpp <- data.table(cbind(id=x$id,obs=x$obs,a=x$a,la=x$la,
                               propa=ifelse(x$a==1,checkrange(predict(GLexp1,type="response")),checkrange(1-predict(GLexp1,type="response"))),
                               propla=ifelse(x$la==1,checkrange(predict(GLexp0,type="response")),checkrange(1-predict(GLexp0,type="response")))))
    GLexpp$p <- GLexpp$propla*GLexpp$propa
    GLexpp$GLwt <- 1/GLexpp$p
    D <- merge(x,GLexpp[,c("id","obs","GLwt")])
    D
  })
  GLiptwimpsumfcs <- lapply(GLiptwimpfcs, function (x) {glm(y~a+la,data=x,weight=GLwt,family=gaussian)})
  GLiptwimpcofcs <- matrix(unlist(lapply(GLiptwimpsumfcs, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLiptwimpsefcs <- matrix(unlist(lapply(GLiptwimpsumfcs, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLiptwfcs <- matrix(unlist(mi.meld(q=GLiptwimpcofcs, se=GLiptwimpsefcs)),nrow=1,ncol=2)
  
  GLdriptwimpsumfcs <- lapply(GLiptwimpfcs, function (x) {glm(outmodel2,data=x,weight=GLwt,family=gaussian)})
  GLdriptwimpcofcs <- matrix(unlist(lapply(GLdriptwimpsumfcs, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLdriptwimpsefcs <- matrix(unlist(lapply(GLdriptwimpsumfcs, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLdriptwfcs <- matrix(unlist(mi.meld(q=GLdriptwimpcofcs, se=GLdriptwimpsefcs)),nrow=1,ncol=2)
  
  # SuperLearner-based IPTW and DR-IPTW
  SLiptwimpfcs <- lapply(impdatafcs, function (x) {
    SLexp0 <- SuperLearner(Y=as.numeric(as.vector(x$la)),X=x[,c("obs","oa","ob","oc","ll")],id=x[,"id"],SL.library=SLlib1,family=binomial)
    SLexp1 <- SuperLearner(Y=as.numeric(as.vector(x$a)),X=x[,c("obs","oa","ob","oc","ll","la","l")],id=x[,"id"],SL.library=SLlib1,family=binomial)
    SLexpp <- data.table(cbind(id=x$id,obs=x$obs,a=x$a,la=x$la,
                               propa=ifelse(x$a==1,checkrange(predict(SLexp1)$pred),checkrange(1-predict(SLexp1)$pred)),
                               propla=ifelse(x$la==1,checkrange(predict(SLexp0)$pred),checkrange(1-predict(SLexp0)$pred))))
    SLexpp$p <- SLexpp$propla*SLexpp$propa
    SLexpp$SLwt <- 1/SLexpp$p
    D <- merge(x,SLexpp[,c("id","obs","SLwt")])
    D
  })
  SLiptwimpsumfcs <- lapply(SLiptwimpfcs, function (x) {glm(y~a+la,data=x,weight=SLwt,family=gaussian)})
  SLiptwimpcofcs <- matrix(unlist(lapply(SLiptwimpsumfcs, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  SLiptwimpsefcs <- matrix(unlist(lapply(SLiptwimpsumfcs, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  SLiptwfcs <- matrix(unlist(mi.meld(q=SLiptwimpcofcs, se=SLiptwimpsefcs)),nrow=1,ncol=2)
  
  SLdriptwimpsumfcs <- lapply(SLiptwimpfcs, function (x) {glm(outmodel2,data=x,weight=SLwt,family=gaussian)})
  SLdriptwimpcofcs <- matrix(unlist(lapply(SLdriptwimpsumfcs, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  SLdriptwimpsefcs <- matrix(unlist(lapply(SLdriptwimpsumfcs, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  SLdriptwfcs <- matrix(unlist(mi.meld(q=SLdriptwimpcofcs, se=SLdriptwimpsefcs)),nrow=1,ncol=2)
  
  # GLM-based TMLE
  GLtmleimpfcs <- lapply(impdatafcs, function (x) {
    ltmle(x[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a","y") else c("obs","oa","ob","oc","ll","lw","la","l","w","a","y")],
          id=x[,"id"],
          Anodes=c("la","a"),
          Lnodes=if (scenario==1) c("l") else c("l","w"),
          Ynodes="y",
          Qform=qform,
          gform=gform1,
          abar=list(c(1,1),c(0,0)),
          estimate.time = FALSE)})
  GLtmleimpsumfcs <- lapply(GLtmleimpfcs, function (x) {summary(x)})
  GLtmleimpcofcs <- matrix(unlist(lapply(GLtmleimpsumfcs, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
  GLtmleimpsefcs <- matrix(unlist(lapply(GLtmleimpsumfcs, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
  GLtmlefcs <- matrix(unlist(mi.meld(q=GLtmleimpcofcs, se=GLtmleimpsefcs)),nrow=1,ncol=2)
  
  # SuperLearner based TMLE
  SLtmleimpfcs <- lapply(impdatafcs, function (x) {
    ltmle(x[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a","y") else c("obs","oa","ob","oc","ll","lw","la","l","w","a","y")],
          id=x[,"id"],
          Anodes=c("la","a"),
          Lnodes=if (scenario==1) c("l") else c("l","w"),
          Ynodes="y",
          Qform=qform,
          gform=gform1,
          abar=list(c(1,1),c(0,0)),
          SL.library = SLlib2,
          estimate.time = FALSE)})
  SLtmleimpsumfcs <- lapply(SLtmleimpfcs, function (x) {summary(x)})
  SLtmleimpcofcs <- matrix(unlist(lapply(SLtmleimpsumfcs, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
  SLtmleimpsefcs <- matrix(unlist(lapply(SLtmleimpsumfcs, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
  SLtmlefcs <- matrix(unlist(mi.meld(q=SLtmleimpcofcs, se=SLtmleimpsefcs)),nrow=1,ncol=2)
  
  # MVN Multiple Imputation Analyses
  # Imputed naive analysis
  # GLM
  GLMimpmvn <- lapply(impdatamvn, function (x) {glm(formula=formula(outmodel),data=x,family=gaussian)})
  GLMimpcomvn <- matrix(unlist(lapply(GLMimpmvn, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLMimpsemvn <- matrix(unlist(lapply(GLMimpmvn, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLMmvn <- matrix(unlist(mi.meld(q=GLMimpcomvn, se=GLMimpsemvn)),nrow=1,ncol=2)
  
  # Random intercept
  RIimpmvn <- lapply(impdatamvn, function (x) {lmer(paste0(outmodel,"+(1|id)"),data=x)})
  RIimpcomvn <- matrix(unlist(lapply(RIimpmvn, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  RIimpsemvn <- matrix(unlist(lapply(RIimpmvn, function (x) {sqrt(vcov(x)[2,2] + vcov(x)[3,3] + 2*vcov(x)[2,3])})),nrow=nimpute,ncol=1)
  RImvn <- matrix(unlist(mi.meld(q=RIimpcomvn, se=RIimpsemvn)),nrow=1,ncol=2)
  
  # GEE
  GEEimpmvn <- lapply(impdatamvn, function (x) {geeglm(formula(outmodel),data=x,id=x$id,waves=x$obs,family=gaussian)})
  GEEimpcomvn <- matrix(unlist(lapply(GEEimpmvn, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GEEimpsemvn <- matrix(unlist(lapply(GEEimpmvn, function (x) {sqrt(summary(x)$cov.scaled[2,2] + summary(x)$cov.scaled[3,3] + 2*summary(x)$cov.scaled[2,3])})),nrow=nimpute,ncol=1)
  GEEmvn <- matrix(unlist(mi.meld(q=GEEimpcomvn, se=GEEimpsemvn)),nrow=1,ncol=2)
  
  # GLM-based IPTW and DR-IPTW
  GLiptwimpmvn <- lapply(impdatamvn, function (x) {
    GLexp0 <- glm(formula=propa0model,data=x,family=binomial)
    GLexp1 <- glm(formula=propa1model,data=x,family=binomial)
    GLexpp <- data.table(cbind(id=x$id,obs=x$obs,a=x$a,la=x$la,
                               propa=ifelse(x$a==1,checkrange(predict(GLexp1,type="response")),checkrange(1-predict(GLexp1,type="response"))),
                               propla=ifelse(x$la==1,checkrange(predict(GLexp0,type="response")),checkrange(1-predict(GLexp0,type="response")))))
    GLexpp$p <- GLexpp$propla*GLexpp$propa
    GLexpp$GLwt <- 1/GLexpp$p
    D <- merge(x,GLexpp[,c("id","obs","GLwt")])
    D
  })
  GLiptwimpsummvn <- lapply(GLiptwimpmvn, function (x) {glm(y~a+la,data=x,weight=GLwt,family=gaussian)})
  GLiptwimpcomvn <- matrix(unlist(lapply(GLiptwimpsummvn, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLiptwimpsemvn <- matrix(unlist(lapply(GLiptwimpsummvn, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLiptwmvn <- matrix(unlist(mi.meld(q=GLiptwimpcomvn, se=GLiptwimpsemvn)),nrow=1,ncol=2)
  
  GLdriptwimpsummvn <- lapply(GLiptwimpmvn, function (x) {glm(outmodel2,data=x,weight=GLwt,family=gaussian)})
  GLdriptwimpcomvn <- matrix(unlist(lapply(GLdriptwimpsummvn, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLdriptwimpsemvn <- matrix(unlist(lapply(GLdriptwimpsummvn, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLdriptwmvn <- matrix(unlist(mi.meld(q=GLdriptwimpcomvn, se=GLdriptwimpsemvn)),nrow=1,ncol=2)
  
  # SuperLearner-based IPTW and DR-IPTW
  SLiptwimpmvn <- lapply(impdatamvn, function (x) {
    SLexp0 <- SuperLearner(Y=as.numeric(as.vector(x$la)),X=x[,c("obs","oa","ob","oc","ll")],id=x[,"id"],SL.library=SLlib1,family=binomial)
    SLexp1 <- SuperLearner(Y=as.numeric(as.vector(x$a)),X=x[,c("obs","oa","ob","oc","ll","la","l")],id=x[,"id"],SL.library=SLlib1,family=binomial)
    SLexpp <- data.table(cbind(id=x$id,obs=x$obs,a=x$a,la=x$la,
                               propa=ifelse(x$a==1,checkrange(predict(SLexp1)$pred),checkrange(1-predict(SLexp1)$pred)),
                               propla=ifelse(x$la==1,checkrange(predict(SLexp0)$pred),checkrange(1-predict(SLexp0)$pred))))
    SLexpp$p <- SLexpp$propla*SLexpp$propa
    SLexpp$SLwt <- 1/SLexpp$p
    D <- merge(x,SLexpp[,c("id","obs","SLwt")])
    D
  })
  SLiptwimpsummvn <- lapply(SLiptwimpmvn, function (x) {glm(y~a+la,data=x,weight=SLwt,family=gaussian)})
  SLiptwimpcomvn <- matrix(unlist(lapply(SLiptwimpsummvn, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  SLiptwimpsemvn <- matrix(unlist(lapply(SLiptwimpsummvn, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  SLiptwmvn <- matrix(unlist(mi.meld(q=SLiptwimpcomvn, se=SLiptwimpsemvn)),nrow=1,ncol=2)
  
  SLdriptwimpsummvn <- lapply(SLiptwimpmvn, function (x) {glm(outmodel2,data=x,weight=SLwt,family=gaussian)})
  SLdriptwimpcomvn <- matrix(unlist(lapply(SLdriptwimpsummvn, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  SLdriptwimpsemvn <- matrix(unlist(lapply(SLdriptwimpsummvn, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  SLdriptwmvn <- matrix(unlist(mi.meld(q=SLdriptwimpcomvn, se=SLdriptwimpsemvn)),nrow=1,ncol=2)
  
  # GLM-based TMLE
  GLtmleimpmvn <- lapply(impdatamvn, function (x) {
    ltmle(x[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a","y") else c("obs","oa","ob","oc","ll","lw","la","l","w","a","y")],
          id=x[,"id"],
          Anodes=c("la","a"),
          Lnodes=if (scenario==1) c("l") else c("l","w"),
          Ynodes="y",
          Qform=qform,
          gform=gform1,
          abar=list(c(1,1),c(0,0)),
          estimate.time = FALSE)})
  GLtmleimpsummvn <- lapply(GLtmleimpmvn, function (x) {summary(x)})
  GLtmleimpcomvn <- matrix(unlist(lapply(GLtmleimpsummvn, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
  GLtmleimpsemvn <- matrix(unlist(lapply(GLtmleimpsummvn, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
  GLtmlemvn <- matrix(unlist(mi.meld(q=GLtmleimpcomvn, se=GLtmleimpsemvn)),nrow=1,ncol=2)
  
  # SuperLearner-based TMLE
  SLtmleimpmvn <- lapply(impdatamvn, function (x) {
    ltmle(x[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a","y") else c("obs","oa","ob","oc","ll","lw","la","l","w","a","y")],
          id=x[,"id"],
          Anodes=c("la","a"),
          Lnodes=if (scenario==1) c("l") else c("l","w"),
          Ynodes="y",
          Qform=qform,
          gform=gform1,
          abar=list(c(1,1),c(0,0)),
          SL.library = SLlib2,
          estimate.time = FALSE)})
  SLtmleimpsummvn <- lapply(SLtmleimpmvn, function (x) {summary(x)})
  SLtmleimpcomvn <- matrix(unlist(lapply(SLtmleimpsummvn, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
  SLtmleimpsemvn <- matrix(unlist(lapply(SLtmleimpsummvn, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
  SLtmlemvn <- matrix(unlist(mi.meld(q=SLtmleimpcomvn, se=SLtmleimpsemvn)),nrow=1,ncol=2)
  
  # JOMO Multiple Imputation Analyses
  # Imputed naive analysis
  #GLM
  GLMimpjomo <- lapply(impdatajomo, function (x) {glm(formula=formula(outmodel),data=x,family=gaussian)})
  GLMimpcojomo <- matrix(unlist(lapply(GLMimpjomo, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLMimpsejomo <- matrix(unlist(lapply(GLMimpjomo, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLMjomo <- matrix(unlist(mi.meld(q=GLMimpcojomo, se=GLMimpsejomo)),nrow=1,ncol=2)
  
  # Random intercept
  RIimpjomo <- lapply(impdatajomo, function (x) {lmer(paste0(outmodel,"+(1|id)"),data=x)})
  RIimpcojomo <- matrix(unlist(lapply(RIimpjomo, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  RIimpsejomo <- matrix(unlist(lapply(RIimpjomo, function (x) {sqrt(vcov(x)[2,2] + vcov(x)[3,3] + 2*vcov(x)[2,3])})),nrow=nimpute,ncol=1)
  RIjomo <- matrix(unlist(mi.meld(q=RIimpcojomo, se=RIimpsejomo)),nrow=1,ncol=2)
  
  # GEE
  GEEimpjomo <- lapply(impdatajomo, function (x) {geeglm(formula(outmodel),data=x,id=x$id,waves=x$obs,family=gaussian)})
  GEEimpcojomo <- matrix(unlist(lapply(GEEimpjomo, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GEEimpsejomo <- matrix(unlist(lapply(GEEimpjomo, function (x) {sqrt(summary(x)$cov.scaled[2,2] + summary(x)$cov.scaled[3,3] + 2*summary(x)$cov.scaled[2,3])})),nrow=nimpute,ncol=1)
  GEEjomo <- matrix(unlist(mi.meld(q=GEEimpcojomo, se=GEEimpsejomo)),nrow=1,ncol=2)
  
  # GLM-based IPTW and DR-IPTW
  GLiptwimpjomo <- lapply(impdatajomo, function (x) {
    GLexp0 <- glm(formula=propa0model,data=x,family=binomial)
    GLexp1 <- glm(formula=propa1model,data=x,family=binomial)
    GLexpp <- data.table(cbind(id=x$id,obs=x$obs,a=x$a,la=x$la,
                               propa=ifelse(x$a==1,checkrange(predict(GLexp1,type="response")),checkrange(1-predict(GLexp1,type="response"))),
                               propla=ifelse(x$la==1,checkrange(predict(GLexp0,type="response")),checkrange(1-predict(GLexp0,type="response")))))
    GLexpp$p <- GLexpp$propla*GLexpp$propa
    GLexpp$GLwt <- 1/GLexpp$p
    D <- merge(x,GLexpp[,c("id","obs","GLwt")])
    D
  })
  GLiptwimpsumjomo <- lapply(GLiptwimpjomo, function (x) {glm(y~a+la,data=x,weight=GLwt,family=gaussian)})
  GLiptwimpcojomo <- matrix(unlist(lapply(GLiptwimpsumjomo, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLiptwimpsejomo <- matrix(unlist(lapply(GLiptwimpsumjomo, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLiptwjomo <- matrix(unlist(mi.meld(q=GLiptwimpcojomo, se=GLiptwimpsejomo)),nrow=1,ncol=2)
  
  GLdriptwimpsumjomo <- lapply(GLiptwimpjomo, function (x) {glm(outmodel2,data=x,weight=GLwt,family=gaussian)})
  GLdriptwimpcojomo <- matrix(unlist(lapply(GLdriptwimpsumjomo, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  GLdriptwimpsejomo <- matrix(unlist(lapply(GLdriptwimpsumjomo, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  GLdriptwjomo <- matrix(unlist(mi.meld(q=GLdriptwimpcojomo, se=GLdriptwimpsejomo)),nrow=1,ncol=2)
  
  # SuperLearner-based IPTW and DR-IPTW
  SLiptwimpjomo <- lapply(impdatajomo, function (x) {
    SLexp0 <- SuperLearner(Y=as.numeric(levels(x$la))[x$la],X=x[,c("obs","oa","ob","oc","ll")],id=x[,"id"],SL.library=SLlib1,family=binomial)
    SLexp1 <- SuperLearner(Y=as.numeric(levels(x$a))[x$a],X=x[,c("obs","oa","ob","oc","ll","la","l")],id=x[,"id"],SL.library=SLlib1,family=binomial)
    SLexpp <- data.table(cbind(id=x$id,obs=x$obs,a=x$a,la=x$la,
                               propa=ifelse(x$a==1,checkrange(predict(SLexp1)$pred),checkrange(1-predict(SLexp1)$pred)),
                               propla=ifelse(x$la==1,checkrange(predict(SLexp0)$pred),checkrange(1-predict(SLexp0)$pred))))
    SLexpp$p <- SLexpp$propla*SLexpp$propa
    SLexpp$SLwt <- 1/SLexpp$p
    D <- merge(x,SLexpp[,c("id","obs","SLwt")])
    D
  })
  SLiptwimpsumjomo <- lapply(SLiptwimpjomo, function (x) {glm(y~a+la,data=x,weight=SLwt,family=gaussian)})
  SLiptwimpcojomo <- matrix(unlist(lapply(SLiptwimpsumjomo, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  SLiptwimpsejomo <- matrix(unlist(lapply(SLiptwimpsumjomo, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  SLiptwjomo <- matrix(unlist(mi.meld(q=SLiptwimpcojomo, se=SLiptwimpsejomo)),nrow=1,ncol=2)
  
  SLdriptwimpsumjomo <- lapply(SLiptwimpjomo, function (x) {glm(outmodel2,data=x,weight=SLwt,family=gaussian)})
  SLdriptwimpcojomo <- matrix(unlist(lapply(SLdriptwimpsumjomo, function (x) {coef(summary(x))[2,1]+coef(summary(x))[3,1]})),nrow=nimpute,ncol=1)
  SLdriptwimpsejomo <- matrix(unlist(lapply(SLdriptwimpsumjomo, function (x) {sqrt(vcovCR(x, cluster=data$id, type = "CR3")[2,2] + vcovCR(x, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(x, cluster=data$id, type = "CR3")[2,3])})),nrow=nimpute,ncol=1)
  SLdriptwjomo <- matrix(unlist(mi.meld(q=SLdriptwimpcojomo, se=SLdriptwimpsejomo)),nrow=1,ncol=2)
  
  # GLM-based TMLE
  GLtmleimpjomo <- lapply(impdatajomo, function (x) {
    x$ll <- as.numeric(levels(x$ll))[x$ll]
    x$la <- as.numeric(levels(x$la))[x$la]
    x$l <- as.numeric(levels(x$l))[x$l]
    x$a <- as.numeric(levels(x$a))[x$a]
    ltmle(x[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a","y") else c("obs","oa","ob","oc","ll","lw","la","l","w","a","y")],
          id=x[,"id"],
          Anodes=c("la","a"),
          Lnodes=if (scenario==1) c("l") else c("l","w"),
          Ynodes="y",
          Qform=qform,
          gform=gform1,
          abar=list(c(1,1),c(0,0)),
          estimate.time = FALSE)})
  GLtmleimpsumjomo <- lapply(GLtmleimpjomo, function (x) {summary(x)})
  GLtmleimpcojomo <- matrix(unlist(lapply(GLtmleimpsumjomo, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
  GLtmleimpsejomo <- matrix(unlist(lapply(GLtmleimpsumjomo, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
  GLtmlejomo <- matrix(unlist(mi.meld(q=GLtmleimpcojomo, se=GLtmleimpsejomo)),nrow=1,ncol=2)
  
  # SuperLearner-based TMLE
  SLtmleimpjomo <- lapply(impdatajomo, function (x) {
    x$ll <- as.numeric(levels(x$ll))[x$ll]
    x$la <- as.numeric(levels(x$la))[x$la]
    x$l <- as.numeric(levels(x$l))[x$l]
    x$a <- as.numeric(levels(x$a))[x$a]
    ltmle(x[,if (scenario==1) c("obs","oa","ob","oc","ll","la","l","a","y") else c("obs","oa","ob","oc","ll","lw","la","l","w","a","y")],
          id=x[,"id"],
          Anodes=c("la","a"),
          Lnodes=if (scenario==1) c("l") else c("l","w"),
          Ynodes="y",
          Qform=qform,
          gform=gform1,
          abar=list(c(1,1),c(0,0)),
          SL.library = SLlib2,
          estimate.time = FALSE)})
  SLtmleimpsumjomo <- lapply(SLtmleimpjomo, function (x) {summary(x)})
  SLtmleimpcojomo <- matrix(unlist(lapply(SLtmleimpsumjomo, function (x) {x$effect.measures$ATE$estimate})),nrow=nimpute,ncol=1)
  SLtmleimpsejomo <- matrix(unlist(lapply(SLtmleimpsumjomo, function (x) {x$effect.measures$ATE$std.dev})),nrow=nimpute,ncol=1)
  SLtmlejomo <- matrix(unlist(mi.meld(q=SLtmleimpcojomo, se=SLtmleimpsejomo)),nrow=1,ncol=2)

  # Combine results into vector
  c(coef(summary(GLGLMc))[2,1]+coef(summary(GLGLMc))[3,1],sqrt(vcovCR(GLGLMc, cluster=data$id, type = "CR3")[2,2] + vcovCR(GLGLMc, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(GLGLMc, cluster=data$id, type = "CR3")[2,3]),
    coef(summary(GLRIc))[2,1]+coef(summary(GLRIc))[3,1],sqrt(vcov(GLRIc)[2,2] + vcov(GLRIc)[3,3] + 2*vcov(GLRIc)[2,3]),
    coef(summary(GLGEEc))[2,1]+coef(summary(GLGEEc))[3,1],sqrt(summary(GLGEEc)$cov.scaled[2,2] + summary(GLGEEc)$cov.scaled[3,3] + 2*summary(GLGEEc)$cov.scaled[2,3]),
    coef(summary(GLiptcw))[2,1]+coef(summary(GLiptcw))[3,1],sqrt(vcovCR(GLiptcw, cluster=data$id, type = "CR3")[2,2] + vcovCR(GLiptcw, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(GLiptcw, cluster=data$id, type = "CR3")[2,3]),
    coef(summary(GLdriptcw))[2,1]+coef(summary(GLdriptcw))[3,1],sqrt(vcovCR(GLdriptcw, cluster=data$id, type = "CR3")[2,2] + vcovCR(GLdriptcw, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(GLdriptcw, cluster=data$id, type = "CR3")[2,3]),
    summary(GLtmlec)$effect.measures$ATE$estimate,summary(GLtmlec)$effect.measures$ATE$std.dev,
    coef(summary(SLGLMc))[2,1]+coef(summary(SLGLMc))[3,1],sqrt(vcovCR(SLGLMc, cluster=data$id, type = "CR3")[2,2] + vcovCR(SLGLMc, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(SLGLMc, cluster=data$id, type = "CR3")[2,3]),
    coef(summary(SLRIc))[2,1]+coef(summary(SLRIc))[3,1],sqrt(vcov(SLRIc)[2,2] + vcov(SLRIc)[3,3] + 2*vcov(SLRIc)[2,3]),
    coef(summary(SLGEEc))[2,1]+coef(summary(SLGEEc))[3,1],sqrt(summary(SLGEEc)$cov.scaled[2,2] + summary(SLGEEc)$cov.scaled[3,3] + 2*summary(SLGEEc)$cov.scaled[2,3]),
    coef(summary(SLiptcw))[2,1]+coef(summary(SLiptcw))[3,1],sqrt(vcovCR(SLiptcw, cluster=data$id, type = "CR3")[2,2] + vcovCR(SLiptcw, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(SLiptcw, cluster=data$id, type = "CR3")[2,3]),
    coef(summary(SLdriptcw))[2,1]+coef(summary(SLdriptcw))[3,1],sqrt(vcovCR(SLdriptcw, cluster=data$id, type = "CR3")[2,2] + vcovCR(SLdriptcw, cluster=data$id, type = "CR3")[3,3] + 2*vcovCR(SLdriptcw, cluster=data$id, type = "CR3")[2,3]),
    summary(SLtmlec)$effect.measures$ATE$estimate,summary(SLtmlec)$effect.measures$ATE$std.dev,
    GLMfcs[1],GLMfcs[2],
    RIfcs[1],RIfcs[2],
    GEEfcs[1],GEEfcs[2],
    GLiptwfcs[1],GLiptwfcs[2],
    GLdriptwfcs[1],GLdriptwfcs[2],
    GLtmlefcs[1],GLtmlefcs[2],
    SLiptwfcs[1],SLiptwfcs[2],
    SLdriptwfcs[1],SLdriptwfcs[2],
    SLtmlefcs[1],SLtmlefcs[2],
    GLMmvn[1],GLMmvn[2],
    RImvn[1],RImvn[2],
    GEEmvn[1],GEEmvn[2],
    GLiptwmvn[1],GLiptwmvn[2],
    GLdriptwmvn[1],GLdriptwmvn[2],
    GLtmlemvn[1],GLtmlemvn[2],
    SLiptwmvn[1],SLiptwmvn[2],
    SLdriptwmvn[1],SLdriptwmvn[2],
    SLtmlemvn[1],SLtmlemvn[2],
    GLMjomo[1],GLMjomo[2],
    RIjomo[1],RIjomo[2],
    GEEjomo[1],GEEjomo[2],
    GLiptwjomo[1],GLiptwjomo[2],
    GLdriptwjomo[1],GLdriptwjomo[2],
    GLtmlejomo[1],GLtmlejomo[2],
    SLiptwjomo[1],SLiptwjomo[2],
    SLdriptwjomo[1],SLdriptwjomo[2],
    SLtmlejomo[1],SLtmlejomo[2])
}

col.names1 <- c("GLMjtco","GLMjtse","RIjtco","RIjtse","GEEjtco","GEEjtse",
               "GLMiptwjtco","GLMiptwjtse","SLiptwjtco","SLiptwjtse",
               "GLMdriptwjtco","GLMdriptwjtse","SLdriptwjtco","SLdriptwjtse",
               "GLMtmlejtco","GLMtmlejtse","SLtmlejtco","SLtmlejtse")

col.names2 <- c("GLMjtco","GLMjtse","RIjtco","RIjtse","GEEjtco","GEEjtse",
               "GLMiptwjtco","GLMiptwjtse","SLiptwjtco","SLiptwjtse",
               "GLMdriptwjtco","GLMdriptwjtse","SLdriptwjtco","SLdriptwjtse",
               "GLMtmlejtco","GLMtmlejtse","SLtmlejtco","SLtmlejtse",
               "GLMIPC_GLMjtco","GLMIPC_GLMjtse","GLMIPC_RIjtco","GLMIPC_RIjtse","GLMIPC_GEEjtco","GLMIPC_GEEjtse",
               "GLMIPC_GLMiptwjtco","GLMIPC_GLMiptwjtse","GLMIPC_GLMdriptwjtco","GLMIPC_GLMdriptwjtse","GLMIPC_GLMtmlejtco","GLMIPC_GLMtmlejtse",
               "SLIPC_GLMjtco","SLIPC_GLMjtse","SLIPC_RIjtco","SLIPC_RIjtse","SLIPC_GEEjtco","SLIPC_GEEjtse",
               "SLIPC_SLiptwjtco","SLIPC_SLiptwjtse","SLIPC_SLdriptwjtco","SLIPC_SLdriptwjtse","SLIPC_SLtmlejtco","SLIPC_SLtmlejtse",
               "FCS_GLMjtco","FCS_GLMjtse","FCS_RIjtco","FCS_RIjtse","FCS_GEEjtco","FCS_GEEjtse",
               "FCS_GLMiptwjtco","FCS_GLMiptwjtse","FCS_SLiptwjtco","FCS_SLiptwjtse",
               "FCS_GLMdriptwjtco","FCS_GLMdriptwjtse","FCS_SLdriptwjtco","FCS_SLdriptwjtse",
               "FCS_GLMtmlejtco","FCS_GLMtmlejtse","FCS_SLtmlejtco","FCS_SLtmlejtse",
               "MVN_GLMjtco","MVN_GLMjtse","MVN_RIjtco","MVN_RIjtse","MVN_GEEjtco","MVN_GEEjtse",
               "MVN_GLMiptwjtco","MVN_GLMiptwjtse","MVN_SLiptwjtco","MVN_SLiptwjtse",
               "MVN_GLMdriptwjtco","MVN_GLMdriptwjtse","MVN_SLdriptwjtco","MVN_SLdriptwjtse",
               "MVN_GLMtmlejtco","MVN_GLMtmlejtse","MVN_SLtmlejtco","MVN_SLtmlejtse",
               "JOMO_GLMjtco","JOMO_GLMjtse","JOMO_RIjtco","JOMO_RIjtse","JOMO_GEEjtco","JOMO_GEEjtse",
               "JOMO_GLMiptwjtco","JOMO_GLMiptwjtse","JOMO_SLiptwjtco","JOMO_SLiptwjtse",
               "JOMO_GLMdriptwjtco","JOMO_GLMdriptwjtse","JOMO_SLdriptwjtco","JOMO_SLdriptwjtse",
               "JOMO_GLMtmlejtco","JOMO_GLMtmlejtse","JOMO_SLtmlejtco","JOMO_SLtmlejtse")

# ############################
# ##       Scenario 1       ##
# ############################
# 
# # Complete data - no missing data
# datacomp1 <- lapply(paste0(filepath1,list.files(path=filepath1))[start:stop],function (x) {
#   data <- data.frame(read_dta(x))
#   data$id <- factor(data$id)
#   data <- data.frame(data[,c(-11,-12,-13,-14,-15,-16)])
#   data
# })
# rescomp1 <- matrix(unlist(lapply(datacomp1, function (x) {simcomp(x)})),nrow=n,ncol=18,byrow=TRUE)
# colnames(rescomp1) <- col.names1
# write_dta(data.frame(rescomp1),path=paste0(paste0(paste0(args[3],"Results/S1comp-"),start),".dta"))
# rm(datacomp1)
# 
# # 25% Missing data scenario
# # MCAR - only complete case analysis run (because it is unbiased)
# datamcara1 <- lapply(paste0(filepath1,list.files(path=filepath1))[start:stop],function (x) {
#   data <- data.frame(read_dta(x))
#   data$id <- factor(data$id)
#   data$y <- ifelse(data$missmcarmod==1,NA,data$y)
#   data <- data.frame(cbind(data[,1:9],c=data[,"missmcarmod"],y=data[,"y"]))
#   data
# })
# resmcara1 <- matrix(unlist(lapply(datamcara1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
# colnames(resmcara1) <- col.names1
# write_dta(data.frame(resmcara1),path=paste0(paste0(paste0(args[3],"Results/S1mod-mcar-"),start),".dta"))
# rm(datamcara1)
# 
# # MAR A - Missingness only related to L
# datamara1 <- lapply(paste0(filepath1,list.files(path=filepath1))[start:stop],function (x) {
#   data <- data.frame(read_dta(x))
#   data$id <- factor(data$id)
#   data$y <- ifelse(data$missmarmod1==1,NA,data$y)
#   data <- data.frame(cbind(data[,1:9],c=data[,"missmarmod1"],y=data[,"y"]))
#   data
# })
# resmara1a <- matrix(unlist(lapply(datamara1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
# resmara1b <- matrix(unlist(lapply(datamara1, function (x) {simmiss(x,nimpute=nimp,scenario=1)})),nrow=n,ncol=78,byrow=TRUE)
# resmara1 <- cbind(resmara1a,resmara1b)
# colnames(resmara1) <- col.names2
# write_dta(data.frame(resmara1),path=paste0(paste0(paste0(args[3],"Results/S1mod-mar1-"),start),".dta"))
# rm(datamara1)
# 
# # MAR B - Missingness related to L and A
# datamarb1 <- lapply(paste0(filepath1,list.files(path=filepath1))[start:stop],function (x) {
#   data <- data.frame(read_dta(x))
#   data$id <- factor(data$id)
#   data$y <- ifelse(data$missmarmod2==1,NA,data$y)
#   data <- data.frame(cbind(data[,1:9],c=data[,"missmarmod2"],y=data[,"y"]))
#   data
# })
# resmarb1a <- matrix(unlist(lapply(datamarb1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
# resmarb1b <- matrix(unlist(lapply(datamarb1, function (x) {simmiss(x,nimpute=nimp,scenario=1)})),nrow=n,ncol=78,byrow=TRUE)
# resmarb1 <- cbind(resmarb1a,resmarb1b)
# colnames(resmarb1) <- col.names2
# write_dta(data.frame(resmarb1),path=paste0(paste0(paste0(args[3],"Results/S1mod-mar2-"),start),".dta"))
# rm(datamarb1)
# 
# # 50% Missing data scenario
# # MCAR - only complete case analysis run (because it is unbiased)
# datamcara2 <- lapply(paste0(filepath1,list.files(path=filepath1))[start:stop],function (x) {
#   data <- data.frame(read_dta(x))
#   data$id <- factor(data$id)
#   data$y <- ifelse(data$missmcarsev==1,NA,data$y)
#   data <- data.frame(cbind(data[,1:9],c=data[,"missmcarsev"],y=data[,"y"]))
#   data
# })
# resmcara2 <- matrix(unlist(lapply(datamcara2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
# colnames(resmcara2) <- col.names1
# write_dta(data.frame(resmcara2),path=paste0(paste0(paste0(args[3],"Results/S1sev-mcar-"),start),".dta"))
# rm(datamcara2)
# 
# # MAR A - Missingness only related to L
# datamara2 <- lapply(paste0(filepath1,list.files(path=filepath1))[start:stop],function (x) {
#   data <- data.frame(read_dta(x))
#   data$id <- factor(data$id)
#   data$y <- ifelse(data$missmarsev1==1,NA,data$y)
#   data <- data.frame(cbind(data[,1:9],c=data[,"missmarsev1"],y=data[,"y"]))
#   data
# })
# resmara2a <- matrix(unlist(lapply(datamara2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
# resmara2b <- matrix(unlist(lapply(datamara2, function (x) {simmiss(x,nimpute=nimp,scenario=1)})),nrow=n,ncol=78,byrow=TRUE)
# resmara2 <- cbind(resmara2a,resmara2b)
# colnames(resmara2) <- col.names2
# write_dta(data.frame(resmara2),path=paste0(paste0(paste0(args[3],"Results/S1sev-mar1-"),start),".dta"))
# rm(datamara2)
# 
# # MAR B - Missingness related to L and A
# datamarb2 <- lapply(paste0(filepath1,list.files(path=filepath1))[start:stop],function (x) {
#   data <- data.frame(read_dta(x))
#   data$id <- factor(data$id)
#   data$y <- ifelse(data$missmarsev2==1,NA,data$y)
#   data <- data.frame(cbind(data[,1:9],c=data[,"missmarsev2"],y=data[,"y"]))
#   data
# })
# resmarb2a <- matrix(unlist(lapply(datamarb2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
# resmarb2b <- matrix(unlist(lapply(datamarb2, function (x) {simmiss(x,nimpute=nimp,scenario=1)})),nrow=n,ncol=78,byrow=TRUE)
# resmarb2 <- cbind(resmarb2a,resmarb2b)
# colnames(resmarb2) <- col.names2
# write_dta(data.frame(resmarb2),path=paste0(paste0(paste0(args[3],"Results/S1sev-mar2-"),start),".dta"))
# rm(datamarb2)

############################
##       Scenario 2       ##
############################

# Complete data - no missing data
datacomp2 <- lapply(paste0(filepath2,list.files(path=filepath2))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data <- data.frame(data[,c(-13,-14,-15,-16,-17,-18)])
  data
})
rescomp2 <- matrix(unlist(lapply(datacomp2, function (x) {simcomp(x)})),nrow=n,ncol=18,byrow=TRUE)
colnames(rescomp2) <- col.names1
write_dta(data.frame(rescomp2),path=paste0(paste0(paste0(args[3],"Results/S2comp-"),start),".dta"))
rm(datacomp2)

# 25% Missing data scenario
# MCAR - only complete case analysis run (because it is unbiased)
datamcarb1 <- lapply(paste0(filepath2,list.files(path=filepath2))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmcarmod==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmcarmod"],y=data[,"y"]))
  data
})
resmcarb1 <- matrix(unlist(lapply(datamcarb1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
colnames(resmcarb1) <- col.names1
write_dta(data.frame(resmcarb1),path=paste0(paste0(paste0(args[3],"Results/S2mod-mcar-"),start),".dta"))
rm(datamcarb1)

# MAR C - Missingness only related to W
datamarc1 <- lapply(paste0(filepath2,list.files(path=filepath2))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmarmod1==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmarmod1"],y=data[,"y"]))
  data
})
resmarc1a <- matrix(unlist(lapply(datamarc1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
resmarc1b <- matrix(unlist(lapply(datamarc1, function (x) {simmiss(x,nimpute=nimp,scenario=2)})),nrow=n,ncol=78,byrow=TRUE)
resmarc1 <- cbind(resmarc1a,resmarc1b)
colnames(resmarc1) <- col.names2
write_dta(data.frame(resmarc1),path=paste0(paste0(paste0(args[3],"Results/S2mod-mar1-"),start),".dta"))
rm(datamarc1)

# MAR D - Missingness related to W and A
datamard1 <- lapply(paste0(filepath2,list.files(path=filepath2))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmarmod2==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmarmod2"],y=data[,"y"]))
  data
})
resmard1a <- matrix(unlist(lapply(datamard1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
resmard1b <- matrix(unlist(lapply(datamard1, function (x) {simmiss(x,nimpute=nimp,scenario=2)})),nrow=n,ncol=78,byrow=TRUE)
resmard1 <- cbind(resmard1a,resmard1b)
colnames(resmard1) <- col.names2
write_dta(data.frame(resmard1),path=paste0(paste0(paste0(args[3],"Results/S2mod-mar2-"),start),".dta"))
rm(datamard1)

# 50% Missing data scenario
# MCAR - only complete case analysis run (because it is unbiased)
datamcarb2 <- lapply(paste0(filepath2,list.files(path=filepath2))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmcarsev==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmcarsev"],y=data[,"y"]))
  data
})
resmcarb2 <- matrix(unlist(lapply(datamcarb2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
colnames(resmcarb2) <- col.names1
write_dta(data.frame(resmcarb2),path=paste0(paste0(paste0(args[3],"Results/S2sev-mcar-"),start),".dta"))
rm(datamcarb2)

# MAR C - Missingness only related to W
datamarc2 <- lapply(paste0(filepath2,list.files(path=filepath2))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmarsev1==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmarsev1"],y=data[,"y"]))
  data
})
resmarc2a <- matrix(unlist(lapply(datamarc2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
resmarc2b <- matrix(unlist(lapply(datamarc2, function (x) {simmiss(x,nimpute=nimp,scenario=2)})),nrow=n,ncol=78,byrow=TRUE)
resmarc2 <- cbind(resmarc2a,resmarc2b)
colnames(resmarc2) <- col.names2
write_dta(data.frame(resmarc2),path=paste0(paste0(paste0(args[3],"Results/S2sev-mar1-"),start),".dta"))
rm(datamarc2)

# MAR D - Missingness related to W and A
datamard2 <- lapply(paste0(filepath2,list.files(path=filepath2))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmarsev2==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmarsev2"],y=data[,"y"]))
  data
})
resmard2a <- matrix(unlist(lapply(datamard2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
resmard2b <- matrix(unlist(lapply(datamard2, function (x) {simmiss(x,nimpute=nimp,scenario=2)})),nrow=n,ncol=78,byrow=TRUE)
resmard2 <- cbind(resmard2a,resmard2b)
colnames(resmard2) <- col.names2
write_dta(data.frame(resmard2),path=paste0(paste0(paste0(args[3],"Results/S2sev-mar2-"),start),".dta"))
rm(datamard2)

############################
##       Scenario 3       ##
############################

# Complete data - no missing data
datacomp3 <- lapply(paste0(filepath3,list.files(path=filepath3))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data <- data.frame(data[,c(-13,-14,-15,-16,-17,-18)])
  data
})
rescomp3 <- matrix(unlist(lapply(datacomp3, function (x) {simcomp(x)})),nrow=n,ncol=18,byrow=TRUE)
colnames(rescomp3) <- col.names1
write_dta(data.frame(rescomp3),path=paste0(paste0(paste0(args[3],"Results/S3comp-"),start),".dta"))
rm(datacomp3)

# 25% Missing data scenario
# MCAR - only complete case analysis run (because it is unbiased)
datamcarc1 <- lapply(paste0(filepath3,list.files(path=filepath3))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmcarmod==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmcarmod"],y=data[,"y"]))
  data
})
resmcarc1 <- matrix(unlist(lapply(datamcarc1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
colnames(resmcarc1) <- col.names1
write_dta(data.frame(resmcarc1),path=paste0(paste0(paste0(args[3],"Results/S3mod-mcar-"),start),".dta"))
rm(datamcarc1)

# MAR E - Missingness only related to W
datamare1 <- lapply(paste0(filepath3,list.files(path=filepath3))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmarmod1==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmarmod1"],y=data[,"y"]))
  data
})
resmare1a <- matrix(unlist(lapply(datamare1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
resmare1b <- matrix(unlist(lapply(datamare1, function (x) {simmiss(x,nimpute=nimp,scenario=2)})),nrow=n,ncol=78,byrow=TRUE)
resmare1 <- cbind(resmare1a,resmare1b)
colnames(resmare1) <- col.names2
write_dta(data.frame(resmare1),path=paste0(paste0(paste0(args[3],"Results/S3mod-mar1-"),start),".dta"))
rm(datamare1)

# MAR F - Missingness related to W and A
datamarf1 <- lapply(paste0(filepath3,list.files(path=filepath3))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmarmod2==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmarmod2"],y=data[,"y"]))
  data
})
resmarf1a <- matrix(unlist(lapply(datamarf1, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
resmarf1b <- matrix(unlist(lapply(datamarf1, function (x) {simmiss(x,nimpute=nimp,scenario=2)})),nrow=n,ncol=78,byrow=TRUE)
resmarf1 <- cbind(resmarf1a,resmarf1b)
colnames(resmarf1) <- col.names2
write_dta(data.frame(resmarf1),path=paste0(paste0(paste0(args[3],"Results/S3mod-mar2-"),start),".dta"))
rm(datamarf1)

# 50% Missing data scenario
# MCAR - only complete case analysis run (because it is unbiased)
datamcarc2 <- lapply(paste0(filepath3,list.files(path=filepath3))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmcarsev==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmcarsev"],y=data[,"y"]))
  data
})
resmcarc2 <- matrix(unlist(lapply(datamcarc2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
colnames(resmcarc2) <- col.names1
write_dta(data.frame(resmcarc2),path=paste0(paste0(paste0(args[3],"Results/S3sev-mcar-"),start),".dta"))
rm(datamcarc2)

# MAR E - Missingness only related to W
datamare2 <- lapply(paste0(filepath3,list.files(path=filepath3))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmarsev1==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmarsev1"],y=data[,"y"]))
  data
})
resmare2a <- matrix(unlist(lapply(datamare2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
resmare2b <- matrix(unlist(lapply(datamare2, function (x) {simmiss(x,nimpute=nimp,scenario=2)})),nrow=n,ncol=78,byrow=TRUE)
resmare2 <- cbind(resmare2a,resmare2b)
colnames(resmare2) <- col.names2
write_dta(data.frame(resmare2),path=paste0(paste0(paste0(args[3],"Results/S3sev-mar1-"),start),".dta"))
rm(datamare2)

# MAR F - Missingness related to W and A
datamarf2 <- lapply(paste0(filepath3,list.files(path=filepath3))[start:stop],function (x) {
  data <- data.frame(read_dta(x))
  data$id <- factor(data$id)
  data$y <- ifelse(data$missmarsev2==1,NA,data$y)
  data <- data.frame(cbind(data[,1:11],c=data[,"missmarsev2"],y=data[,"y"]))
  data
})
resmarf2a <- matrix(unlist(lapply(datamarf2, function (x) {simcomp(x[,!(names(x) %in% "c")])})),nrow=n,ncol=18,byrow=TRUE)
resmarf2b <- matrix(unlist(lapply(datamarf2, function (x) {simmiss(x,nimpute=nimp,scenario=2)})),nrow=n,ncol=78,byrow=TRUE)
resmarf2 <- cbind(resmarf2a,resmarf2b)
colnames(resmarf2) <- col.names2
write_dta(data.frame(resmarf2),path=paste0(paste0(paste0(args[3],"Results/S3sev-mar2-"),start),".dta"))
rm(datamarf2)