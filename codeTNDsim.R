#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------
# Simulations DAGs for test-negative design
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

rm(list=ls())
set.seed(1)

#number of digits to round to in means, sd's etc
ndig <- 3

#size population simulated
N <- 10^5
#number of times to simulate population
nrep <- 500
#number of cases and controls in case-control and cohort design (inside the minMod function)
K <- 5000

#grid for coefficients
incr <- log(c(0.1, 0.5, 1, 1.5, 2, 3)) #for c2
incr01 <- log(c(0.01, 0.2, 0.4, 0.6, 0.8, 1)) #for c3, c4
ngrid <- length(incr)
#------------------

#make log odds ratio out of two probabilities 
# input: p=P[Y=1|X=1], q=P[Y=1|X=0]
# output: log(OR_{XY})
funlogOR <- function(p,q){log((p/(1-p))/(q/(1-q)))}

#make logit function
funlogit <- function(p){log(p/(1-p))}

#make expit function
funexpit <- function(x){1/(1+exp(-x))}

#-----------------------------------------------------------------------------------
# Fixed parameters
#-----------------------------------------------------------------------------------
#C
#P[C=1]
pC <- 0.5 

#---------------------
#H
#c(P[H=1|C=0], P(H=1|C=1))
aHaux <- c(0.3,0.7) 
#c(logit(P[H=1|C=0]), log(OR_{HC})) = c(a1,a2)
aH <- c(funlogit(aHaux[1]), funlogOR(aHaux[2],aHaux[1]))
#compute marginal: P[H=1]
pH <- funexpit(aH[1]+aH[2])*pC + funexpit(aH[1])*(1-pC)

#---------------------
#V
#c(P[V=1|C=0,H=0], P[V=1|C=1,H=0], P[V=1|C=0,H=1])
bVaux <- c(0.1,0.2,0.3)
#c(logit(P[V=1|C=0,V=0]), log(OR_{CV}(H=0)), log(OR_{HV}(C=0))) = c(b1,b2,b3)
bV <- c(funlogit(bVaux[1]),funlogOR(bVaux[2],bVaux[1]),funlogOR(bVaux[3],bVaux[1]))
#compute marginal: P[V=1]: from parameters
pV <- funexpit(bV[1])*(1-funexpit(aH[1]))*(1-pC) + #(C=0,H=0)
  funexpit(bV[1]+bV[2])*(1-funexpit(aH[1]+aH[2]))*pC + #(C=1,H=0)
  funexpit(bV[1]+bV[3])*funexpit(aH[1])*(1-pC) + #(C=0,H=1)
  funexpit(bV[1]+bV[2]+bV[3])*funexpit(aH[1]+aH[2])*pC #(C=1,H=1)
#Sanity check: compute P[V=1|C=1,H=1] #note: this is fixed (does not change with grid)
pVc1h1 <- funexpit(bV[1]+bV[2]+bV[3])

#---------------------
#I
#c(P[I=1|C=0,H=0,V=0], P[I=1|C=1,H=0,V=0], P[I=1|C=0,H=1,V=0], P[I=1|C=0,H=0,V=1])
cIaux <- c(0.3,0.9,0.5,0.2) #just test-values, because you generate grid for these coefficients!!
#c(logit(P[I=1|C=0,H=0,V=0]), log(OR_{CI}(H=0,V=0)), log(OR_{HI}(C=0,V=0)), log(OR_{VI}(C=0,H=0))) = c(c1,c2,c3,c4)
cI <- c(funlogit(cIaux[1]), funlogOR(cIaux[2],cIaux[1]), funlogOR(cIaux[3],cIaux[1]), funlogOR(cIaux[4],cIaux[1]))
#compute P[I=1|V=0] from parameters
pInovacc <- (1/(1-pV)) * ( funexpit(cI[1]) * (1-funexpit(aH[1])) * (1-pC) * (1-funexpit(bV[1])) + #(C=0,H=0)
                             funexpit(cI[1]+cI[2]) * (1-funexpit(aH[1]+aH[2])) *pC * (1-funexpit(bV[1]+bV[2])) + #(C=1,H=0)
                             funexpit(cI[1]+cI[3]) * funexpit(aH[1]) * (1-pC) * (1-funexpit(bV[1]+bV[3])) + #(C=0,H=1)
                             funexpit(cI[1]+cI[2]+cI[3]) * funexpit(aH[1]+aH[2]) * pC * (1-funexpit(bV[1]+bV[2]+bV[3])) ) #(C=1,H=1)
#compute P[I=1|V=1] from parameters
pIvacc <- (1/pV) * ( funexpit(cI[1]+cI[4]) * (1-funexpit(aH[1])) * (1-pC) * funexpit(bV[1]) +
                       funexpit(cI[1]+cI[2]+cI[4]) * (1-funexpit(aH[1]+aH[2])) * pC * funexpit(bV[1]+bV[2]) +
                       funexpit(cI[1]+cI[3]+cI[4]) * funexpit(aH[1]) * (1-pC) * funexpit(bV[1]+bV[3]) +
                       funexpit(cI[1]+cI[2]+cI[3]+cI[4]) * funexpit(aH[1]+aH[2]) * pC * funexpit(bV[1]+bV[2]+bV[3]) )
#compute marginal: P[I=1]
pI <- pIvacc * pV + pInovacc * (1-pV)

#fix a value for P[I=1|V=0] such that, together with given grid for c2, c3, c4, can compute c1
pInovaccSET <- 0.2

#---------------------
#S
#c(P[S=1|H=0,I=0],P[S=1|H=1,I=0],P[S=1|H=0,I=1])
dSaux <- c(0.01,0.1,0.4)
#c(logit(P[S=1|H=0,I=0]), log(OR_{HS}(I=0)), log(OR_{IS}(H=0))) = c(d1,d2,d3)
dS <- c(funlogit(dSaux[1]), funlogOR(dSaux[2],dSaux[1]), funlogOR(dSaux[3],dSaux[1]))
#compute some auxiliary quantities:
#compute P[I=0|H=0]
pi0h0 <- (1/(1-pH)) * ( (1-funexpit(cI[1])) * (1-funexpit(bV[1])) * (1-funexpit(aH[1])) * (1-pC) +
                          (1-funexpit(cI[1]+cI[4])) * funexpit(bV[1]) * (1-funexpit(aH[1])) * (1-pC) +
                          (1-funexpit(cI[1]+cI[2])) * (1-funexpit(bV[1]+bV[2])) * (1-funexpit(aH[1]+aH[2])) * pC +
                          (1-funexpit(cI[1]+cI[2]+cI[4])) * funexpit(bV[1]+bV[2]) * (1-funexpit(aH[1]+aH[2])) * pC ) 
#compute P[I=0|H=1]
pi0h1 <- (1/pH) * ( (1-funexpit(cI[1]+cI[3])) * (1-funexpit(bV[1]+bV[3])) * funexpit(aH[1]) * (1-pC) +
                      (1-funexpit(cI[1]+cI[4]+cI[3])) * funexpit(bV[1]+bV[3]) * funexpit(aH[1]) * (1-pC) +
                      (1-funexpit(cI[1]+cI[2]+cI[3])) * (1-funexpit(bV[1]+bV[2]+bV[3])) * funexpit(aH[1]+aH[2]) * pC +
                      (1-funexpit(cI[1]+cI[2]+cI[4]+cI[3])) * funexpit(bV[1]+bV[2]+bV[3]) * funexpit(aH[1]+aH[2]) * pC ) 
pi1h0 <- 1-pi0h0
pi1h1 <- 1-pi0h1

#SANITY CHECK!:
#recompute P[I=1] by conditioning on H
pISanity <- pi1h0 * (1-pH) + pi1h1 * pH

#compute P[S=1|I=0]
pSnoinf <- funexpit(dS[1])*(1-pi1h0)*(1-pH)/(1-pI) + funexpit(dS[1]+dS[2])*(1-pi1h1)*pH/(1-pI)
#compute P[S=1|I=1]
pSinf <- funexpit(dS[1]+dS[3])*pi1h0*(1-pH)/pI + funexpit(dS[1]+dS[2]+dS[3])*pi1h1*pH/pI
#compute marginal: P[S=1]
pSSanity <- pSinf * pI + pSnoinf * (1-pI)

#SANITY CHECK!:
#compute P[S=1] again
pS <- funexpit(dS[1]) * (1-pi1h0) * (1-pH) +
  funexpit(dS[1]+dS[2]) * (1-pi1h1) * pH +
  funexpit(dS[1]+dS[3]) * pi1h0 * (1-pH) +
  funexpit(dS[1]+dS[2]+dS[3]) * pi1h1 * pH  
#-------------------------

#S' = selection variable for interaction model
#P[S'=1|H=1,I=0],P[S'=1|H=1,I=1])
dSauxint <- c(0.1,0.88) #plogis(sum(dS))=0.88, for comparison
#c(logit(P[S'=1|H=1,I=0]), log(OR_{IS'}(H=1))) = c(d1,d2)
dSint <- c(funlogit(dSauxint[1]), funlogOR(dSauxint[2],dSauxint[1]))

#---------------------
#A
#A|I=0 ~ Beta(alpha0,beta0)
alpha0 <- 2
beta0 <- 4
#A|I=1 ~ Beta(alpha1,beta1)
alpha1 <- 4
beta1 <- 2
#parameter vector
betap <- c(alpha0,beta0,alpha1,beta1)


#---------------------

#for scenario with another disease (Scenario III)

#O
#c(P[O=1|C=0,H=0], P[O=1|C=1,H=0], P[O=1|C=0,H=1])
coIaux <- c(0.2,0.3,0.1)
#c(logit(P[O=1|C=0,H=0,V=0]), log(OR_{CO}(H=0,V=0)), log(OR_{HO}(C=0,V=0))) = c(co1,co2,co3)
coI <- c(funlogit(coIaux[1]), funlogOR(coIaux[2],coIaux[1]), funlogOR(coIaux[3],coIaux[1]))

#SO
#P[S=1|H=0,I=0,O=1]
dSaux4 <- 0.2
dS4 <- funlogOR(dSaux4,dSaux[1])

#---------------------

#for scenario with arrow V->S (Scenario IV)

#SO
#P[S=1|H=0,I=0,V=1]
dSVaux <- 0.005
dSV <- funlogOR(dSVaux,dSaux[1])


#-------------------------
#compute marginal odds ratio: OR_{VI} diretly from parameters
trueMargOR <- exp(funlogOR(pIvacc,pInovacc))
trueMargVE <- 1-trueMargOR

#-------------------------
#compute marginal quantities from logistic regression coefficients: only for Baseline scenario!
getMarg <- function(pC,aH,bV,cI,dS){
  
  #compute P[H=1]
  pH <- funexpit(aH[1]+aH[2])*pC + funexpit(aH[1])*(1-pC)
  
  #compute P[V=1]
  pV <- funexpit(bV[1])*(1-funexpit(aH[1]))*(1-pC) +
    funexpit(bV[1]+bV[2])*(1-funexpit(aH[1]+aH[2]))*pC +
    funexpit(bV[1]+bV[3])*funexpit(aH[1])*(1-pC) +
    funexpit(bV[1]+bV[2]+bV[3])*funexpit(aH[1]+aH[2])*pC
  #compute P[V=1|C=1,H=1]
  pVc1h1 <- funexpit(bV[1]+bV[2]+bV[3])
  
  #compute P[I=1|V=0] from parameters
  pInovacc <- (1/(1-pV)) * ( funexpit(cI[1]) * (1-funexpit(aH[1])) * (1-pC) * (1-funexpit(bV[1])) + #(C=0,H=0)
                               funexpit(cI[1]+cI[2]) * (1-funexpit(aH[1]+aH[2])) *pC * (1-funexpit(bV[1]+bV[2])) + #(C=1,H=0)
                               funexpit(cI[1]+cI[3]) * funexpit(aH[1]) * (1-pC) * (1-funexpit(bV[1]+bV[3])) + #(C=0,H=1)
                               funexpit(cI[1]+cI[2]+cI[3]) * funexpit(aH[1]+aH[2]) * pC * (1-funexpit(bV[1]+bV[2]+bV[3])) ) #(C=1,H=1)
  #compute P[I=1|V=1] from parameters
  pIvacc <- (1/pV) * ( funexpit(cI[1]+cI[4]) * (1-funexpit(aH[1])) * (1-pC) * funexpit(bV[1]) +
                         funexpit(cI[1]+cI[2]+cI[4]) * (1-funexpit(aH[1]+aH[2])) * pC * funexpit(bV[1]+bV[2]) +
                         funexpit(cI[1]+cI[3]+cI[4]) * funexpit(aH[1]) * (1-pC) * funexpit(bV[1]+bV[3]) +
                         funexpit(cI[1]+cI[2]+cI[3]+cI[4]) * funexpit(aH[1]+aH[2]) * pC * funexpit(bV[1]+bV[2]+bV[3]) )
  #compute marginal: P[I=1]
  pI <- pIvacc * pV + pInovacc * (1-pV)
  
  
  #compute some auxiliary quantities:
  #compute P[I=0|H=0]
  pi0h0 <- (1/(1-pH)) * ( (1-funexpit(cI[1])) * (1-funexpit(bV[1])) * (1-funexpit(aH[1])) * (1-pC) +
                            (1-funexpit(cI[1]+cI[4])) * funexpit(bV[1]) * (1-funexpit(aH[1])) * (1-pC) +
                            (1-funexpit(cI[1]+cI[2])) * (1-funexpit(bV[1]+bV[2])) * (1-funexpit(aH[1]+aH[2])) * pC +
                            (1-funexpit(cI[1]+cI[2]+cI[4])) * funexpit(bV[1]+bV[2]) * (1-funexpit(aH[1]+aH[2])) * pC ) 
  #compute P[I=0|H=1]
  pi0h1 <- (1/pH) * ( (1-funexpit(cI[1]+cI[3])) * (1-funexpit(bV[1]+bV[3])) * funexpit(aH[1]) * (1-pC) +
                        (1-funexpit(cI[1]+cI[4]+cI[3])) * funexpit(bV[1]+bV[3]) * funexpit(aH[1]) * (1-pC) +
                        (1-funexpit(cI[1]+cI[2]+cI[3])) * (1-funexpit(bV[1]+bV[2]+bV[3])) * funexpit(aH[1]+aH[2]) * pC +
                        (1-funexpit(cI[1]+cI[2]+cI[4]+cI[3])) * funexpit(bV[1]+bV[2]+bV[3]) * funexpit(aH[1]+aH[2]) * pC ) 
  pi1h0 <- 1-pi0h0
  pi1h1 <- 1-pi0h1
  
  #compute P[S=1|I=0]
  pSnoinf <- funexpit(dS[1])*(1-pi1h0)*(1-pH)/(1-pI) + funexpit(dS[1]+dS[2])*(1-pi1h1)*pH/(1-pI)
  #compute P[S=1|I=1]
  pSinf <- funexpit(dS[1]+dS[3])*pi1h0*(1-pH)/pI + funexpit(dS[1]+dS[2]+dS[3])*pi1h1*pH/pI
  #compute marginal: P[S=1]
  pS <- pSinf * pI + pSnoinf * (1-pI)
  
  #compute other quantities of interest
  #compute P[I=1|S=1]
  pi1s1 <- pSinf * pI / pS
  #compute P[S=1,H=1]
  ps1ANDh1 <- funexpit(dS[1]+dS[2]) * pi0h1 * pH + funexpit(dS[1]+dS[2]+dS[3]) * pi1h1 * pH
  #compute P[I=1|S=1,H=1]
  pi1s1h1 <- funexpit(dS[1]+dS[2]+dS[3]) * pi1h1 * pH / ps1ANDh1
  #compute P[H=1|S=1]
  ph1s1 <- ps1ANDh1/pS
  
  #compute marginal odds ratio: OR_{VI} directly from parameters
  trueMargOR <- exp(funlogOR(pIvacc,pInovacc))
  trueMargVE <- 1-trueMargOR
  
  list(pMarg = c(pC=pC,pH=pH,pV=pV,pVc1h1=pVc1h1,pI=pI,pS=pS,
                 pSinf=pSinf,pSnoinf=pSnoinf,
                 pi1s1=pi1s1,
                 ps1ANDh1=ps1ANDh1,
                 pi1s1h1=pi1s1h1,
                 ph1s1=ph1s1,
                 pInovacc = pInovacc,
                 pIvacc=pIvacc,
                 trueMargOR=trueMargOR,trueMargVE=trueMargVE,
                 trueCondOR=exp(cI[4]),trueCondVE=1-exp(cI[4])),
       params=list(pC,aH,bV,cI,dS))
}


#-------------------------
#get cI[1] from cI[2:4] and P[I=1|V=0]
getC1 <- function(c1,c2,c3,c4,pC,aH,bV,pInovacc){
  pV <- funexpit(bV[1])*(1-funexpit(aH[1]))*(1-pC) +
    funexpit(bV[1]+bV[2])*(1-funexpit(aH[1]+aH[2]))*pC +
    funexpit(bV[1]+bV[3])*funexpit(aH[1])*(1-pC) +
    funexpit(bV[1]+bV[2]+bV[3])*funexpit(aH[1]+aH[2])*pC
  
  pInovacc - (1/(1-pV)) * ( funexpit(c1) * (1-funexpit(aH[1])) * (1-pC) * (1-funexpit(bV[1])) + #(C=0,H=0)
                              funexpit(c1+c2) * (1-funexpit(aH[1]+aH[2])) *pC * (1-funexpit(bV[1]+bV[2])) + #(C=1,H=0)
                              funexpit(c1+c3) * funexpit(aH[1]) * (1-pC) * (1-funexpit(bV[1]+bV[3])) + #(C=0,H=1)
                              funexpit(c1+c2+c3) * funexpit(aH[1]+aH[2]) * pC * (1-funexpit(bV[1]+bV[2]+bV[3])) ) #(C=1,H=1)
}



#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#generate parameter combinations from grid for parameters that will vary
#for fixed parameters, use values already defined above
c234 <- expand.grid(c2=incr,c3=incr01,c4=incr01)
c234$c1 <- sapply(c(1:ngrid^3),FUN=function(x){
  uniroot(getC1,c2=c234[x,1],c3=c234[x,2],c4=c234[x,3],pC=pC,aH=aH,bV=bV,pInovacc=pInovaccSET,interval=c(-10, 10))$root})
c234 <- c234[,c(4,1:3)]

#make a data frame to see what marginal probabilities you get with a certain set of parameters
fullsimMarg <- data.frame(matrix(0,ncol=18,nrow=ngrid^3))
colnames(fullsimMarg) <- c("pC","pH","pV","pVc1h1","pI","pS",
                           "pSinf","pSnoinf",
                           "pi1s1",
                           "ps1ANDh1",
                           "pi1s1h1",
                           "ph1s1",
                           "pInovacc",
                           "pIvacc",
                           "trueMargOR","trueMargVE","trueCondOR","trueCondVE")
for(i in 1:ngrid^3){
  #print(i)
  fullsimMarg[i,] <- getMarg(pC=pC,aH=aH,bV=bV,cI=unlist(c234[i,]),dS=dS)$pMarg
}

#-------------------------------------------------------------------------

# function to simulate population based on given parameters, and run logistic regression
# input: N, K, parameters
modMin <- function(N,K,pC,aH,bV,cI,dS,dSint,betap,coI,dS4,dSV){
  
  C <- rbinom(N,1,pC)
  H <- rbinom(N,1,plogis(aH[1]+aH[2]*C))
  V <- rbinom(N,1,plogis(bV[1]+bV[2]*C+bV[3]*H))
  I <- rbinom(N,1,plogis(cI[1]+cI[2]*C+cI[3]*H+cI[4]*V))
  S <- rbinom(N,1,plogis(dS[1]+dS[2]*H+dS[3]*I))
  
  A <- ifelse(I==0,rbeta(N,betap[1],betap[2]),rbeta(N,betap[3],betap[4]))
  SA <- rbinom(N,1,plogis(dS[1]+dS[2]*H+dS[3]*A))
  
  #selection variable for interaction
  SI <- ifelse(H==0,0,rbinom(N,1,plogis(dSint[1]*H+dSint[2]*H*I)))
  
  #model with another disease
  O <- rbinom(N,1,plogis(coI[1]+coI[2]*C+coI[3]*H))
  SO <-  rbinom(N,1,plogis(dS[1]+dS[2]*H+dS[3]*I+dS4*O))
  
  #model with arrow V->S
  SV <- rbinom(N,1,plogis(dS[1]+dS[2]*H+dS[3]*I+dSV*V))
  
  #------------------
  
  #concatenate data
  dat <- data.frame(C=C,H=H,V=V,I=I,S=S,SI=SI,A=A,SA=SA,O=O,SO=SO,SV=SV) 
  
  #------------------
  
  #tnd selection
  tndSel <- which(dat$S==1)
  #logistic regression within S=1, adjusted for C, unadjusted for H
  fit_unadj <- glm(I ~ V + C, family="binomial", data=dat[tndSel,]) 
  #logistic regression within S=1, adjusted for C, H
  fit_adj <- glm(I ~ V + C + H, family="binomial", data=dat[tndSel,]) 
  
  #---------------------
  
  #logistic regression for interaction model
  tndSelI <- which(dat$SI==1)
  fit_int <- glm(I ~ V + C , family="binomial", data=dat[tndSelI,]) 
  
  #-----------------------
  
  #tnd selection, severity affects probability of selection
  tndSelA <- which(dat$SA==1)
  #logistic regression within SA=1, adjusted for C, A, unadjusted for H
  fit_unadjA <- glm(I ~ V + C + A, family="binomial", data=dat[tndSelA,]) 
  #logistic regression within S=1, adjusted for C, A, H
  fit_adjA <- glm(I ~ V + C + H + A, family="binomial", data=dat[tndSelA,]) 
  
  #-----------------------
  
  #tnd selection, with another disease
  tndSelO <- which(dat$SO==1)
  #logistic regression within SA=1, adjusted for C, A, unadjusted for H
  fit_unadjO <- glm(I ~ V + C, family="binomial", data=dat[tndSelO,]) 
  #logistic regression within S=1, adjusted for C, A, H
  fit_adjO <- glm(I ~ V + C + H, family="binomial", data=dat[tndSelO,]) 
  
  #-----------------------
  
  #tnd selection, with another disease
  tndSelV <- which(dat$SV==1)
  #logistic regression within SA=1, adjusted for C, A, unadjusted for H
  fit_unadjV <- glm(I ~ V + C, family="binomial", data=dat[tndSelV,]) 
  #logistic regression within S=1, adjusted for C, A, H
  fit_adjV <- glm(I ~ V + C + H, family="binomial", data=dat[tndSelV,]) 
  
  #-----------------------
  #save input K
  Kin <- K
  #-----------------------
  
  #lower required K if it turns out you don't have enough cases with S=1
  K <- min(Kin,length(which(dat$I==1 & dat$S==1)))
  
  #logistic regression adjusted for C in case-control selection
  ccSel <- c(sample(x=which(dat$I==1 & dat$S==1),size=K),sample(x=which(dat$I==0),size=K))
  fit_cc <- glm(I ~ V + C , family="binomial", data=dat[ccSel,]) 
  
  #logistic regression adjusted for C in case-control selection, with cases selected from general population
  ccmSel <- c(sample(x=which(dat$I==1),size=K),sample(x=which(dat$I==0),size=K)) #sample cases from general population
  fit_ccm <- glm(I ~ V + C , family="binomial", data=dat[ccmSel,]) 
  
  #-----------------------
  
  K <- min(Kin,sum(dat$V==1),sum(dat$V==0))
  
  #logistic regression adjusted for C in cohort selection (no loss to follow-up)
  chSel <- c(sample(x=which(dat$V==1),size=K),sample(x=which(dat$V==0),size=K))
  fit_ch <- glm(I ~ V + C, family="binomial", data=dat[chSel,]) 
  
  #-----------------------
  
  list(pemp = c(sum(dat[,"H"])/N,sum(dat[,"V"])/N,sum(dat[,"I"])/N,sum(dat[,"S"])/N,
                sum(dat[,"SI"])/N,
                sum(dat[,"SA"])/N,
                sum(dat[,"O"])/N,sum(dat[,"SO"])/N,
                sum(dat[,"SV"])/N),
       
       coef_unadj = coef(summary(fit_unadj))["V",],
       coef_adj = coef(summary(fit_adj))["V",],
       tndSeltabl = table(dat$I[tndSel]),
       
       #------------
       
       coef_int = coef(summary(fit_int))["V",],
       intSeltab = table(dat$I[tndSelI]),
       
       #------------
       coef_unadjA = coef(summary(fit_unadjA))["V",],
       coef_adjA = coef(summary(fit_adjA))["V",],
       tndSeltablA = table(dat$I[tndSelA]),
       
       #------------
       coef_unadjO = coef(summary(fit_unadjO))["V",],
       coef_adjO = coef(summary(fit_adjO))["V",],
       tndSeltablO = table(dat$I[tndSelO]),
       
       #------------
       coef_unadjV = coef(summary(fit_unadjV))["V",],
       coef_adjV = coef(summary(fit_adjV))["V",],
       tndSeltablV = table(dat$I[tndSelV]),
       
       #--------------
       coef_cc = coef(summary(fit_cc))["V",],
       ccSeltab = table(dat$I[ccSel]),
       
       coef_ccm = coef(summary(fit_ccm))["V",],
       ccmSeltab = table(dat$I[ccmSel]),
       
       coef_ch = coef(summary(fit_ch))["V",],
       chSeltab = table(dat$I[chSel]),
       
       params=list(N,pC,aH,bV,cI,dS,dSint,betap,coI,dS4,dSV)
  )
}


#-----------------------------------------------------------------------------------
# Vary parameters and run simulation for grid
# Parameters to vary: c2 ("C->I"), c3 ("H->I") and c4 ("V->I")
#-----------------------------------------------------------------------------------
#simulate nrep populations of size N for each grid entry
fullsimGrid <- vector("list",ngrid^3)
for(i in 1:ngrid^3){
  print(i)
  fullsimGrid[[i]] <- lapply(c(1:nrep),FUN=function(x){modMin(N=N,K=K,pC=pC,aH=aH,bV=bV,cI=unlist(c234[i,]),dS=dS,dSint=dSint,betap=betap,
                                                              coI=coI,dS4=dS4,dSV=dSV)})
}

