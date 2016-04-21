install.packages("googlesheets")
install.packages("metafor")
library(googlesheets)
library(metafor)

###read in Google sheets data
gs_auth(new_user=TRUE)
my.sheet <- gs_title("Data Extraction")

my.data <- data.frame(gs_read_csv(ss=my.sheet, ws="Data", is.na(TRUE)))
names(my.data)


###Subset data into mammals and abundance
mammals<-rbind(my.data[which(my.data$Taxa=="Mammals" 
                             & my.data$Response=="Abundance"),],
               my.data[which(my.data$Taxa=="Rodents" 
                             & my.data$Response=="Abundance"),])


head(mammals)
dim(mammals)

wt2<-escalc(measure="SMD",m1i=TrtMean,m2i=CntrlMean,sd1i=sqrt(TrtVar),
            sd2i=sqrt(CntrlVar),n1i=n.Trt,n2i=n.Cntrl, 
            data=mammals,var.names=c("SMD","SMD_var"),digits=4)
head(wt2)
wt2$Taxa
wt2<-wt2[order(wt2$AuthorYear),]


###Histogram
hist(wt2$SMD, breaks=15, xlab="hedge's d", col=2)
abline(v=0,col=4,lty=3,lwd=5)
forest(wt2$SMD,wt2$SMD_var,slab=wt2$Genus,pch=19,main="Genus")


fsn(wt2$SMD,wt2$SMD_var,type="Rosenberg",alpha=.05,digits=4)
#returns a value of 0

##Check how many papers are used in the forest plot

forest(wt2$SMD,wt2$SMD_var,slab=wt2$AuthorYear,pch=19,main="PaperID")
Papers<-unique(wt2$AuthorYear)
write.csv(Papers,"MammalAbundancePapers.csv")

####Fixed Effects Model-excludes studies with "NA" Values

fixef.model <- rma(SMD, SMD_var, data=wt2, method = "FE")
summary(fixef.model)
plot(fixef.model)

forest(fixef.model,slab=wt2$Genus)

##Funnel plots were checked using a subset of wt2 points to ensure that they were working

funnel(fixef.model)
fixef.trim<-trimfill(fixef.model)
funnel(fixef.trim, main="Mammal Abundance Fixed Effect")  

####Random Effects Model-excludes studies with "NA" Values
ranef.model <- rma(SMD, SMD_var, data=wt2, method = "HE")

summary(ranef.model)
plot(ranef.model)

forest(ranef.model,slab=wt2$Genus)
funnel(ranef.model)

ranef.trim<-trimfill(ranef.model)
funnel(ranef.trim,main="Mammal Abundance Random Effect")


#####SDM as functions proposed:

##PaperID
##Genus
##VegetationType?


##instal MCMCglmm
install.packages("MCMCglmm")
require(MCMCglmm)


fixef.model <- rma(SMD ~ #PARAMETER, SMD_var, data=wt2, method = "FE")
                     summary(fixef.model)
                   plot(fixef.model,slab=wt2$Genus)
                   
                   
                   ranef.model <- rma(SMD~#PARAMETER, SMD_var, data=wt2, method = "HE")
                                        summary(ranef.model)
                                      plot(ranef.model)
                                      forest(ranef.model,slab=wt2$PaperID)
                                      ##back to problem at hand... comparing alternative approaches to inference
                                      
                                      ##maximum likelihood (also, I-T model selection)
                                      
                                      ml.model.reduced <- rma(SMD~ #PARAMETERS, SMD_var, data=wt2, method = "REML") #maximum-likelihood
                                                                summary(ml.model.reduced)
                                                              ml.model.full    <- rma(SMD ~ #ALLPARAMETERS, var.d, data=wt2, method = "REML") 
                                                                                        
                                                                                        ##Bayesian inference (not to be confused with empirical Bayes option in rma, not same thing)
                                                                                        
                                                                                        ##define priors (these are "minimally informative" priors for the covariances, see MCMCglmm man pages for help)
                                                                                        ##the only one of interest here is the "R" covariance, which ends of "translating" to a prior for the among-study
                                                                                        ##variance in effect sizes
                                                                                        
                                                                                        prior = list(R = list(V = 1, nu=0), G = list(list ( V = 1, n = 1, fix=1)))
                                                                                      
                                                                                      ##we're using MCMCglmm() to approximate the posterior probability distributions of the parameters of the model
                                                                                      ##that we're interested in... easier than WinBUGS or JAGS, less flexible
                                                                                      
                                                                                      
                                                                                      ###Will have to remove NAs before running this model
                                                                                      
                                                                                      Selected<-wt2[!is.na(wt2$SMD_var)& !is.na(wt2$Genus),]
                                                                                      names(Selected)
                                                                                      
                                                                                      bayes.model<-MCMCglmm(SMD~Genus, mev=Selected$SMD_var, data=Selected, nitt=15000, thin=10, burnin=5000, prior=prior)
                                                                                      summary(bayes.model)
                                                                                      
                                                                                      ##plots MCMC chains and posterior approximations
                                                                                      plot(bayes.model)
                                                                                      ##posterior probability distributions... use summary()
