## opening pleasantries: packages

#install.packages("googlesheets")
#install.packages("forestplot")
#install.packages("metafor")
#install.packages("maps")
#install.packages("extrafont")
library(googlesheets)
library(forestplot)
library(metafor)
library(maps)
library(ggplot2)
library(extrafont)

#font_import()
loadfonts(device="win")

## read in Google sheets data

#gs_auth(new_user=TRUE)  ##get the token and authenticate
#my.sheet <- gs_title("Data Extraction")  ##
#my.data <- data.frame(gs_read_csv(ss=my.sheet, ws="Data", is.na(TRUE)))

setwd("C:/Users/db179/Google Drive/meta/")
#write.csv(my.data, file="meta20170206.csv")
my.data <- read.csv("meta20170924.csv")



## take care of NA/string problem

my.data$TrtMean[my.data$TrtMean=="NA"] <- NA
my.data$CntrlMean[my.data$CntrlMean=="NA"] <- NA
my.data$TrtVar[my.data$TrtVar=="NA"] <- NA
my.data$CntrlVar[my.data$CntrlVar=="NA"] <- NA
my.data$nTrt[my.data$nTrt=="NA"] <- NA
my.data$nCntrl[my.data$nCntrl=="NA"] <- NA

## treat appropriate fields as numeric

my.data$TrtMean <- as.numeric(as.character(as.factor(my.data$TrtMean)))
my.data$CntrlMean <- as.numeric(as.character(as.factor(my.data$CntrlMean)))
my.data$TrtVar <- as.numeric(as.character(as.factor(my.data$TrtVar)))
my.data$CntrlVar <- as.numeric(as.character(as.factor(my.data$CntrlVar)))
my.data$nTrt <- as.numeric(as.character(as.factor(my.data$nTrt)))
my.data$nCntrl <- as.numeric(as.character(as.factor(my.data$nCntrl)))



## store to csv and then extract again for road access!

#setwd("C:/Users/db179/Google Drive/meta/")
#write.csv(my.data, file="meta20170206.csv")
#my.data <- read.csv("meta20170924.csv")



######################
##MAPS

plot.new()
map("state")
map.axes()
points(my.data$LocationY ~ my.data$LocationX, pch=20, cex=1.25, col="blue")
text(my.data$LocationY ~ my.data$LocationX, labels=my.data$AuthorYear)
abline(v=-100)

################################################################################################

##Subset to richness data
richness<- subset(my.data,Response=="Richness")
head(richness)

#size of dataset
dim(richness)

#Calculate SD from Var
richness$sd_cont
richness$sd_cont<-(sqrt(richness$CntrlVar))
head(richness$sd_cont)

richness$sd_treat
richness$sd_treat<-(sqrt(richness$TrtVar))
richness$TrtMean <- as.numeric(richness$TrtMean)
richness$CntrlMean <- as.numeric(richness$CntrlMean)
head(richness$sd_treat)

#Remove NAs
richness <- richness[!is.na(richness$CntrlVar),]
richness <- richness[!is.na(richness$TrtMean),]
richness <- richness[!is.na(richness$nTrt),]
nrow(richness)


#Take standared mean difference and SMD variance
rich<-escalc(measure="SMD",m1i=(TrtMean+1), m2i=(CntrlMean+1),sd1i=sd_treat
              ,sd2i=sd_cont,n1i=nTrt,n2i=nCntrl, data=richness,
              var.names=c("SMD","SMD_var"),digits=4)

#Remove NAs
rich <- rich[!is.na(rich$SMD),]

#sort
rich <- rich[order(rich$Taxa),]

#include cis for graphing
#rich <- summary(rich)

##Histogram of headge's d
hist(rich$SMD, breaks=15, xlab="Hedge's d", col=2, main="Species Richness")
abline(v=0,col=4,lty=3,lwd=5)



##Forest plot
forest(rich$SMD,rich$SMD_var,slab=rich$Taxa,main="Species Richness",pch=19)

fixef.model <- rma(SMD, SMD_var, data=rich, method = "FE", mods= ~SampleDesign-1)
ranef.model <- rma(SMD, SMD_var, data=rich, method="HE", mods= ~SampleDesign-1)

##Fixed effects model of species richness with NA omitted
fixef.model <- rma(SMD, SMD_var, data=rich, method = "FE", mods= ~Taxa-1)

ranef.model <- rma(SMD, SMD_var, data=rich, method="HE", mods= ~Taxa-1)

##funnel plot with trimandfill
fixef.model <- rma(SMD, SMD_var, data=rich, method = "FE")
trimandfill <- trimfill(fixef.model)
par(family = "Times New Roman", cex=1.5)
funnel(fixef.model)
funnel(trimandfill)

res.mv <- rma.mv(SMD, SMD_var, mods= ~Taxa-1, random = ~ factor(StudyYear) | PaperID, data=rich)
res.mv <- rma.mv(SMD, SMD_var, mods= ~Taxa-1, random = ~ PaperID, data=rich)
res.mv <- rma.mv(SMD, SMD_var, mods= ~SampleDesign-1, random = ~ PaperID, data=rich)

forest(ranef.model,slab=rich$Taxa, cex=0.75)

dat <- data.frame(cite=c("Amphibians (2)","Birds (18)","Mammals (13)","Reptiles (1)"),yi=ranef.model$b,lowerci=ranef.model$ci.lb,upperci=ranef.model$ci.ub)

plot(rich$nTrt, rich$SMD, ylab="Hedge's d", xlab="sample size (sites)", pch=19)
abline(h=0, lty=2)
abline(h=-0.2619)

fsn(yi=SMD,vi=SMD_var,data=rich, type="Rosenberg")


##this yields a decent summary figure using ggplot2

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family="Times New Roman", size=24),
        legend.position='none')

ggplot(data=rich, aes(rich$SMD)) + 
	geom_histogram(binwidth=0.5) + 
	geom_vline(xintercept = 0, size=2) + 
	scale_x_continuous(name='Standardized Mean Difference (Hedges\' d)') +
	ylab('Frequency') +
	apatheme

p = ggplot(dat, aes(y=cite, x=yi, xmin=lowerci, xmax=upperci)) +
	geom_point(color = 'black', shape=18, size=12) + 
	geom_errorbarh(height=.2, size=1) + 
	ylab('Taxonomic Class') +
	scale_x_continuous(limits=c(-5,2), name='Standardized Mean Difference (Hedges\' d)') +
	geom_vline(xintercept=0, color='black', linetype='dashed' )+
	apatheme

p

##forest plot overall

richsum <- summary(rich)

dat2 <- data.frame(cite=richsum$AuthorYear,yi=richsum$SMD,lowerci=richsum$ci.lb,upperci=richsum$ci.ub,tester=richsum$Taxa)

p = ggplot(dat2, aes(y=cite, x=yi, xmin=lowerci, xmax=upperci, shape=tester)) +
	geom_point(color = 'black', shape=18, size=12) + 
	geom_errorbarh(height=.2, size=1) + 
	ylab('Taxonomic Class') +
	scale_x_continuous(limits=c(-8,4), breaks = c(-8,-4,0,4), name='Standardized Mean Difference (g)') +
	geom_vline(xintercept=0, color='black', linetype='dashed' ) +
	facet_grid(tester~., scales= 'free')+
	apatheme

p



###############################################################################






## create a new table "complete" which contains all complete records

complete <- my.data[which(is.na(my.data$TrtMean)==F &
		  is.na(my.data$CntrlMean)==F &
		  is.na(my.data$TrtVar)==F &
		  is.na(my.data$CntrlVar)==F &
		  is.na(my.data$nTrt)==F &
		  is.na(my.data$nCntrl)==F),]

## create a new table "abundance" which only contains abundance studies

abundance <- complete[which(complete$Response=="Abundance" | 
			     complete$Response=="Density" |
			     complete$Response=="Nest Density" | 
			     complete$Response=="Relative Abundance"),]

#abundance$TrtMean <- as.numeric(abundance$TrtMean)
#abundance$CntrlMean <- as.numeric(abundance$CntrlMean)
#abundance$TrtVar <- as.numeric(abundance$TrtVar)
#abundance$CntrlVar <- as.numeric(abundance$CntrlVar)


wt2<-escalc(measure="SMD",m1i=TrtMean,m2i=CntrlMean,sd1i=sqrt(TrtVar),
            sd2i=sqrt(CntrlVar),n1i=nTrt,n2i=nCntrl, 
            data=abundance,var.names=c("SMD","SMD_var"),digits=4, na.action=na.exclude)

#Remove NAs
wt2 <- wt2[!is.na(wt2$SMD),]

#sort
wt2 <- wt2[order(wt2$Taxa),]

fixef.model <- rma(SMD, SMD_var, data=wt2, method = "FE")

plot(wt2$nTrt, wt2$SMD, ylab="Hedge's d", xlab="sample size (sites)", pch=19)
abline(h=0, lty=2)
abline(h=-0.1892)

fsn(yi=SMD,vi=SMD_var,data=wt2, type="Rosenberg")

wt2sum <- summary(wt2)

dat2 <- data.frame(cite=wt2sum$AuthorYear,yi=wt2sum$SMD,lowerci=wt2sum$ci.lb,upperci=wt2sum$ci.ub,tester=wt2sum$Taxa)

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family="Times New Roman", size=18),
        legend.position='none')

p = ggplot(dat2, aes(y=cite, x=yi, xmin=lowerci, xmax=upperci, shape=tester)) +
	geom_point(color = 'black', shape=18, size=8) + 
	geom_errorbarh(height=.2, size=1) + 
	ylab('Taxonomic Class') +
	scale_x_continuous(limits=c(-15,10), breaks = c(-15,-10,-5,0,5,10), name='Standardized Mean Difference (g)') +
	geom_vline(xintercept=0, color='black', linetype='dashed' ) +
	facet_grid(tester~., scales= 'free')+
	apatheme

p


wt2$Riparian[which(is.na(wt2$Riparian))] <- 0


forest(wt2$SMD,wt2$SMD_var,slab=wt2$AuthorYear,pch=19,cex=.5,main="Genus")
ranef.model <- rma(SMD, SMD_var, data=wt2, mods= ~SampleDesign*Taxa, method = "HE")
res.mv <- rma(SMD, SMD_var, data=wt2, mods= ~Riparian-1, method = "HE")
null <- rma(SMD~1, mod=PaperID, SMD_var, data=wt2, method = "ML")
res.mv <- rma.mv(SMD, SMD_var, mods= ~Taxa/Genus-1, random = ~ PaperID, data=wt2, method="ML")
res.mv <- rma(SMD, SMD_var, data=wt2, mods= ~SampleDesign*Taxa, method = "HE")



datx <- data.frame(cite=c("Non-riparian","Riparian"),yi=res.mv$b,lowerci=res.mv$ci.lb,upperci=res.mv$ci.ub)








dat <- data.frame(cite=c("Amphibians (15)","Birds (129)","Mammals (72)","Reptiles (29)"),yi=res.mv$b,lowerci=res.mv$ci.lb,upperci=res.mv$ci.ub)


##this yields a decent summary figure using ggplot2

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family="Times New Roman", size=24),
        legend.position='none')

ggplot(data=wt2, aes(wt2$SMD)) + 
	geom_histogram(binwidth=0.5) + 
	geom_vline(xintercept = 0, size=2) + 
	scale_x_continuous(name='Standardized Mean Difference (Hedges\' d)') +
	ylab('Frequency') +
	apatheme

p = ggplot(datx, aes(y=cite, x=yi, xmin=lowerci, xmax=upperci)) +
	geom_point(color = 'black', shape=18, size=12) + 
	geom_errorbarh(height=.2, size=1) + 
	ylab('') +
	scale_x_continuous(limits=c(-1.6,1), breaks = c(-1,0,1), name='Standardized Mean Difference (Hedges\' d)') +
	geom_vline(xintercept=0, color='black', linetype='dashed' )+
	apatheme

p




















