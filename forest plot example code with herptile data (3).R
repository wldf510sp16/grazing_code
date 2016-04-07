## code taken from Nathan's mammal code to generate forest plots

install.packages("googlesheets")
install.packages("metafor")
library(googlesheets)
library(metafor)

###read in Google sheets data
gs_auth(new_user=TRUE)
my.sheet <- gs_title("Data Extraction")

my.data <- data.frame(gs_read_csv(ss=my.sheet, ws="Data", is.na(TRUE)))
names(my.data)


###Subset data into reptiles and amphibians
herptiles<-rbind(my.data[which(my.data$Taxa=="Amphibians" 
                             & my.data$Response=="Abundance"),],
               my.data[which(my.data$Taxa=="Reptiles" 
                             & my.data$Response=="Abundance"),])


head(herptiles)
dim(herptiles)

wt2<-escalc(measure="SMD",m1i=TrtMean,m2i=CntrlMean,sd1i=sqrt(TrtVar),
            sd2i=sqrt(CntrlVar),n1i=n.Trt,n2i=n.Cntrl, 
            data=herptiles,var.names=c("SMD","SMD_var"),digits=4)
head(wt2)
wt2$Taxa
wt2<-wt2[order(wt2$Genus),]
names(wt2)

###Histogram
hist(wt2$SMD, breaks=15, xlab="hedge's d", col=2)
abline(v=0,col=4,lty=3,lwd=5)
forest(wt2$SMD,wt2$SMD_var,slab=wt2$Species,pch=19)



####Fixed Effects Model-excludes studies with "NA" Values

fixef.model <- rma(SMD, SMD_var, data=wt2, method = "FE")
summary(fixef.model)
dev.new()
plot(fixef.model)
forest(fixef.model)

####Random Effects Model-excludes studies with "NA" Values
ranef.model <- rma(SMD, SMD_var, data=wt2, method = "HE")

summary(ranef.model)
dev.new()
plot(ranef.model)
forest(ranef.model)