### Some sample approaches to dealing with publication bias

fixef.model <- rma(hedge.d, var.d, data=dat, method = "FE")

trimandfill <- trimfill(fixef.model)

funnel(trimandfill) 

funnel(fixef.model)

fsn(yi=hedge.d,vi=var.d, data=dat, type="Rosenberg") #fail-safe N

plot(dat$year, dat$hedge.d) #scatterplot
