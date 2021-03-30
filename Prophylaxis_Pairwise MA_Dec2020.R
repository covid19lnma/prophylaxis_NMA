#############################################################################################
### Random effects pairwise Bayesian Meta-Analysis, using Turner Priors, 0.5 and 1.0 half normal priors
### Comparison also with FE and RE frequentist models
###
### Original code from LG, with modifications by AQ to facilitate analysis.
### Outputs true to Long's for ease of GRADING

### Jan 13, 2021
#############################################################################################

# load necessary libraries:
require("metafor")
require("bayesmeta")
library(dplyr)
library(forestplot)
library(meta)
library(scales)

setwd("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Data")

#############################################################################################
#### Mortality
#############################################################################################

# load data:
data <- read.csv("mortality - wide data format.csv",
                 header = TRUE)

#### Hydroxychloroquine vs. Standard care/Placebo
####################################################

data_hcq <- filter(data,
	t1 == "hydroxychloroquine" &
	t2 == "a_placebo/standard care")

effsize1_mort <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize1_mort2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_hcq)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("all-cause mortality", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner1 <- bayesmeta(effsize1_mort, tau.prior=TP$dprior)
bm.hnorm110 <- bayesmeta(effsize1_mort, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105 <- bayesmeta(effsize1_mort, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner1.2 <- bayesmeta(effsize1_mort2, tau.prior=TP$dprior)
bm.hnorm110.2 <- bayesmeta(effsize1_mort2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105.2 <- bayesmeta(effsize1_mort2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed1 <- rma.uni(effsize1_mort, method="FE")
rma.random.DL1 <- rma.uni(effsize1_mort, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed1.2 <- rma.uni(effsize1_mort2, method="FE")
rma.random.DL1.2 <- rma.uni(effsize1_mort2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates1 <- rbind("HNorm(0.5)" = c(bm.hnorm105$summary[2,1:2],   bm.hnorm105$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1$summary[2,1:2],bm.Turner1$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1$tau2.fix),rma.fixed1$b,rma.fixed1$ci.lb,rma.fixed1$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1$tau2), rma.random.DL1$b, rma.random.DL1$ci.lb,rma.random.DL1$ci.ub))

# estimates where RD is the measure of effects
estimates1.2 <- rbind("HNorm(0.5)" = c(bm.hnorm105.2$summary[2,1:2],   bm.hnorm105.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1.2$summary[2,1:2],bm.Turner1.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1.2$tau2.fix),rma.fixed1.2$b,rma.fixed1$ci.lb,rma.fixed1.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1.2$tau2), rma.random.DL1.2$b, rma.random.DL1.2$ci.lb,rma.random.DL1.2$ci.ub))


#### Ivermectin vs. Standard care/Placebo
####################################################

data_iv <- filter(data,
	t1 == "ivermectin" &
	t2 == "a_placebo/standard care")

effsize2_mort <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_iv)

effsize2_mort2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_iv)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("all-cause mortality", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner2 <- bayesmeta(effsize2_mort, tau.prior=TP$dprior)
bm.hnorm210 <- bayesmeta(effsize2_mort, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205 <- bayesmeta(effsize2_mort, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner2.2 <- bayesmeta(effsize2_mort2, tau.prior=TP$dprior)
bm.hnorm210.2 <- bayesmeta(effsize2_mort2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205.2 <- bayesmeta(effsize2_mort2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed2 <- rma.uni(effsize2_mort, method="FE")
rma.random.DL2 <- rma.uni(effsize2_mort, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed2.2 <- rma.uni(effsize2_mort2, method="FE")
rma.random.DL2.2 <- rma.uni(effsize2_mort2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates2 <- rbind("HNorm(0.5)" = c(bm.hnorm205$summary[2,1:2],   bm.hnorm205$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2$summary[2,1:2],bm.Turner2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2$tau2.fix),rma.fixed2$b,rma.fixed2$ci.lb,rma.fixed2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2$tau2), rma.random.DL2$b, rma.random.DL2$ci.lb,rma.random.DL2$ci.ub))

# estimates where RD is the measure of effects
estimates2.2 <- rbind("HNorm(0.5)" = c(bm.hnorm205.2$summary[2,1:2],   bm.hnorm205.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2.2$summary[2,1:2],bm.Turner2.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2.2$tau2.fix),rma.fixed2.2$b,rma.fixed2$ci.lb,rma.fixed2.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2.2$tau2), rma.random.DL2.2$b, rma.random.DL2.2$ci.lb,rma.random.DL.2$ci.ub))

#### Results & Forest plot
####################################################

sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/Mortality - All Comparisons.txt")
print("   ")
print("   ")
print("   ")
print("Hydroxychloroquine vs. Standard care/Placebo")
print("   ")
print("   ")
print("Summary measure = OR")
print(estimates1)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates1.2)
print("   ")
print("   ")
print("Ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates2)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates2.2)
print("   ")
print("   ")
sink()


pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/Mortality - All Comparisons.pdf",width=8,height = 5)

### Hydroxychloroquine vs. Standard Care/Placebo, OR
yrange1 <- c(-7 - nrow(effsize1_mort), 1)
forest.default(effsize1_mort$yi, vi = effsize1_mort$vi, refline = 0,
               rows = seq(-2, -length(effsize1_mort$yi) - 1, by = -1),width=0,
               xlim = c(-10,10), ylim = yrange1, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize1_mort[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("Mortality: Hydroxychloroquine vs. Standard Care/Placebo - OR", line = 0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1))
  addpoly(estimates1[i,"mu"], ci.lb=estimates1[i,3], ci.ub=estimates1[i,4], atransf=exp,
          mlab=row.names(estimates1)[i], rows=yrange1[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Hydroxychloroquine vs. Standard Care/Placebo, RD
yrange1.2 <- c(-7 - nrow(effsize1.2), 1)
forest.default(effsize1.2$yi, vi = effsize1.2$vi, refline = 0,
               rows = seq(-2, -length(effsize1.2$yi) - 1, by = -1),width=0,
               xlim = c(-10,10), ylim = yrange1.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize1.2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               #atransf=exp,
               digits=2)

title("Mortality: Hydroxychloroquine vs. Standard Care/Placebo - RD", line=-1)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1.2))
  addpoly(estimates1.2[i,"mu"], ci.lb=estimates1.2[i,3], ci.ub=estimates1.2[i,4], #atransf=exp,
          mlab=row.names(estimates1.2)[i], rows=yrange1.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, OR
yrange2 <- c(-6 - nrow(effsize2_mort), 2)
forest.default(effsize2_mort$yi, vi = effsize2_mort$vi, refline = 0,
               rows = seq(-2, -length(effsize2_mort$yi) - 1, by = -1),width=0,
               xlim = c(-10,10), ylim = yrange2, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize2_mort[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("Mortality: Ivermectin vs. Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2))
  addpoly(estimates2[i,"mu"], ci.lb=estimates2[i,3], ci.ub=estimates2[i,4], atransf=exp,
          mlab=row.names(estimates2)[i], rows=yrange2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, RD
yrange2.2 <- c(-6 - nrow(effsize2_mort2), 2)
forest.default(effsize2_mort2$yi, vi = effsize2_mort2$vi, refline = 0,
               rows = seq(-2, -length(effsize2_mort2$yi) - 1, by = -1),width=0,
               xlim = c(-10, 10),
               #alim = c(-5, 5),
               ylim = yrange2.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize2_mort2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               #atransf=exp,
               digits=2)

title("Mortality: Ivermectin vs. Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2.2))
  addpoly(estimates2.2[i,"mu"], ci.lb=estimates2.2[i,3], ci.ub=estimates2.2[i,4], #atransf=exp,
          mlab=row.names(estimates2.2)[i], rows=yrange2.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

dev.off()

#############################################################################################
#### Lab-confirmed COVID-19
#############################################################################################

# load data:
data <- read.csv("COVID-19 lab confirmed - wide data format.csv",
                 header = TRUE)

#### Hydroxychloroquine vs. Standard care/Placebo
####################################################

data_hcq <- filter(data,
	t1 == "hydroxychloroquine" &
	t2 == "a_placebo/standard care")

effsize1_labinf <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize1_labinf2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_hcq)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("infection / onset of new disease", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner1 <- bayesmeta(effsize1_labinf, tau.prior=TP$dprior)
bm.hnorm110 <- bayesmeta(effsize1_labinf, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105 <- bayesmeta(effsize1_labinf, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner1.2 <- bayesmeta(effsize1_labinf2, tau.prior=TP$dprior)
bm.hnorm110.2 <- bayesmeta(effsize1_labinf2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105.2 <- bayesmeta(effsize1_labinf2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed1 <- rma.uni(effsize1_labinf, method="FE")
rma.random.DL1 <- rma.uni(effsize1_labinf, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed1.2 <- rma.uni(effsize1_labinf2, method="FE")
rma.random.DL1.2 <- rma.uni(effsize1_labinf2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates1 <- rbind("HNorm(0.5)" = c(bm.hnorm105$summary[2,1:2],   bm.hnorm105$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1$summary[2,1:2],bm.Turner1$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1$tau2.fix),rma.fixed1$b,rma.fixed1$ci.lb,rma.fixed1$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1$tau2), rma.random.DL1$b, rma.random.DL1$ci.lb,rma.random.DL1$ci.ub))

# estimates where RD is the measure of effects
estimates1.2 <- rbind("HNorm(0.5)" = c(bm.hnorm105.2$summary[2,1:2],   bm.hnorm105.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1.2$summary[2,1:2],bm.Turner1.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1.2$tau2.fix),rma.fixed1.2$b,rma.fixed1$ci.lb,rma.fixed1.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1.2$tau2), rma.random.DL1.2$b, rma.random.DL1.2$ci.lb,rma.random.DL1.2$ci.ub))


#### Ivermectin vs. Standard care/Placebo
####################################################

data_iv <- filter(data,
	t1 == "ivermectin" &
	t2 == "a_placebo/standard care")

effsize2_labinf <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_iv)

effsize2_labinf2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_iv)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("infection / onset of new disease", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner2 <- bayesmeta(effsize2_labinf, tau.prior=TP$dprior)
bm.hnorm210 <- bayesmeta(effsize2_labinf, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205 <- bayesmeta(effsize2_labinf, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner2.2 <- bayesmeta(effsize2_labinf2, tau.prior=TP$dprior)
bm.hnorm210.2 <- bayesmeta(effsize2_labinf2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205.2 <- bayesmeta(effsize2_labinf2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed2 <- rma.uni(effsize2_labinf, method="FE")
rma.random.DL2 <- rma.uni(effsize2_labinf, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed2.2 <- rma.uni(effsize2_labinf2, method="FE")
rma.random.DL2.2 <- rma.uni(effsize2_labinf2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates2 <- rbind("HNorm(0.5)" = c(bm.hnorm205$summary[2,1:2],   bm.hnorm205$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2$summary[2,1:2],bm.Turner2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2$tau2.fix),rma.fixed2$b,rma.fixed2$ci.lb,rma.fixed2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2$tau2), rma.random.DL2$b, rma.random.DL2$ci.lb,rma.random.DL2$ci.ub))

# estimates where RD is the measure of effects
estimates2.2 <- rbind("HNorm(0.5)" = c(bm.hnorm205.2$summary[2,1:2],   bm.hnorm205.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2.2$summary[2,1:2],bm.Turner2.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2.2$tau2.fix),rma.fixed2.2$b,rma.fixed2$ci.lb,rma.fixed2.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2.2$tau2), rma.random.DL2.2$b, rma.random.DL2.2$ci.lb,rma.random.DL2.2$ci.ub))


#### Iota-carrageenan, ivermectin vs. Standard care/Placebo
####################################################

data_iciv <- filter(data,
	t1 == "iota-carrageenan, ivermectin" &
	t2 == "a_placebo/standard care")

effsize3_labinf <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_iciv)

effsize3_labinf2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_iciv)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("infection / onset of new disease", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner3 <- bayesmeta(effsize3_labinf, tau.prior=TP$dprior)
bm.hnorm310 <- bayesmeta(effsize3_labinf, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm305 <- bayesmeta(effsize3_labinf, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner3.2 <- bayesmeta(effsize3_labinf2, tau.prior=TP$dprior)
bm.hnorm310.2 <- bayesmeta(effsize3_labinf2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm305.2 <- bayesmeta(effsize3_labinf2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed3 <- rma.uni(effsize3_labinf, method="FE")
rma.random.DL3 <- rma.uni(effsize3_labinf, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed3.2 <- rma.uni(effsize3_labinf2, method="FE")
rma.random.DL3.2 <- rma.uni(effsize3_labinf2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates3 <- rbind("HNorm(0.5)" = c(bm.hnorm305$summary[2,1:2],   bm.hnorm305$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner3$summary[2,1:2],bm.Turner3$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed3$tau2.fix),rma.fixed3$b,rma.fixed3$ci.lb,rma.fixed3$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL3$tau2), rma.random.DL3$b, rma.random.DL3$ci.lb,rma.random.DL3$ci.ub))

# estimates where RD is the measure of effects
estimates3.2 <- rbind("HNorm(0.5)" = c(bm.hnorm305.2$summary[2,1:2],   bm.hnorm305.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner3.2$summary[2,1:2],bm.Turner3.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed3.2$tau2.fix),rma.fixed3.2$b,rma.fixed3$ci.lb,rma.fixed3.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL3.2$tau2), rma.random.DL3.2$b, rma.random.DL3.2$ci.lb,rma.random.DL3.2$ci.ub))


#### Results & Forest plot
####################################################

sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/COVID-19 Lab confirmed - All Comparisons.txt")
print("   ")
print("   ")
print("   ")
print("Hydroxychloroquine vs. Standard care/Placebo")
print("   ")
print("   ")
print("Summary measure = OR")
print(estimates1)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates1.2)
print("   ")
print("   ")
print("Ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates2)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates2.2)
print("   ")
print("   ")
print("   ")
print("   ")
print("Iota-carrageenan, ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates3)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates3.2)
print("   ")
print("   ")
sink()

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/COVID-19 Lab confirmed - All Comparisons.pdf",width=8,height = 5)

### Hydroxychloroquine vs. Standard Care/Placebo, OR
yrange1 <- c(-7 - nrow(effsize1_labinf), 1)
forest.default(effsize1_labinf$yi, vi = effsize1_labinf$vi, refline = 0,
               rows = seq(-2, -length(effsize1_labinf$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize1_labinf[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("COVID-19 Lab Confirmed: Hydroxychloroquine vs. Standard Care/Placebo - OR", line= 0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1))
  addpoly(estimates1[i,"mu"], ci.lb=estimates1[i,3], ci.ub=estimates1[i,4], atransf=exp,
          mlab=row.names(estimates1)[i], rows=yrange1[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Hydroxychloroquine vs. Standard Care/Placebo, RD
yrange1.2 <- c(-7 - nrow(effsize1_labinf2), 1)
forest.default(effsize1_labinf2$yi, vi = effsize1_labinf2$vi, refline = 0,
               rows = seq(-2, -length(effsize1_labinf2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize1_labinf2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("COVID-19 Lab Confirmed: Hydroxychloroquine vs. Standard Care/Placebo - RD", line= 0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1.2))
  addpoly(estimates1.2[i,"mu"], ci.lb=estimates1.2[i,3], ci.ub=estimates1.2[i,4],
          mlab=row.names(estimates1.2)[i], rows=yrange1.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, OR
yrange2 <- c(-6 - nrow(effsize2_labinf), 2)
forest.default(effsize2_labinf$yi, vi = effsize2_labinf$vi, refline = 0,
               rows = seq(-2, -length(effsize2_labinf$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange2, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize2_labinf[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("COVID-19 Lab Confirmed: Ivermectin vs. Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2))
  addpoly(estimates2[i,"mu"], ci.lb=estimates2[i,3], ci.ub=estimates2[i,4], atransf=exp,
          mlab=row.names(estimates2)[i], rows=yrange2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, RD
yrange2.2 <- c(-6 - nrow(effsize2_labinf2), 2)
forest.default(effsize2_labinf2$yi, vi = effsize2_labinf2$vi, refline = 0,
               rows = seq(-2, -length(effsize2_labinf2$yi) - 1, by = -1),width=0,
               xlim = c(-10,10),
               alim = c(-2, 2),
               ylim = yrange2.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize2_labinf2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("COVID-19 Lab Confirmed: Ivermectin vs. Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2.2))
  addpoly(estimates2.2[i,"mu"], ci.lb=estimates2.2[i,3], ci.ub=estimates2.2[i,4],
          mlab=row.names(estimates2.2)[i], rows=yrange2.2[1]+5-i, col=colvec[i],cex=1.5,width =0)


### Iota-carrageenan, ivermectin  vs. Standard Care/Placebo, OR
yrange3 <- c(-6 - nrow(effsize3_labinf), 2)
forest.default(effsize3_labinf$yi, vi = effsize3_labinf$vi, refline = 0,
               rows = seq(-2, -length(effsize3_labinf$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange3, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize3_labinf[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("COVID-19 Lab Confirmed: Iota-carrageenan, ivermectin vs. Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates3))
  addpoly(estimates3[i,"mu"], ci.lb=estimates3[i,3], ci.ub=estimates3[i,4], atransf=exp,
          mlab=row.names(estimates3)[i], rows=yrange3[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Iota-carrageenan, ivermectin vs. Standard Care/Placebo, RD
yrange3.2 <- c(-6 - nrow(effsize3_labinf2), 1)
forest.default(effsize3_labinf2$yi, vi = effsize3_labinf2$vi, refline = 0,
               rows = seq(-2, -length(effsize3_labinf2$yi) - 1, by = -1),width=0,
               xlim = c(-4,4),
               alim = c(-1, 1),
               ylim = yrange3.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize3_labinf2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("COVID-19 Lab Confirmed: Iota-carrageenan, ivermectin vs. Standard Care/Placebo, RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates3.2))
  addpoly(estimates3.2[i,"mu"], ci.lb=estimates3.2[i,3], ci.ub=estimates3.2[i,4],
          mlab=row.names(estimates3.2)[i], rows=yrange3.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

dev.off()

#############################################################################################
#### Lab-confirmed and suspected COVID-19
#############################################################################################

# load data:
data <- read.csv("COVID-19 lab confirmed & suspected - wide data format.csv",
                 header = TRUE)

#### Hydroxychloroquine vs. Standard care/Placebo
####################################################

data_hcq <- filter(data,
	t1 == "hydroxychloroquine" &
	t2 == "a_placebo/standard care")

effsize1_labsinf <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize1_labsinf2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_hcq)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("infection / onset of new disease", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner1 <- bayesmeta(effsize1_labsinf, tau.prior=TP$dprior)
bm.hnorm110 <- bayesmeta(effsize1_labsinf, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105 <- bayesmeta(effsize1_labsinf, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner1.2 <- bayesmeta(effsize1_labsinf2, tau.prior=TP$dprior)
bm.hnorm110.2 <- bayesmeta(effsize1_labsinf2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105.2 <- bayesmeta(effsize1_labsinf2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed1 <- rma.uni(effsize1_labsinf, method="FE")
rma.random.DL1 <- rma.uni(effsize1_labsinf, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed1.2 <- rma.uni(effsize1_labsinf2, method="FE")
rma.random.DL1.2 <- rma.uni(effsize1_labsinf2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates1 <- rbind("HNorm(0.5)" = c(bm.hnorm105$summary[2,1:2],   bm.hnorm105$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1$summary[2,1:2],bm.Turner1$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1$tau2.fix),rma.fixed1$b,rma.fixed1$ci.lb,rma.fixed1$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1$tau2), rma.random.DL1$b, rma.random.DL1$ci.lb,rma.random.DL1$ci.ub))

# estimates where RD is the measure of effects
estimates1.2 <- rbind("HNorm(0.5)" = c(bm.hnorm105.2$summary[2,1:2],   bm.hnorm105.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1.2$summary[2,1:2],bm.Turner1.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1.2$tau2.fix),rma.fixed1.2$b,rma.fixed1$ci.lb,rma.fixed1.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1.2$tau2), rma.random.DL1.2$b, rma.random.DL1.2$ci.lb,rma.random.DL1.2$ci.ub))


#### Ivermectin vs. Standard care/Placebo
####################################################

data_iv <- filter(data,
	t1 == "ivermectin" &
	t2 == "a_placebo/standard care")

effsize2_labsinf <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_iv)

effsize2_labsinf2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_iv)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("infection / onset of new disease", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner2 <- bayesmeta(effsize2_labsinf, tau.prior=TP$dprior)
bm.hnorm210 <- bayesmeta(effsize2_labsinf, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205 <- bayesmeta(effsize2_labsinf, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner2.2 <- bayesmeta(effsize2_labsinf2, tau.prior=TP$dprior)
bm.hnorm210.2 <- bayesmeta(effsize2_labsinf2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205.2 <- bayesmeta(effsize2_labsinf2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed2 <- rma.uni(effsize2_labsinf, method="FE")
rma.random.DL2 <- rma.uni(effsize2_labsinf, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed2.2 <- rma.uni(effsize2_labsinf2, method="FE")
rma.random.DL2.2 <- rma.uni(effsize2_labsinf2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates2 <- rbind("HNorm(0.5)" = c(bm.hnorm205$summary[2,1:2],   bm.hnorm205$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2$summary[2,1:2],bm.Turner2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2$tau2.fix),rma.fixed2$b,rma.fixed2$ci.lb,rma.fixed2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2$tau2), rma.random.DL2$b, rma.random.DL2$ci.lb,rma.random.DL2$ci.ub))

# estimates where RD is the measure of effects
estimates2.2 <- rbind("HNorm(0.5)" = c(bm.hnorm205.2$summary[2,1:2],   bm.hnorm205.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2.2$summary[2,1:2],bm.Turner2.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2.2$tau2.fix),rma.fixed2.2$b,rma.fixed2$ci.lb,rma.fixed2.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2.2$tau2), rma.random.DL2.2$b, rma.random.DL2.2$ci.lb,rma.random.DL2.2$ci.ub))


#### Iota-carrageenan, ivermectin vs. Standard care/Placebo
##############################################################

data_iciv <- filter(data,
	t1 == "iota-carrageenan, ivermectin" &
	t2 == "a_placebo/standard care")

effsize3_labsinf <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize3_labsinf2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_hcq)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("infection / onset of new disease", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner3 <- bayesmeta(effsize3_labsinf, tau.prior=TP$dprior)
bm.hnorm310 <- bayesmeta(effsize3_labsinf, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm305 <- bayesmeta(effsize3_labsinf, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner3.2 <- bayesmeta(effsize3_labsinf2, tau.prior=TP$dprior)
bm.hnorm310.2 <- bayesmeta(effsize3_labsinf2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm305.2 <- bayesmeta(effsize3_labsinf2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed3 <- rma.uni(effsize3_labsinf, method="FE")
rma.random.DL3 <- rma.uni(effsize3_labsinf, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed3.2 <- rma.uni(effsize3_labsinf2, method="FE")
rma.random.DL3.2 <- rma.uni(effsize3_labsinf2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates3 <- rbind("HNorm(0.5)" = c(bm.hnorm305$summary[2,1:2],   bm.hnorm305$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner3$summary[2,1:2],bm.Turner3$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed3$tau2.fix),rma.fixed3$b,rma.fixed3$ci.lb,rma.fixed3$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL3$tau2), rma.random.DL3$b, rma.random.DL3$ci.lb,rma.random.DL3$ci.ub))

# estimates where RD is the measure of effects
estimates3.2 <- rbind("HNorm(0.5)" = c(bm.hnorm305.2$summary[2,1:2],   bm.hnorm305.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner3.2$summary[2,1:2],bm.Turner3.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed3.2$tau2.fix),rma.fixed3.2$b,rma.fixed3$ci.lb,rma.fixed3.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL3.2$tau2), rma.random.DL3.2$b, rma.random.DL3.2$ci.lb,rma.random.DL3.2$ci.ub))


#### Results & Forest plot
####################################################

sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/COVID-19 Lab confirmed & suspected - All Comparisons.txt")
print("   ")
print("   ")
print("   ")
print("Hydroxychloroquine vs. Standard care/Placebo")
print("   ")
print("   ")
print("Summary measure = OR")
print(estimates1)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates1.2)
print("   ")
print("   ")
print("Ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates2)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates2.2)
print("   ")
print("   ")
print("   ")
print("   ")
print("Iota-carrageenan, ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates3)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates3.2)
print("   ")
print("   ")
sink()

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/COVID-19 Lab confirmed & suspected - All Comparisons.pdf",width=8,height = 5)

### Hydroxychloroquine vs. Standard Care/Placebo, OR
yrange1 <- c(-7 - nrow(effsize1_labsinf), 1)
forest.default(effsize1_labsinf$yi, vi = effsize1_labsinf$vi, refline = 0,
               rows = seq(-2, -length(effsize1_labsinf$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize1_labsinf[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("COVID-19 Lab Confirmed & Suspected: Hydroxychloroquine vs. Standard Care/Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1))
  addpoly(estimates1[i,"mu"], ci.lb=estimates1[i,3], ci.ub=estimates1[i,4], atransf=exp,
          mlab=row.names(estimates1)[i], rows=yrange1[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Hydroxychloroquine vs. Standard Care/Placebo, RD
yrange1.2 <- c(-7 - nrow(effsize1_labsinf2), 1)
forest.default(effsize1_labsinf2$yi, vi = effsize1_labsinf2$vi, refline = 0,
               rows = seq(-2, -length(effsize1_labsinf2$yi) - 1, by = -1),width=0,
               xlim = c(-5,5),
               alim = c(-2,2),
               ylim = yrange1.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize1_labsinf2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
              digits=2)

title("COVID-19 Lab Confirmed & Suspected: Hydroxychloroquine vs. Standard Care/Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1.2))
  addpoly(estimates1.2[i,"mu"], ci.lb=estimates1.2[i,3], ci.ub=estimates1.2[i,4],
          mlab=row.names(estimates1.2)[i], rows=yrange1.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, OR
yrange2 <- c(-6 - nrow(effsize2_labsinf), 2)
forest.default(effsize2_labsinf$yi, vi = effsize2_labsinf$vi, refline = 0,
               rows = seq(-2, -length(effsize2_labsinf$yi) - 1, by = -1),width=0,
               xlim = c(-10,10),
               alim = c(-4, 1),
               ylim = yrange2, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize2_labsinf[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("COVID-19 Lab Confirmed & Suspected: Ivermectin vs. Placebo - OR", line=-1)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2))
  addpoly(estimates2[i,"mu"], ci.lb=estimates2[i,3], ci.ub=estimates2[i,4], atransf=exp,
          mlab=row.names(estimates2)[i], rows=yrange2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, RD
yrange2.2 <- c(-6 - nrow(effsize2_labsinf2), 2)
forest.default(effsize2_labsinf2$yi, vi = effsize2_labsinf2$vi, refline = 0,
               rows = seq(-2, -length(effsize2_labsinf2$yi) - 1, by = -1),width=0,
               xlim = c(-6, 6),
               alim = c(-1, 1),
               ylim = yrange2.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize2_labsinf2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("COVID-19 Lab Confirmed & Suspected: Ivermectin vs. Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2.2))
  addpoly(estimates2.2[i,"mu"], ci.lb=estimates2.2[i,3], ci.ub=estimates2.2[i,4],
          mlab=row.names(estimates2.2)[i], rows=yrange2.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Iota-carrageenan, ivermectin  vs. Standard Care/Placebo, OR
yrange3 <- c(-6 - nrow(effsize3), 2)
forest.default(effsize3$yi, vi = effsize3$vi, refline = 0,
               rows = seq(-2, -length(effsize3$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange3, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize3[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("COVID-19 Lab Confirmed & Suspected: Iota-carrageenan, ivermectin vs. Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates3))
  addpoly(estimates3[i,"mu"], ci.lb=estimates3[i,3], ci.ub=estimates3[i,4], atransf=exp,
          mlab=row.names(estimates3)[i], rows=yrange3[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Iota-carrageenan, ivermectin vs. Standard Care/Placebo, RD
yrange3.2 <- c(-6 - nrow(effsize3.2), 1)
forest.default(effsize3.2$yi, vi = effsize3.2$vi, refline = 0,
               rows = seq(-2, -length(effsize3.2$yi) - 1, by = -1),width=0,
               xlim = c(-10,10), ylim = yrange3.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize3.2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("COVID-19 Lab Confirmed & Suspected: Iota-carrageenan, ivermectin vs. Standard Care/Placebo", line=-1)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates3.2))
  addpoly(estimates3.2[i,"mu"], ci.lb=estimates3.2[i,3], ci.ub=estimates3.2[i,4],
          mlab=row.names(estimates3.2)[i], rows=yrange3.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

dev.off()

#############################################################################################
#### Admission to hospital
#############################################################################################

# load data:
data <- read.csv("admission to hospital - wide data format.csv",
                 header = TRUE)

#### Hydroxychloroquine vs. Standard care/Placebo
####################################################

data_hcq <- filter(data,
	t1 == "hydroxychloroquine" &
	t2 == "a_placebo/standard care")

effsize1_hospadmin <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize1_hospadmin2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_hcq)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("resource use / hospital stay / process", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner1 <- bayesmeta(effsize1_hospadmin, tau.prior=TP$dprior)
bm.hnorm110 <- bayesmeta(effsize1_hospadmin, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105 <- bayesmeta(effsize1_hospadmin, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner1.2 <- bayesmeta(effsize1_hospadmin2, tau.prior=TP$dprior)
bm.hnorm110.2 <- bayesmeta(effsize1_hospadmin2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105.2 <- bayesmeta(effsize1_hospadmin2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed1 <- rma.uni(effsize1_hospadmin, method="FE")
rma.random.DL1 <- rma.uni(effsize1_hospadmin, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed1.2 <- rma.uni(effsize1_hospadmin2, method="FE")
rma.random.DL1.2 <- rma.uni(effsize1_hospadmin2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates1 <- rbind("HNorm(0.5)" = c(bm.hnorm105$summary[2,1:2],   bm.hnorm105$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1$summary[2,1:2],bm.Turner1$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1$tau2.fix),rma.fixed1$b,rma.fixed1$ci.lb,rma.fixed1$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1$tau2), rma.random.DL1$b, rma.random.DL1$ci.lb,rma.random.DL1$ci.ub))

# estimates where RD is the measure of effects
estimates1.2 <- rbind("HNorm(0.5)" = c(bm.hnorm105.2$summary[2,1:2],   bm.hnorm105.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1.2$summary[2,1:2],bm.Turner1.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1.2$tau2.fix),rma.fixed1.2$b,rma.fixed1$ci.lb,rma.fixed1.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1.2$tau2), rma.random.DL1.2$b, rma.random.DL1.2$ci.lb,rma.random.DL1.2$ci.ub))


#### Ivermectin vs. Standard care/Placebo
####################################################

data_iv <- filter(data,
	t1 == "ivermectin" &
	t2 == "a_placebo/standard care")

effsize2_hospadmin <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_iv)

effsize2_hospadmin2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_iv)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("resource use / hospital stay / process", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner2 <- bayesmeta(effsize2_hospadmin, tau.prior=TP$dprior)
bm.hnorm210 <- bayesmeta(effsize2_hospadmin, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205 <- bayesmeta(effsize2_hospadmin, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner2.2 <- bayesmeta(effsize2_hospadmin2, tau.prior=TP$dprior)
bm.hnorm210.2 <- bayesmeta(effsize2_hospadmin2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205.2 <- bayesmeta(effsize2_hospadmin2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed2 <- rma.uni(effsize1, method="FE")
rma.random.DL2 <- rma.uni(effsize1, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed2.2 <- rma.uni(effsize2_hospadmin2, method="FE")
rma.random.DL2.2 <- rma.uni(effsize2_hospadmin2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates2 <- rbind("HNorm(0.5)" = c(bm.hnorm205$summary[2,1:2],   bm.hnorm205$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2$summary[2,1:2],bm.Turner2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2$tau2.fix),rma.fixed2$b,rma.fixed2$ci.lb,rma.fixed2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2$tau2), rma.random.DL2$b, rma.random.DL2$ci.lb,rma.random.DL2$ci.ub))

# estimates where RD is the measure of effects
estimates2.2 <- rbind("HNorm(0.5)" = c(bm.hnorm205.2$summary[2,1:2],   bm.hnorm205.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2.2$summary[2,1:2],bm.Turner2.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2.2$tau2.fix),rma.fixed2.2$b,rma.fixed2$ci.lb,rma.fixed2.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2.2$tau2), rma.random.DL2.2$b, rma.random.DL2.2$ci.lb,rma.random.DL2.2$ci.ub))


#### Iota-carrageenan, ivermectin vs. Standard care/Placebo
#############################################################

data_iciv <- filter(data,
	t1 == "iota-carrageenan, ivermectin" &
	t2 == "a_placebo/standard care")

effsize3_hospadmin <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize3_hospadmin2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_hcq)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("resource use / hospital stay / process", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner3 <- bayesmeta(effsize3_hospadmin, tau.prior=TP$dprior)
bm.hnorm310 <- bayesmeta(effsize3_hospadmin, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm305 <- bayesmeta(effsize3_hospadmin, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner3.2 <- bayesmeta(effsize3_hospadmin2, tau.prior=TP$dprior)
bm.hnorm310.2 <- bayesmeta(effsize3_hospadmin2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm305.2 <- bayesmeta(effsize3_hospadmin2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed3 <- rma.uni(effsize3_hospadmin, method="FE")
rma.random.DL3 <- rma.uni(effsize3_hospadmin, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed3.2 <- rma.uni(effsize3_hospadmin2, method="FE")
rma.random.DL3.2 <- rma.uni(effsize3_hospadmin2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates3 <- rbind("HNorm(0.5)" = c(bm.hnorm305$summary[2,1:2],   bm.hnorm305$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner3$summary[2,1:2],bm.Turner3$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed3$tau2.fix),rma.fixed3$b,rma.fixed3$ci.lb,rma.fixed3$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL3$tau2), rma.random.DL3$b, rma.random.DL3$ci.lb,rma.random.DL3$ci.ub))

# estimates where RD is the measure of effects
estimates3.2 <- rbind("HNorm(0.5)" = c(bm.hnorm305.2$summary[2,1:2],   bm.hnorm305.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner3.2$summary[2,1:2],bm.Turner3.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed3.2$tau2.fix),rma.fixed3.2$b,rma.fixed3$ci.lb,rma.fixed3.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL3.2$tau2), rma.random.DL3.2$b, rma.random.DL3.2$ci.lb,rma.random.DL3.2$ci.ub))


#### Results & Forest plot
####################################################

sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/Admission to Hospital - All Comparisons.txt")
print("   ")
print("   ")
print("   ")
print("Hydroxychloroquine vs. Standard care/Placebo")
print("   ")
print("   ")
print("Summary measure = OR")
print(estimates1)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates1.2)
print("   ")
print("   ")
print("Ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates2)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates2.2)
print("   ")
print("   ")
print("   ")
print("   ")
print("Iota-carrageenan, ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates3)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates3.2)
print("   ")
print("   ")
sink()

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/Admission to Hospital - All Comparisons.pdf",width=8,height = 5)

### Hydroxychloroquine vs. Standard Care/Placebo, OR
yrange1 <- c(-7 - nrow(effsize1_hospadmin), 1)
forest.default(effsize1_hospadmin$yi, vi = effsize1_hospadmin$vi, refline = 0,
               rows = seq(-2, -length(effsize1_hospadmin$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize1_hospadmin[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("Admission to Hospital: Hydroxychloroquine vs. Standard Care/Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1))
  addpoly(estimates1[i,"mu"], ci.lb=estimates1[i,3], ci.ub=estimates1[i,4], atransf=exp,
          mlab=row.names(estimates1)[i], rows=yrange1[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Hydroxychloroquine vs. Standard Care/Placebo, RD
yrange1.2 <- c(-7 - nrow(effsize1_hospadmin2), 1)
forest.default(effsize1_hospadmin2$yi, vi = effsize1_hospadmin2$vi, refline = 0,
               rows = seq(-2, -length(effsize1_hospadmin2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize1_hospadmin2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("Admission to Hospital: Hydroxychloroquine vs. Standard Care/Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1.2))
  addpoly(estimates1.2[i,"mu"], ci.lb=estimates1.2[i,3], ci.ub=estimates1.2[i,4],
          mlab=row.names(estimates1.2)[i], rows=yrange1.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, OR
yrange2 <- c(-6 - nrow(effsize2), 2)
forest.default(effsize2$yi, vi = effsize2$vi, refline = 0,
               rows = seq(-2, -length(effsize2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange2, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("Admission to Hospital: Ivermectin vs. Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2))
  addpoly(estimates2[i,"mu"], ci.lb=estimates2[i,3], ci.ub=estimates2[i,4], atransf=exp,
          mlab=row.names(estimates2)[i], rows=yrange2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, RD
yrange2.2 <- c(-6 - nrow(effsize2.2), 2)
forest.default(effsize2.2$yi, vi = effsize2.2$vi, refline = 0,
               rows = seq(-2, -length(effsize2.2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10), ylim = yrange2.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize2.2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("Admission to Hospital: Ivermectin vs. Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2.2))
  addpoly(estimates2.2[i,"mu"], ci.lb=estimates2.2[i,3], ci.ub=estimates2.2[i,4],
          mlab=row.names(estimates2.2)[i], rows=yrange2.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Iota-carrageenan, ivermectin  vs. Standard Care/Placebo, OR
yrange3 <- c(-6 - nrow(effsize3), 2)
forest.default(effsize3$yi, vi = effsize3$vi, refline = 0,
               rows = seq(-2, -length(effsize3$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange3, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize3[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("Admission to Hospital: Iota-carrageenan, ivermectin vs. Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates3))
  addpoly(estimates3[i,"mu"], ci.lb=estimates3[i,3], ci.ub=estimates3[i,4], atransf=exp,
          mlab=row.names(estimates3)[i], rows=yrange3[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Iota-carrageenan, ivermectin vs. Standard Care/Placebo, RD
yrange3.2 <- c(-6 - nrow(effsize3.2), 1)
forest.default(effsize3.2$yi, vi = effsize3.2$vi, refline = 0,
               rows = seq(-2, -length(effsize3.2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10), ylim = yrange3.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize3.2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("Admission to Hospital: Iota-carrageenan, ivermectin vs. Standard Care/Placebo", line=-1)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates3.2))
  addpoly(estimates3.2[i,"mu"], ci.lb=estimates3.2[i,3], ci.ub=estimates3.2[i,4], atransf=exp,
          mlab=row.names(estimates3.2)[i], rows=yrange3.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

dev.off()


#############################################################################################
#### Adverse events leading to discontinuation
#############################################################################################

# load data:
data <- read.csv("AEs - wide data format.csv",
                 header = TRUE)

#### Hydroxychloroquine vs. Standard care/Placebo
####################################################

data_hcq <- filter(data,
	t1 == "hydroxychloroquine" &
	t2 == "a_placebo/standard care")

effsize1_AE <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize1_AE2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_hcq)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("adverse events", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner1 <- bayesmeta(effsize1_AE, tau.prior=TP$dprior)
bm.hnorm110 <- bayesmeta(effsize1_AE, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105 <- bayesmeta(effsize1_AE, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner1.2 <- bayesmeta(effsize1_AE2, tau.prior=TP$dprior)
bm.hnorm110.2 <- bayesmeta(effsize1_AE2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105.2 <- bayesmeta(effsize1_AE2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed1 <- rma.uni(effsize1_AE, method="FE")
rma.random.DL1 <- rma.uni(effsize1_AE, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed1.2 <- rma.uni(effsize1_AE2, method="FE")
rma.random.DL1.2 <- rma.uni(effsize1_AE2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates1 <- rbind("HNorm(0.5)" = c(bm.hnorm105$summary[2,1:2],   bm.hnorm105$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1$summary[2,1:2],bm.Turner1$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1$tau2.fix),rma.fixed1$b,rma.fixed1$ci.lb,rma.fixed1$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1$tau2), rma.random.DL1$b, rma.random.DL1$ci.lb,rma.random.DL1$ci.ub))

# estimates where RD is the measure of effects
estimates1.2 <- rbind("HNorm(0.5)" = c(bm.hnorm105.2$summary[2,1:2],   bm.hnorm105.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner1.2$summary[2,1:2],bm.Turner1.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed1.2$tau2.fix),rma.fixed1.2$b,rma.fixed1$ci.lb,rma.fixed1.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL1.2$tau2), rma.random.DL1.2$b, rma.random.DL1.2$ci.lb,rma.random.DL1.2$ci.ub))


#### Ivermectin vs. Standard care/Placebo
####################################################

data_iv <- filter(data,
	t1 == "ivermectin" &
	t2 == "a_placebo/standard care")

effsize2_AE <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_iv)

effsize2_AE2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_iv)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("adverse events", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner2 <- bayesmeta(effsize2_AE, tau.prior=TP$dprior)
bm.hnorm210 <- bayesmeta(effsize2_AE, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205 <- bayesmeta(effsize2_AE, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner2.2 <- bayesmeta(effsize2_AE2, tau.prior=TP$dprior)
bm.hnorm210.2 <- bayesmeta(effsize2_AE2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm205.2 <- bayesmeta(effsize2_AE2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed2 <- rma.uni(effsize1, method="FE")
rma.random.DL2 <- rma.uni(effsize1, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed2.2 <- rma.uni(effsize2_AE2, method="FE")
rma.random.DL2.2 <- rma.uni(effsize2_AE2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates2 <- rbind("HNorm(0.5)" = c(bm.hnorm205$summary[2,1:2],   bm.hnorm205$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2$summary[2,1:2],bm.Turner2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2$tau2.fix),rma.fixed2$b,rma.fixed2$ci.lb,rma.fixed2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2$tau2), rma.random.DL2$b, rma.random.DL2$ci.lb,rma.random.DL2$ci.ub))

# estimates where RD is the measure of effects
estimates2.2 <- rbind("HNorm(0.5)" = c(bm.hnorm205.2$summary[2,1:2],   bm.hnorm205.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner2.2$summary[2,1:2],bm.Turner2.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed2.2$tau2.fix),rma.fixed2.2$b,rma.fixed2$ci.lb,rma.fixed2.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL2.2$tau2), rma.random.DL2.2$b, rma.random.DL2.2$ci.lb,rma.random.DL2.2$ci.ub))


#### Iota-carrageenan, ivermectin vs. Standard care/Placebo
#############################################################

data_iciv <- filter(data,
	t1 == "iota-carrageenan, ivermectin" &
	t2 == "a_placebo/standard care")

effsize3_AE <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize3_AE2 <- escalc(measure = "RD",
                        ai = e.events,  n1i = e.total,
                        ci = c.events, n2i = c.total,
                        slab = study,
                        data = data_hcq)

# determine corresponding prior parameters(?TurnerEtAlPrior to help):
TP <- TurnerEtAlPrior("adverse events", "pharma", "placebo / control")
print(TP)

# a prior 95 percent interval for tau:
TP$qprior(c(0.025,0.975))

# perform bayesian meta analysis with TurnerEtAl Prior and half-normal:
bm.Turner3 <- bayesmeta(effsize3_AE, tau.prior=TP$dprior)
bm.hnorm310 <- bayesmeta(effsize3_AE, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm305 <- bayesmeta(effsize3_AE, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner3.2 <- bayesmeta(effsize3_AE2, tau.prior=TP$dprior)
bm.hnorm310.2 <- bayesmeta(effsize3_AE2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm305.2 <- bayesmeta(effsize3_AE2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed3 <- rma.uni(effsize3_AE, method="FE")
rma.random.DL3 <- rma.uni(effsize3_AE, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed3.2 <- rma.uni(effsize3_AE2, method="FE")
rma.random.DL3.2 <- rma.uni(effsize3_AE2, method="DL")

# assemble estimates (on log-OR scale);
# estimates for tau and mu (posterior medians for Bayesian analyses)
# and 95% interval for mu:
estimates3 <- rbind("HNorm(0.5)" = c(bm.hnorm305$summary[2,1:2],   bm.hnorm305$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner3$summary[2,1:2],bm.Turner3$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed3$tau2.fix),rma.fixed3$b,rma.fixed3$ci.lb,rma.fixed3$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL3$tau2), rma.random.DL3$b, rma.random.DL3$ci.lb,rma.random.DL3$ci.ub))

# estimates where RD is the measure of effects
estimates3.2 <- rbind("HNorm(0.5)" = c(bm.hnorm305.2$summary[2,1:2],   bm.hnorm305.2$summary[5:6,2]),
                    "Turner Prior"=c(bm.Turner3.2$summary[2,1:2],bm.Turner3.2$summary[5:6,2]),
                    "Frequentist.fixed" = c(sqrt(rma.fixed3.2$tau2.fix),rma.fixed3.2$b,rma.fixed3$ci.lb,rma.fixed3.2$ci.ub),
                    "Frequentist.random.DL" = c(sqrt(rma.random.DL3.2$tau2), rma.random.DL3.2$b, rma.random.DL3.2$ci.lb,rma.random.DL3.2$ci.ub))


#### Results & Forest plot
####################################################

sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/AEs leading to Discontinuation - All Comparisons.txt")
print("   ")
print("   ")
print("   ")
print("Hydroxychloroquine vs. Standard care/Placebo")
print("   ")
print("   ")
print("Summary measure = OR")
print(estimates1)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates1.2)
print("   ")
print("   ")
print("Ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates2)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates2.2)
print("   ")
print("   ")
print("   ")
print("   ")
print("Iota-carrageenan, ivermectin vs. Standard care/Placebo")
print("Summary measure = OR")
print(estimates3)
print("   ")
print("   ")
print("Summary measure = RD")
print(estimates3.2)
print("   ")
print("   ")
sink()

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/7. Pairwise MA/AEs Leading to Discontinuation - All Comparisons.pdf",width=8,height = 5)

### Hydroxychloroquine vs. Standard Care/Placebo, OR
yrange1 <- c(-7 - nrow(effsize1_AE), 1)
forest.default(effsize1_AE$yi, vi = effsize1_AE$vi, refline = 0,
               rows = seq(-2, -length(effsize1_AE$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize1_AE[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("AEs Leading to Discontinuation: Hydroxychloroquine vs. Standard Care/Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1))
  addpoly(estimates1[i,"mu"], ci.lb=estimates1[i,3], ci.ub=estimates1[i,4], atransf=exp,
          mlab=row.names(estimates1)[i], rows=yrange1[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Hydroxychloroquine vs. Standard Care/Placebo, RD
yrange1.2 <- c(-7 - nrow(effsize1_AE2), 1)
forest.default(effsize1_AE2$yi, vi = effsize1_AE2$vi, refline = 0,
               rows = seq(-2, -length(effsize1_AE2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize1_AE2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("AEs Leading to Discontinuation: Hydroxychloroquine vs. Standard Care/Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1.2))
  addpoly(estimates1.2[i,"mu"], ci.lb=estimates1.2[i,3], ci.ub=estimates1.2[i,4],
          mlab=row.names(estimates1.2)[i], rows=yrange1.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, OR
yrange2 <- c(-6 - nrow(effsize2), 2)
forest.default(effsize2$yi, vi = effsize2$vi, refline = 0,
               rows = seq(-2, -length(effsize2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange2, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("AEs Leading to Discontinuation: Ivermectin vs. Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2))
  addpoly(estimates2[i,"mu"], ci.lb=estimates2[i,3], ci.ub=estimates2[i,4], atransf=exp,
          mlab=row.names(estimates2)[i], rows=yrange2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Ivermectin vs. Standard Care/Placebo, RD
yrange2.2 <- c(-6 - nrow(effsize2.2), 2)
forest.default(effsize2.2$yi, vi = effsize2.2$vi, refline = 0,
               rows = seq(-2, -length(effsize2.2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange2.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize2.2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("AEs Leading to Discontinuation: Ivermectin vs. Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates2.2))
  addpoly(estimates2.2[i,"mu"], ci.lb=estimates2.2[i,3], ci.ub=estimates2.2[i,4],
          mlab=row.names(estimates2.2)[i], rows=yrange2.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Iota-carrageenan, ivermectin  vs. Standard Care/Placebo, OR
yrange3 <- c(-6 - nrow(effsize3), 2)
forest.default(effsize3$yi, vi = effsize3$vi, refline = 0,
               rows = seq(-2, -length(effsize3$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange3, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize3[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("AEs Leading to Discontinuation: Iota-carrageenan, ivermectin vs. Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates3))
  addpoly(estimates3[i,"mu"], ci.lb=estimates3[i,3], ci.ub=estimates3[i,4], atransf=exp,
          mlab=row.names(estimates3)[i], rows=yrange3[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Iota-carrageenan, ivermectin vs. Standard Care/Placebo, RD
yrange3.2 <- c(-6 - nrow(effsize3.2), 1)
forest.default(effsize3.2$yi, vi = effsize3.2$vi, refline = 0,
               rows = seq(-2, -length(effsize3.2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange3.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize3.2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("AEs Leading to Discontinuation: Iota-carrageenan, ivermectin vs. Standard Care/Placebo", line=-1)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates3.2))
  addpoly(estimates3.2[i,"mu"], ci.lb=estimates3.2[i,3], ci.ub=estimates3.2[i,4],
          mlab=row.names(estimates3.2)[i], rows=yrange3.2[1]+5-i, col=colvec[i],cex=1.5,width =0)

dev.off()
