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

setwd("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Subgroups Data/Adverse event data")


#############################################################################################
#### Cardiac toxicity
#############################################################################################

# load data:
data <- read.csv("cardiac toxicity_wide.csv",
                 header = TRUE)

#### Hydroxychloroquine vs. Standard care/Placebo
####################################################

data_hcq <- filter(data,
	t1 == "hydroxychloroquine" &
	t2 == "a_SCPlacebo")

effsize1_ct <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize1_ct2 <- escalc(measure = "RD",
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
bm.Turner1 <- bayesmeta(effsize1_ct, tau.prior=TP$dprior)
bm.hnorm110 <- bayesmeta(effsize1_ct, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105 <- bayesmeta(effsize1_ct, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner1.2 <- bayesmeta(effsize1_ct2, tau.prior=TP$dprior)
bm.hnorm110.2 <- bayesmeta(effsize1_ct2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105.2 <- bayesmeta(effsize1_ct2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed1 <- rma.uni(effsize1_ct, method="FE")
rma.random.DL1 <- rma.uni(effsize1_ct, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed1.2 <- rma.uni(effsize1_ct2, method="FE")
rma.random.DL1.2 <- rma.uni(effsize1_ct2, method="DL")

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


#### Results & Forest plot
####################################################

sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Subgroups Analysis/Cardiac toxicity - Hydroxychloroquine.txt")
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
sink()

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Subgroups Analysis/Cardiac toxicity - Hydroxychloroquine.pdf",width=8,height = 5)

### Hydroxychloroquine vs. Standard Care/Placebo, OR
yrange1 <- c(-7 - nrow(effsize1_ct), 1)
forest.default(effsize1_ct$yi, vi = effsize1_ct$vi, refline = 0,
               rows = seq(-2, -length(effsize1_ct$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize1_ct[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("Cardiac Toxicity: Hydroxychloroquine vs. Standard Care/Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1))
  addpoly(estimates1[i,"mu"], ci.lb=estimates1[i,3], ci.ub=estimates1[i,4], atransf=exp,
          mlab=row.names(estimates1)[i], rows=yrange1[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Hydroxychloroquine vs. Standard Care/Placebo, RD
yrange1.2 <- c(-7 - nrow(effsize1_ct2), 1)
forest.default(effsize1_ct2$yi, vi = effsize1_ct2$vi, refline = 0,
               rows = seq(-2, -length(effsize1_ct2$yi) - 1, by = -1),width=0,
               xlim = c(-5, 5),
               alim = c(-2, 2),
               ylim = yrange1.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize1_ct2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("Cardiac toxicity: Hydroxychloroquine vs. Standard Care/Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1.2))
  addpoly(estimates1.2[i,"mu"], ci.lb=estimates1.2[i,3], ci.ub=estimates1.2[i,4],
          mlab=row.names(estimates1.2)[i], rows=yrange1.2[1]+5-i, col=colvec[i],cex=1.5,width =0)


#############################################################################################
#### Non-serious GI effects
#############################################################################################

# load data:
data <- read.csv("GI-related AEs_wide.csv",
                 header = TRUE)

#### Hydroxychloroquine vs. Standard care/Placebo
####################################################

data_hcq <- filter(data,
	t1 == "hydroxychloroquine" &
	t2 == "a_SCPlacebo")

effsize1_gi <- escalc(measure = "OR",
                    ai = e.events,  n1i = e.total,
                    ci = c.events, n2i = c.total,
                    slab = study,
                    data = data_hcq)

effsize1_gi2 <- escalc(measure = "RD",
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
bm.Turner1 <- bayesmeta(effsize1_gi, tau.prior=TP$dprior)
bm.hnorm110 <- bayesmeta(effsize1_gi, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105 <- bayesmeta(effsize1_gi, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform MA with RD as the measure of effect
bm.Turner1.2 <- bayesmeta(effsize1_gi2, tau.prior=TP$dprior)
bm.hnorm110.2 <- bayesmeta(effsize1_gi2, tau.prior = function(x){dhalfnormal(x,scale=1.0)})
bm.hnorm105.2 <- bayesmeta(effsize1_gi2, tau.prior = function(x){dhalfnormal(x,scale=0.5)})

# perform FE MA
rma.fixed1 <- rma.uni(effsize1_gi, method="FE")
rma.random.DL1 <- rma.uni(effsize1_gi, method="DL")

# perform FE MA with RD as the measure of effect
rma.fixed1.2 <- rma.uni(effsize1_gi2, method="FE")
rma.random.DL1.2 <- rma.uni(effsize1_gi2, method="DL")

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


#### Results & Forest plot
####################################################

sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Subgroups Analysis/GI-related AEs - Hydroxychloroquine.txt")
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
sink()

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Subgroups Analysis/GI-related AEs - Hydroxychloroquine.pdf",width=8,height = 5)

### Hydroxychloroquine vs. Standard Care/Placebo, OR
yrange1 <- c(-7 - nrow(effsize1_gi), 1)
forest.default(effsize1_gi$yi, vi = effsize1_gi$vi, refline = 0,
               rows = seq(-2, -length(effsize1_gi$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1, top=2, steps=5, level=95,
               xlab="Odds ratio", slab = effsize1_gi[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               atransf=exp, digits=2)

title("GI-related adverse effects: Hydroxychloroquine vs. Standard Care/Placebo - OR", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1))
  addpoly(estimates1[i,"mu"], ci.lb=estimates1[i,3], ci.ub=estimates1[i,4], atransf=exp,
          mlab=row.names(estimates1)[i], rows=yrange1[1]+5-i, col=colvec[i],cex=1.5,width =0)

### Hydroxychloroquine vs. Standard Care/Placebo, RD
yrange1.2 <- c(-7 - nrow(effsize1_gi2), 1)
forest.default(effsize1_gi2$yi, vi = effsize1_gi2$vi, refline = 0,
               rows = seq(-2, -length(effsize1_gi2$yi) - 1, by = -1),width=0,
               #xlim = c(-10,10),
               ylim = yrange1.2, top=2, steps=5, level=95,
               xlab="Risk difference", slab = effsize1_gi2[,"study"],efac=1, pch=15,cex=1.5,cex.lab=1.5,
               digits=2)

title("GI-related adverse effects: Hydroxychloroquine vs. Standard Care/Placebo - RD", line=0)
colvec <- c("green","blue","red","yellow")
par(mar=c(0.5,0.5,0.5,0.5))
par("mai")
for (i in 1:nrow(estimates1.2))
  addpoly(estimates1.2[i,"mu"], ci.lb=estimates1.2[i,3], ci.ub=estimates1.2[i,4],
          mlab=row.names(estimates1.2)[i], rows=yrange1.2[1]+5-i, col=colvec[i],cex=1.5,width =0)
