#############################################################################################
### Bayesian NMA in gemtc
### Original code from LG, with modifications by AQ to faciliate analysis.

### Jan 13, 2021
#############################################################################################

# load necessary libraries:
require("metafor")
require("bayesmeta")
library(dplyr)
library(forestplot)
library(meta)
library(scales)
library(gemtc)

setwd("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Data")

#############################################################################################
## MORTALITY

data <- read.csv("mortality - gemtc.csv",
                 header = TRUE)

network <- mtc.network(data)

# network plot
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/2. Network plot/mortality.pdf",
  width = 8,
  height = 3)

M=data$node
plot(network,layout=igraph::layout.circle, dynamic.edge.width=T, margin=0,
     edge.color="black",vertex.color="red",vertex.size=M,
     vertex.shape="circle",vertex.label.dist=-0,
     vertex.label.cex=1.5,vertex.label.color="blue",vertex.label.degree=pi/1,
     use.description=FALSE)
dev.off()

### set informative prior based on pairwise MA
model <- mtc.model(network,
                    type = "consistency",
                    likelihood = "binom",
                    link="logit",
                    linearModel = "random",
                    n.chain = 3,
                    powerAdjust=NA,
                    dic = TRUE,
                    hy.prior=mtc.hy.prior("var", "dlnorm", -3.95, 0.5569))

results <- mtc.run(model, sampler ="JAGS", n.adapt=50000, n.iter=2000000, thin=1)

### save network summary, posterior summary, convergence diagnostics
sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/3. Network summary, Posterior summary, convergence diagnostic/mortality.csv")
summary(network)
summary(results)
gelman.diag(results)
sink()

### Creates a forest plot of the relative effects
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/4. NMA forest plot/mortality.pdf",
    width = 8,
    height = 3)

forest(relative.effect(results,"a_StandardCarePlacebo"),
    digits = 4,
    use.description = TRUE)
dev.off()

###generate league table
mtcresults = as.data.frame(round(exp(relative.effect.table(results)),2))
write.csv(mtcresults, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/5. League table/mortality.csv")

### calculating probability of rank, preferredDirection=1 for benefit outcomes
rank.prob <- rank.probability(results, preferredDirection=-1)

### robability for each treatment to be best, second best, etc.
print(rank.prob)
write.csv(rank.prob, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/Prbest_mortality.csv")

##Generate quantile ranks
rank.quantiles <- rank.quantiles(rank.prob, probs=c("2.5%"=0.025, "50%"=0.5, "97.5%"=0.975))
write.csv(rank.quantiles, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/Rank.quantiles_mortality.csv")

# plot a 'rankogram'
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/rankogram_mortality.pdf")
plot(rank.prob, beside=TRUE)
dev.off()

### calculating SUCRA
cumrank.prob <- apply(t(rank.prob), 2, cumsum)
sucra <- round(colMeans(cumrank.prob[-nrow(cumrank.prob),]),4)
write.csv(sucra, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/SUCRA_mortality.csv")

####Nodesplit analysis
result.node <- mtc.nodesplit(network,
                             likelihood="binom",
                             link="logit",
                             linearModel="random",
                             n.chain =3,
                             n.adapt=50000, n.iter=2000000, thin=1,
                             hy.prior=mtc.hy.prior("var", "dlnorm", -3.95,0.5569))

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/8. Nodesplit plots/mortality.pdf")
plot(summary(result.node),xlim=log(c(0.1,5)),digits=4)
dev.off()

#############################################################################################
## COVID-19 LAB CONFIRMED

data <- read.csv("COVID-19 lab confirmed - gemtc.csv",
                 header = TRUE)

network <- mtc.network(data)

# network plot
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/2. Network plot/COVID-19 lab confirmed.pdf",
  width = 8,
  height = 3)

M=data$node
plot(network,layout=igraph::layout.circle, dynamic.edge.width=T, margin=0,
     edge.color="black",vertex.color="red",vertex.size=M,
     vertex.shape="circle",vertex.label.dist=-0,
     vertex.label.cex=1.5,vertex.label.color="blue",vertex.label.degree=pi/1,
     use.description=FALSE)
dev.off()

### set informative prior based on pairwise MA
model <- mtc.model(network,
                    type = "consistency",
                    likelihood = "binom",
                    link="logit",
                    linearModel = "random",
                    n.chain = 3,
                    powerAdjust=NA,
                    dic = TRUE,
                    hy.prior=mtc.hy.prior("var", "dlnorm", -2.49, 0.4328))

results <- mtc.run(model, sampler ="JAGS", n.adapt=50000, n.iter=2000000, thin=1)

### save network summary, posterior summary, convergence diagnostics
sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/3. Network summary, Posterior summary, convergence diagnostic/COVID-19 lab confirmed.csv")
summary(network)
summary(results)
gelman.diag(results)
sink()

### Creates a forest plot of the relative effects
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/4. NMA forest plot/COVID-19 lab confirmed.pdf",
    width = 8,
    height = 3)

forest(relative.effect(results,"a_SCPlacebo"),
    digits = 4,
    use.description = TRUE)
dev.off()

###generate league table
mtcresults = as.data.frame(round(exp(relative.effect.table(results)),2))
write.csv(mtcresults, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/5. League table/COVID-19 lab confirmed.csv")

### calculating probability of rank, preferredDirection=1 for benefit outcomes
rank.prob <- rank.probability(results, preferredDirection=-1)

### robability for each treatment to be best, second best, etc.
print(rank.prob)
write.csv(rank.prob, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/Prbest_COVID-19 lab confirmed.csv")

##Generate quantile ranks
rank.quantiles <- rank.quantiles(rank.prob, probs=c("2.5%"=0.025, "50%"=0.5, "97.5%"=0.975))
write.csv(rank.quantiles, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/Rank.quantiles_COVID-19 lab confirmed.csv")

# plot a 'rankogram'
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/rankogram_COVID-19 lab confirmed.pdf")
plot(rank.prob, beside=TRUE)
dev.off()

### calculating SUCRA
cumrank.prob <- apply(t(rank.prob), 2, cumsum)
sucra <- round(colMeans(cumrank.prob[-nrow(cumrank.prob),]),4)
write.csv(sucra, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/SUCRA_COVID-19 lab confirmed.csv")

####Nodesplit analysis
result.node <- mtc.nodesplit(network,
                             likelihood="binom",
                             link="logit",
                             linearModel="random",
                             n.chain =3,
                             n.adapt=50000, n.iter=2000000, thin=1,
                             hy.prior=mtc.hy.prior("var", "dlnorm", -2.49, 0.4328))

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/8. Nodesplit plots/COVID-19 lab confirmed.pdf")
plot(summary(result.node),xlim=log(c(0.1,5)),digits=4)
dev.off()


#############################################################################################
## COVID-19 LAB CONFIRMED & SUSPECTED

data <- read.csv("COVID-19 lab confirmed & suspected - gemtc.csv",
                 header = TRUE)

## network created using RD
network <- mtc.network(data)

# network plot
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/2. Network plot/COVID-19 lab confirmed & suspected.pdf",
    width = 8,
    height = 3)

M=data$node
plot(network,layout=igraph::layout.circle, dynamic.edge.width=T, margin=0,
     edge.color="black",vertex.color="red",vertex.size=M,
     vertex.shape="circle",vertex.label.dist=-0,
     vertex.label.cex=1.5,vertex.label.color="blue",vertex.label.degree=pi/1,
     use.description=FALSE)
dev.off()

### set informative prior based on pairwise MA
model <- mtc.model(network,
                   type = "consistency",
                   likelihood = "binom",
                   link="logit",
                   linearModel = "random",
                   n.chain = 3,
                   powerAdjust=NA,
                   dic = TRUE,
                   hy.prior=mtc.hy.prior("var", "dlnorm", -2.49, 0.4328))

results <- mtc.run(model, sampler ="JAGS", n.adapt=50000, n.iter=2000000, thin=1)

### save network summary, posterior summary, convergence diagnostics
sink("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/3. Network summary, Posterior summary, convergence diagnostic/COVID-19 lab confirmed & suspected.csv")
summary(network)
summary(results)
gelman.diag(results)
sink()

### Creates a forest plot of the relative effects
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/4. NMA forest plot/COVID-19 lab confirmed & suspected.pdf",
    width = 8,
    height = 3)

forest(relative.effect(results,"a_placebo_standardcare"),
       digits = 4,
       use.description = TRUE)
dev.off()

###generate league table
mtcresults = as.data.frame(round(exp(relative.effect.table(results)),4))
write.csv(mtcresults, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/5. League table/COVID-19 lab confirmed & suspected.csv")

### calculating probability of rank, preferredDirection=1 for benefit outcomes
rank.prob <- rank.probability(results, preferredDirection=-1)

### probability for each treatment to be best, second best, etc.
print(rank.prob)
write.csv(rank.prob, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/Prbest_COVID-19 lab confirmed & suspected.csv")

##Generate quantile ranks
rank.quantiles <- rank.quantiles(rank.prob, probs=c("2.5%"=0.025, "50%"=0.5, "97.5%"=0.975))
write.csv(rank.quantiles, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/Rank.quantiles_COVID-19 lab confirmed & suspected.csv")

# plot a 'rankogram'
pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/rankogram_COVID-19 lab confirmed & suspected.pdf")
plot(rank.prob, beside=TRUE)
dev.off()

### calculating SUCRA
cumrank.prob <- apply(t(rank.prob), 2, cumsum)
sucra <- round(colMeans(cumrank.prob[-nrow(cumrank.prob),]),4)
write.csv(sucra, file="/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/6. Rank/SUCRA_COVID-19 lab confirmed & suspected.csv")

####Nodesplit analysis
result.node <- mtc.nodesplit(network,
                             likelihood="normal",
                             link="binom",
                             linearModel="clogclog",
                             n.chain =3,
                             n.adapt=50000, n.iter=2000000, thin=1,
                             hy.prior=mtc.hy.prior("var", "dlnorm", -2.49, 0.4328))

pdf("/Users/aqasim/Dropbox/COVID-19 LNMA/Prophylaxis Analysis/January 13, 2021/8. Nodesplit plots/COVID-19 lab confirmed & suspected.pdf")
plot(summary(result.node2),xlim=log(c(0.1,5)),digits=4)
dev.off()
