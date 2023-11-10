# 5/59
require(geoR)
data(elevation)
help(elevation)
points(elevation, cex.min=2,cex.max=2,col="gray")
summary(elevation)
plot(elevation, lowess=T)

# 6/59
library(akima)
int=interp(elevation$coords[,1],elevation$coords[,2],elevation$data)
image(int,xlim=range(elevation$coords[,1]),ylim=range(elevation$coords[,2]))
contour(int,add=T)

# 8/59
par(mfrow=c(1,3))
points(elevation, cex.max=2.5)
points(elevation,trend="1st",pt.div=2,abs=T,cex.max=2.5)
points(elevation,trend="2nd",pt.div=2,abs=T,cex.max=2.5)

# 11/59
plot(variog(elevation,uvec=seq(0,5,by=0.5)),type="b")
res1.v <- variog(elevation,trend="1st",uvec=seq(0,5,by=0.5))
plot(res1.v,type="b")

# 15/59
m10=likfit(elevation, ini=c(3000,2), cov.model="matern",kappa=1.5)
m10
# likfit: estimated model parameters:
#   beta tausq sigmasq phi
# " 848.316" "48.138" "3510.233" "1.198"
# Practical Range with cor=0.05 for
# asymptotic range: 5.685412
# likfit: maximised log-likelihood = -242.1

# 16/59
m11=likfit(elevation,trend="1st",ini=c(1300,2),cov.model="matern",kappa=1.5)
m11
# likfit: estimated model parameters:
#   beta0 beta1 beta2 tausq sigmasq phi
# " 912.48" "-4.99" "-16.46" "34.89" "1693.13" "0.81"
# Practical Range with cor=0.05 for asymptotic
# range: 3.824244
# likfit: maximised log-likelihood = -240.1

# 17/59
plot(variog(elevation,uvec=seq(0,5,by=0.5)),type="b")
lines.variomodel(seq(0,5,by=0.5),
   cov.pars=c(m10$sigmasq,m10$phi),
   cov.model="mat",kap=1.5,nug=m10$tausq)

# 22/59
locs=pred_grid(c(0,6.3),c(0,6.3),by=0.1)
KC=krige.control(type="sk",obj.mod=m10)
skt=krige.conv(elevation, krige=KC,loc=locs)
# spatial maps of predictions
pred.lim=range(c(sk$pred,skt$pred))
sd.lim=range(c(sk$kr,skt$kr))
image(sk,col=gray(seq(1,0,l=51)),zlim=pred.lim)
contour(sk,add=T,nlev=6)
points(elevation,add=TRUE,cex.max=2)
               
# 25/59
x=seq(0,1,l=101)
plot(x,cov.spatial(x,cov.model="mat",kappa=0.5,cov.pars=c(1,0.25)),type="l")

# 27/59
x=seq(0,1,l=101)
plot(x,cov.spatial(x,cov.model="spherical",cov.pars=c(1,0.75)),type="l")
plot(x,cov.spatial(x,cov.model="wave",cov.pars=c(1,0.05)),type="l")

# 29/59
image(grf(20^2,grid="reg",cov.pars=c(1,0.13),cov.model="mat",kappa=2.5),col=gray(seq(1,0,l=51)),xlab=" ", ylab=" ")

# 36/59
args(krige.bayes)
function(geodata, coords = geodata$coords,
    data = geodata$data,
    locations = "no", borders, model,prior,output)
  
# 37/59
MC=model.control(trend.d="1st",trend.l="1st",kappa=1.5)
PC=prior.control(phi.discrete=seq(0,6,l=21),
    phi.prior="reciprocal",
    tausq.rel.prior="unif",tausq.rel.discrete=seq(0,1,l=11))
OC=output.control(n.post=1000,moments=T)

# Now krige.bayes can be used:
set.seed(268)
locs=pred_grid(c(0,6.3),c(0,6.3),by=0.3)
skb = krige.bayes(elevation,loc=locs,model=MC,prior=PC,output=OC)
par(mfrow=c(2,1));plot(skb)
par(mfrow=c(3,2));hist(skb)

# 40/59
# Map of the mean of the predictive dist.
image(skb,col=gray(seq(1,0,l=51)))
points(elevation,add=TRUE,cex.max=2)
# Map of the variance of the predictive dist.
image(skb,val="variance",col=gray(seq(1,0,l=51)),main="prediccion var")
contour(skb,val="variance",add=T)
points(elevation$coords,pch="+")
# Sample of the predictive dist.
image(skb,val="simulation",col=gray(seq(1,0,l=51)),number.col=1,main="simulacion")
contour(skb,val="simulation",number.col=1,add=T)

# 43/59
# Predictions for some points
names(skb$predictive)
dim(skb$predictive$simulations)
pred=skb$predictive$simulations
par(mfrow=c(3,2))
for(i in 1:6)
  {hist(pred[,i],prob=T,xlab=" ",main=paste("Prediction on",i)) }
names(skb$predictive)

# 46/59
sim <- grf(grid =expand.grid(x=seq(1,10,l = 10),y =seq(1, 10, l = 10)),cov.pars = c(0.1, 0.2))
attr(sim,"class") <- "geodata"
sim$units.m <- c(rep(5, 100))
sim$data <- rpois(100, lambda =sim$units.m*exp(sim$data))
plot(sim$coords[,1],sim$coords[,2],type = "n")
text(sim$coords[,1],sim$coords[,2],format(sim$data))

# 48/59
library(geoRglm)#!
model2=krige.glm.control(cov.pars =c(1,1),beta = 1)
test2.tune=pois.krige(p50,krige=model2,mcmc.input=list(S.scale=0.2,thin = 1))

# 49/59
test2.tune=pois.krige(p50, krige = model2, mcmc.input = list(S.scale = 0.5, thin = 1))
par(mfrow=c(1,2))
plot(log(test2.tune$intensity[45,]), type = "l")
acf(log(test2.tune$intensity[45,]), type = "correlation", plot = TRUE)

# 51/59
test2=pois.krige(p50,locations=cbind(c(0.5,0.5),c(1,0.4)),
    krige=model2, mcmc.input=mcmc.control(S.scale = 0.5),
    output=output.glm.control(sim.predict = TRUE))

# 52/59
test2$predict
test2$krige.var
test2$mcmc.error

par(mfrow = c(1,2))
hist(test2$simulations[1,],main="(0.5, 0.5)")
hist(test2$simulations[2,],main ="(1, 0.4)")

model2.u=krige.glm.control(cov.pars = c(1,1),beta = 1,type.krige = "ok")
test2.unif.beta =pois.krige(p50,krige =model2.u,mcmc.input = list(S.scale = 0.5))

# 54/59
prior5=prior.glm.control(phi.prior="fixed", phi = 0.1)
mcmc5.tune=mcmc.control(S.scale = 0.01, thin = 1)
test5.tune=pois.krige.bayes(p50, prior = prior5, mcmc.input = mcmc5.tune)
# S.scale parameter adjusted to improve acceptance rate.
mcmc5=mcmc.control(S.scale = 0.075, thin = 100)
test5=pois.krige.bayes(p50, locations = t(cbind(c(2.5,3),c(-6050,-3270))),
    prior = prior5, mcmc.input = mcmc5,
    output = list(threshold=10,
    quantile = c(0.05,0.99)))

# 55/59
par(mfrow = c(1,3))
hist(test5$posterior$simulations[10,],main="(9, 0)")
hist(test5$posterior$simulations[23,],main= "(2,2)")
hist(test5$posterior$simulations[36,],main="(5,3)")

# 57/59
# Example with uncertainty in other parameters
mcmc6.tune<-mcmc.control(S.scale=0.075,n.iter=2000,thin=100,phi.scale = 0.01)
prior6<-prior.glm.control(phi.prior="uniform",phi.discrete=seq(0.02,1, 0.02),tausq.rel = 0.05)
test6.tune <- pois.krige.bayes(p50,prior = prior6,mcmc.input= mcmc6.tune)
#This may take some time.
mcmc6 <-mcmc.control(S.scale=0.075,n.iter=400000,thin=200,burn.in=5000,phi.scale=0.12,phi.start=0.5)
test6 <- pois.krige.bayes(p50,locations = t(cbind(c(2.5,3.5),c(-60,-37))),prior=prior6,mcmc.input=mcmc6)
# some posterior distributions.
par(mfrow=c(1,3))
hist(test6$posterior$beta$sample,main="beta")
hist(test6$posterior$sigmasq$sample,main="sigmasq")

