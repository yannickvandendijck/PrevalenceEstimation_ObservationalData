###---------------------------------------------###
###								###
### Prevalence and Trend Estimation from 		###
### Observational Data with Highly Variable 	###
### Post-Stratification Weights			###
###								###
### Y. Vandendijck, C. Faes, N. Hens		###
###								###
### Supplementary Materials: Appendix C		###
###								###
###---------------------------------------------###

### Libraries
library(mvtnorm)
### Functions
expit <- function(x){exp(x)/(1+exp(x))}


### Some input parameters
sigma.sq = 0.02
w0 = 3
beta0 = -2
beta1 = 0.2
seed = 4321


### Define population and sample size
Nh = c(7500, 7500, 8000, 8000, 8500, 8500, 9000, 9000, 9000,
		9000,9000,9000,8500,8500,8000,8000,7500,7500)
nh =  c(5,15,30,75,125,150,200,275,375,
		375,275,200,150,100,80,40,20,10)
N =sum(Nh); n = sum(nh)
Ph = Nh/N; ph = nh/n
wh = Ph/ph
Xh = 1:18


### Simulate data according to LIN model
set.seed(seed)
D = sigma.sq*diag(length(Xh))
delta.h = beta0 + beta1*Xh
delta.star.h = rmvnorm(1, mean = delta.h, sigma=D)
true.prop.h = expit(delta.star.h)
Y.ih = rbinom(length(Xh), Nh, true.prop.h)
population.Yh = Y.ih/Nh
population.prevalence = sum(population.Yh*Nh)/N

y.ih = rbinom(length(Xh), nh, population.Yh)
simulated.data = data.frame(cbind(y.ih, nh, Xh))
names(simulated.data) = c("resp", "nh", "xh")


### Export to SAS to fit model
dir = "F:/PhD/Weighting/Paper/Paper major revision 1 files/Appendix C/"
filename = "PrevalenceExample.txt"
write.table(simulated.data, paste(dir,filename,sep=""),row.names=F)


### Calculate Points Estimates
filename = "InputSAS1.txt"
import1 = read.table(paste(dir,filename,sep=""), header=T)

nh.hat = rep(0,length(Xh))
nh.hat[wh > w0] = (Nh[wh > w0]/N)*(n/w0)
gamma = (n-sum(nh.hat)) / sum(nh[wh <= w0])
nh.hat[wh <= w0] = gamma*nh[wh <= w0]
pi.hat = nh.hat/Nh

sample.prop = y.ih / nh
pred.prop = import1$prob
y.ws = 1/N * sum(nh*sample.prop + (Nh-nh)*pred.prop)
y.ws.GREG = 1/N * sum((nh/pi.hat)*sample.prop + (Nh-(nh/pi.hat))*pred.prop)


### Calculate vairance y.ws analytically
filename = "InputSAS2.txt"
import2 = read.table(paste(dir,filename,sep=""), header=T)

X = cbind(rep(1,length(Xh)), Xh)
Z = diag(length(Xh))
C = cbind(X, Z)
G = import2$Estimate*diag(18)
B = matrix(0, length(Xh)+2, length(Xh)+2)
B[3:(length(Xh)+2), 3:(length(Xh)+2)] = solve(G)
eta = import1$eta
Delta = diag( as.vector(exp(eta)/((1+exp(eta))^2)) )
A = diag( as.vector(expit(eta)*(1-expit(eta))) )
Gamma = diag( 1/nh )
Sigma.p = solve(Delta)%*%A%*%solve(Delta)%*%Gamma
Theta = Delta%*%C%*% solve(t(C)%*%solve(Sigma.p,C)+B) %*%t(C)%*%Delta
V1 = Nh-nh
var.analytical.y.ws = t(V1)%*%Theta%*%V1 / (N^2)





### Calculate variance y.ws with bootstrap
# generate bootstrap samples
D = import2$Estimate*diag(length(Xh))
V = Z%*%G%*%t(Z) + Sigma.p
pseudo.data = solve(Delta)%*%(sample.prop - import1$prob) + eta
beta.hat = solve( t(X)%*%solve(V,X) , t(X)%*%solve(V,pseudo.data) )

y.boot = vector()
y.pop.mean = vector()
B = 250
set.seed(seed)
for (b in 1:B){
	tilde.u = rmvnorm(1, mean=rep(0,length(Xh)), sigma=D)
	tilde.delta.star = X%*%beta.hat + Z%*%as.vector(tilde.u)
	population.b = rbinom(18,round(Nh),expit(tilde.delta.star))
	y.pop.mean[b] = sum(population.b)/sum(Nh)
	sample.b = rbinom(length(Xh), nh, (population.b / Nh))
	y.boot = c(y.boot, sample.b)
}

n.boot <- rep(nh,times=B)
xh.boot <- rep(Xh,times=B)
b.boot <- rep(1:B,each=length(Xh))

bootstrap.data <- as.data.frame(cbind(y.boot,n.boot,h.boot,b.boot))
names(bootstrap.data) <- c("y","nh","xh","b")
filename = "BootstrapData.txt"
write.table(bootstrap.data, paste(dir,filename,sep=""),row.names=F)

# Analyse results bootstrap samples
filename = "BootstrapInputSAS1.txt"
import3 = read.table(paste(dir,filename,sep=""), header=T)

bootstrap.y.ws <- vector()
for (b in 1:B){
	set.b <- import3[import3$b==b,]
	ph.boot.fit <- set.b$prob
	ph.boot <- set.b$y/set.b$nh
	if (set.b$prob[1]==999){bootstrap.y.ws[b]=NA}
	else {bootstrap.y.ws[b]=(1/N)*(sum( nh*ph.boot ) + sum( (Nh-nh)*ph.boot.fit ))}
}

v.k = (bootstrap.y.ws - y.pop.mean)^2
var.y.ws.boot = mean(v.k, na.rm=T)





### Calculate jackknife variance y.ws.greg
# construct bootstrap samples
G=250
s=n/G
y.jack = NULL
set.seed(seed)
for (h in 1:length(Xh)){
	yy=c(rep(1,round(y.ih[h])), rep(0,round(nh[h]-y.ih[h])))
	y.jack = c(y.jack, sample(yy))
}
group.ind = rep(1:G, times=ceiling(s))[1:n]
jack.data.full = cbind(y.jack, group.ind, rep(1:length(Xh), times=nh))

y.boot=n.boot=h.boot=b.boot=NULL
for (g in 1:G){
	set = jack.data.full[group.ind!=g,]
	y.boot = c(y.boot, as.vector(tapply(set[,1],set[,3],sum)))
	n.boot = c(n.boot, as.vector(tapply(set[,1],set[,3],length)))
	h.boot = c(h.boot, Xh)
	b.boot = c(b.boot, rep(g,length(Xh)))
}

jackknife.data <- as.data.frame(cbind(y.boot,n.boot,h.boot,b.boot))
names(jackknife.data) <- c("y","nh","xh","b")
filename = "JackknifeData.txt"
write.table(jackknife.data, paste(dir,filename,sep=""),row.names=F)

# Analyse results jackknife samples
filename = "JackknifeInputSAS1.txt"
import4 = read.table(paste(dir,filename,sep=""), header=T)

jackknife.y.ws.greg <- vector()
for (g in 1:G){
	set.g <- import4[import4$b==g,]
	ph.boot <- set.g$y/set.g$nh
	ph.boot.fit <- set.g$prob
	jackknife.y.ws.greg[g] <- (1/N)*(sum( (set.g$nh/pi.hat)*ph.boot ) + sum( (Nh-set.g$nh/pi.hat)*ph.boot.fit ))
}

v.k = (jackknife.y.ws.greg - mean(jackknife.y.ws.greg))^2
var.y.ws.greg = (G-1)/G*sum(v.k)



