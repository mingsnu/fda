## This practice code is modified from the 'fdarm-ch**.R' in the 'fda' package:
## system.file('scripts', package='fda')

##  load the fda package
library(fda)

## ---- ex1
unitRng = c(0,1)
const.basis = create.constant.basis(unitRng)
monom.basis = create.monomial.basis(unitRng, nbasis=1)
fourier.basis = create.fourier.basis(unitRng, nbasis=5, period=1)
bspline.basis = create.bspline.basis(unitRng,
    nbasis=5, norder=2, breaks=seq(0, 1, length=5) )
const.basis
bspline.basis
plot(const.basis)
plot(monom.basis)
plot(fourier.basis)
plot(bspline.basis)

## ---- ex2-1
yearRng = c(0,365)
daybasis65 = create.fourier.basis(yearRng, 65)
## same as
daybasis65 = create.fourier.basis(c(0,365), 65, 365)

## ---- ex2-2
zerobasis  = create.fourier.basis(c(0, 365), 65, dropind=1)

## ---- ex2-3
str(daybasis65)
str(zerobasis)

## ---- ex3
##  order 4 spline, one interior knot
## there are many ways to specify one particula basis
b1 = create.bspline.basis(c(0, 1), 5)
b2 = create.bspline.basis(breaks=c(0, .5, 1))
identical(b1,b2)
## check knots by using the 'knots' function
knots(b1, interior=FALSE)
## draw the basis by using plot
plot(b1, lwd=2)

## Order 4 with 3 knots at 0.5 allows the first derivative to change
## abruptly at 0.5 while the function remains continuous
## Note: the basis lines in the plot can only use 6 distinguish colors,
## colors are recursively used if there are more than 6 basis lines.
b3 = create.bspline.basis(breaks=c(0, 0.5, 0.5, 0.5, 1))
knots(b3, interior=FALSE)
plot(b3, lwd=2)

## ---- ex4-1
splinebasis = create.bspline.basis(c(0,10), 13)
plot(splinebasis, xlab='t', ylab='Bspline basis functions B(t)',
     las=1, lwd=2)

## ---- ex4-2
basis2 = create.bspline.basis(c(0,2*pi), 5, 2)
basis3 = create.bspline.basis(c(0,2*pi), 6, 3)
basis4 = create.bspline.basis(c(0,2*pi), 7, 4)

theta = seq(0, 2*pi, length=201)
sin.theta = sin(theta)

sin2 = Data2fd(theta, sin.theta, basis2)
sin3 = Data2fd(theta, sin.theta, basis3)
sin4 = Data2fd(theta, sin.theta, basis4)

sin2.theta = predict(sin2, theta)
sin3.theta = predict(sin3, theta)
sin4.theta = predict(sin4, theta)

sinRng = range(sin2.theta)
pi3 = ((1:3)*pi/2)

par(mfrow=c(3,2), mar=c(3,4,2,2)+.1)

plot(theta, sin2.theta, type='l', ylim=sinRng, xlab='', ylab='Order = 2',
     main='sine(t)' )
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin2.theta = predict(sin2, theta, 1)
plot(theta, Dsin2.theta, type='l', ylim=sinRng, xlab='', ylab='',
     main='D sine(t)')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

plot(theta, sin3.theta, type='l', ylim=sinRng, xlab='', ylab='Order = 3')
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin3.theta = predict(sin3, theta, 1)
plot(theta, Dsin3.theta, type='l', ylim=sinRng, xlab='', ylab='')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

plot(theta, sin4.theta, type='l', ylim=sinRng, xlab='t', ylab='Order = 4')
lines(theta, sin.theta, lty='dashed')
abline(v=pi3, lty='dotted')

Dsin4.theta <- predict(sin4, theta, 1)
plot(theta, Dsin4.theta, type='l', ylim=sinRng, xlab='t', ylab='')
lines(theta, cos(theta), lty='dashed')
abline(v=pi3, lty='dotted')

## ---- ex5
## sum of basis function values (all equal to 1)
t = seq(0, 10, 0.5)

## calculate the function values corresponding to each basis at each
## time point t. Or use 'p = eval.basis(t, splinebasis)'.
p = predict(splinebasis, t)
rowSums(p)

## ---- ex6
daybasis65 = create.fourier.basis(c(0,365), 65)
## dummy coefmat
coefmat = matrix(0, 65, 35, dimnames=list(
     daybasis65$names, CanadianWeather$place) )
tempfd = fd(coefmat, daybasis65)
class(tempfd)
str(tempfd)

## ---- ex7
##  Two order 2 splines over unit interval
unitRng = c(0, 1)
bspl2 = create.bspline.basis(unitRng, norder=2)
plot(bspl2, lwd=2)

f1 = fd(c(-1, 1), bspl2)
f2 = fd(c(1, 3), bspl2)
fsum = f1 + f2
fdif = f2 - f1
fpro = f1*f2
fsqr = f1^2

plot(f1,   lwd=2, xlab="", ylab="Line 1")
plot(f2,   lwd=2, xlab="", ylab="Line 2")
plot(fsum, lwd=2, xlab="", ylab="Line 1 + Line2")
plot(fdif, lwd=2, xlab="", ylab="Line 2 - Line1")
plot(fpro, lwd=2, xlab="", ylab="Line 1 * Line2")
plot(fsqr, lwd=2, xlab="", ylab="Line 1 ^ 2")
## sqare root for f2; can't do it for f1
plot(f2^0.5, lwd=2, xlab="", ylab="Line 2 ^ 0.5")
## reciprocal of a positive function with no values near zero
## can't do it for f1
plot(f2^(-1), lwd=2, xlab="", ylab="Line 2 ^ (-1)")

## use fdnames argument
f3 = fd(c(-1, 1), bspl2, list("aaa", "bbb", "ccc"))
## x axis label should be 'aaa', y axis label should be 'ccc'
plot(f3)

## ---- ex8
## mean and sum method
## mean & sum can be applied to a set of functions
## We'll talk about 'smooth' function later. Here the smooth function
## is used to calculate the coefficient matrix
daybasis65 = create.fourier.basis(c(0,365), 65)
tempfd = smooth.basis(day.5,
          CanadianWeather$dailyAv[,,'Temperature.C'], daybasis65)$fd
meanTempfd = mean(tempfd)
sumTempfd  = sum(tempfd)
plot(tempfd)
lines(meanTempfd, lwd = 3, col = "red")
## can not use 'sumTempfd/35'
lines(sumTempfd*(1/35), lwd = 3, lty=2, col = "blue")

## `[` operator
## you can use `[` operator to choose which line to plot
plot(tempfd[1:2])

## predict method
t = day.5[sample(1:length(day.5), 10)]
## can also use 'p = eval.fd(t, meanTempfd)'
p = predict(meanTempfd, t)
points(t, p, pch=19, col="green")

## boxplot
boxplot(tempfd)

## ---- ex9-1
## Smoothing Using Regression Analysis
MtlDaily = as.matrix(MontrealTemp)
thawdata = t(MtlDaily[,16:47])

## only use part of the data
daytime = ((16:47)+0.5)
plot(daytime, apply(thawdata,1,mean), "b", lwd=2,
     xlab="Day", ylab="Temperature (deg C)", cex=1.2)
## create basis
thawbasis = create.bspline.basis(c(16,48),7)
## the 'phi' matrix
thawbasismat = predict(thawbasis, daytime)
## calculate the coefficient matrix
## crossprod(x, y) is the same as 't(x) %*% y'
## crossprod(x) is the same as 't(x) %*% x'
## same as using 'tmp=solve(t(thawbasismat)%*%thawbasismat)%*%t(thawbasismat)%*%thawdata'
thawcoef = solve(crossprod(thawbasismat),
    crossprod(thawbasismat,thawdata))

## build the functional data object
thawfd = fd(thawcoef, thawbasis,
    list("Day", "Year", "Temperature (deg C)"))
plot(thawfd, lty=1, lwd=2, col=1)
## fitted 'mean line' on the mean scatterplot
plot(daytime, apply(thawdata,1,mean), xlab="Day", ylab="Temperature (deg C)")
lines(mean(thawfd), lty=1, lwd=2, col=1)

## compare a curve to the data from which it was estimated by using
## the 'plotfit.fd' function.
plotfit.fd(thawdata[,1], daytime, thawfd[1],
           lty=1, lwd=2, main='')

## ---- ex9-2
thawSmooth = smooth.basis(daytime, thawdata, thawbasis)
## compare the difference between 'thawSmooth$y2cMap%*%thawdata'
## and 'thawcoef' we calculated earlier, 
sum((thawSmooth$y2cMap%*%thawdata-thawcoef)^2)

## ---- ex10-1
Lbasis  = create.constant.basis(c(0,365));  #  create a constant basis
Lcoef   = matrix(c(0,(2*pi/365)^2,0),1,3)   #  set up three coefficients
wfdobj  = fd(Lcoef,Lbasis)      # define an FD object for weight functions
wfdlist = fd2list(wfdobj)       # convert the FD object to a cell object
harmaccelLfd0 = Lfd(3, wfdlist)  #  define the operator object

## ---- ex10-2
## use 'int2Lfd' to create the acceleration object
accelLfd = int2Lfd(2)

## use 'vec2Lfd' to create the harmonic acceleration object and
## the result is the same with 'harmaccelLfd0'
## Lcoef can be a vector
harmaccelLfd <- vec2Lfd(Lcoef, c(0,365))

## back to the temperature example introduced above
LmeanTempVec = eval.fd(day.5, meanTempfd, harmaccelLfd)
plot(day.5, LmeanTempVec, type="l", cex=1.2,
     xlab="Day", ylab="Harmonic Acceleration")
abline(h=0)

## ---- ex10-3
## Second derivative to meanTempfd
D2tempfd = deriv.fd(meanTempfd, 2)
## Apply harmonic acceleration operator to meanTempfd
Ltempfd  = deriv.fd(meanTempfd, harmaccelLfd)
plot(Ltempfd)

## ---- ex11
