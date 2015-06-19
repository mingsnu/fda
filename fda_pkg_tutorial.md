

# FDA package tutorial

Recall that a function $x(t)$ can be defined as a linear combination of $K$ known *basis functions*:

$$x(t) = \sum_{k=1}^K c_k\phi_k(t) = \mathbf{c}^{'}\mathbf{\phi}(t)$$

We often want to consider a sample of $N$ functions, $x_i(t) =
\sum_{k=1}^K c_{ik}\phi_k(t), i=1,...,N$, in matrix notation:

$$\mathbf{x}(t) = \mathbf{C}\mathbf{\phi}(t)$$

where $\mathbf{x}(t)$ is a vector of length $N$ containing the functions $x_i(t)$, and the coefficient matrix $C$ has $N$ rows and $K$ columns.

So functions are built in two stages:

1. Define a set of basis functions $\mathbf{\phi}_k$.

2. Then set up a vector, matrix, or array of coefficients to define the function as
a linear combination of these basis functions.

## Specify Basis Systems

The Fourier basis system is the usual choice for periodic functions,
and the spline basis system (and bsplines in particular) tends to
serve well for nonperiodic functions.

### `create` functions

The `create` functions in R that set up constant, monomial, fourier,
bspline, exponential, polygonal and power basis systems: (many
optional arguments are omitted here)

    basisObj = create.constant.basis(rangeval)
	basisObj = create.monomial.basis(rangeval, nbasis)
	basisObj = create.fourier.basis(rangeval, nbasis, period)
	basisObj = create.bspline.basis(rangeval, nbasis, norder, breaks)
	basisObj = create.exponential.basis(rangeval, nbasis)
	basisObj = create.polygonal.basis(rangeval)
	basisObj = create.power.basis(rangeval, nbasis)

The `basisObj` is called a **functional basis object** and it has a
`basisfd` class attribute.

**Some common arguments**

*rangeval*: a numeric vector of length 2 defining the interval over which
the functional data object can be evaluated

*nbasis*: an integer variable specifying the number of basis
 functions ($K$).

**Return value**

The `create.**.basis` function returns a `list` object with a
`basisfd` class attribute.

**Examples**

<div class="chunk" id="ex1"><div class="rcode"><div class="source"><pre class="knitr r">unitRng = c(0,1)
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
</pre></div>
</div></div>

### Fourier basis

Here is the complete calling sequence in R for the `create.fourier.basis` function:

    create.fourier.basis(rangeval = c(0, 1), nbasis = 3, period = diff(rangeval), 
                         dropind = NULL, quadvals = NULL, values = NULL, basisvalues = NULL, 
                         names = NULL, axes = NULL) 

where only the first three arguments are required to define a Fourier basis system. Some of the arguments are explained as follows:

- *nbasis*: must be an odd integer (or will be increased by 1 automatically).

- *period*: the width of any interval over which the Fourier functions
repeat themselves or are periodic, default is `diff(rangeval)`

- *dropind*: an optional vector of integers specifiying basis functions to be dropped.

**Examples**

1. The following code sets up a Fourier basis with $K=65$
basis functions with a period of 365 (days):
<div class="chunk" id="ex2-1"><div class="rcode"><div class="source"><pre class="knitr r">yearRng = c(0,365)
daybasis65 = create.fourier.basis(yearRng, 65)
## same as
daybasis65 = create.fourier.basis(c(0,365), 65, 365)
</pre></div>
</div></div>

2. The following code sets up a Fourier basis that does not include the
initial constant term:
<div class="chunk" id="ex2-2"><div class="rcode"><div class="source"><pre class="knitr r">zerobasis  = create.fourier.basis(c(0, 365), 65, dropind=1)
</pre></div>
</div></div>

3. What's the difference between `daybasis65` and `zerobasis`? 
<div class="chunk" id="ex2-3"><div class="rcode"><div class="source"><pre class="knitr r">str(daybasis65)
str(zerobasis)
</pre></div>
</div></div>

### Spline basis

Recall the concepts introduced before, `break points`, `knots` and `order`.

> By default, and in the large majority of applications, there will be only a single knot at every break point except for the boundary values.

> The end points, however, are assigned **as many knots as** the order of the spline, implying that the function value will, typically, drop to zero outside of the interval over which the function is defined.

For example, if we define a function over [0,1] with a single interior break point at 0.5 with cubic splinebasis (order 4), then the knots are (0, 0, 0, 0, 0.5, 1, 1, 1, 1). (**Question**: What if we want to allow the first derivative to change abruptly at 0.5 while the function remains continuous?)

> You will not have to worry about those multiple knots at the end points; the code takes care of this automatically. You will be typically constructing spline functions where you will only have to supply break points, and if these break points are equally spaced, you will not even have to supply these.

The number of basis functions is determined by the relation:

$$number\ of\ basis\ functions = order + number\ of\ interior\ knots$$

**Examples**

<div class="chunk" id="ex3"><div class="rcode"><div class="source"><pre class="knitr r">##  order 4 spline, one interior knot
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
</pre></div>
</div></div>

#### B-splines

Here is the complete calling sequence in R for the `create.bspline.basis` function:

    create.bspline.basis(rangeval = NULL, nbasis = NULL, norder = 4, breaks = NULL, 
                         dropind = NULL, quadvals = NULL, values = NULL, basisvalues = NULL, 
                         names = "bspl") 

Some of the arguments are explained as follows:
  
- *norder*: an integer specifying the order of b-splines, which is one
higher than their degree. The default of 4 gives cubic splines.

- *breaks*: a vector specifying the break points defining the b-spline.
Also called knots, these are a strictly increasing sequence of
junction points between piecewise polynomial segments.

**Examples**

The 13 spline basis functions defined over the interval [0,10] by nine
interior knots.

<div class="chunk" id="ex4-1"><div class="rcode"><div class="source"><pre class="knitr r">splinebasis = create.bspline.basis(c(0,10), 13)
plot(splinebasis, xlab='t', ylab='Bspline basis functions B(t)',
     las=1, lwd=2)
</pre></div>
<div class="rimage default"><img src="figure/ex4-1-1.png" title="plot of chunk ex4-1" alt="plot of chunk ex4-1" class="plot" /></div>
</div></div>

The following graph is used to illustrate the role of the order of a spline.

<div class="chunk" id="ex4-2"><div class="rcode"><div class="rimage default"><img src="figure/ex4-2-1.png" title="plot of chunk ex4-2" alt="plot of chunk ex4-2" class="plot" /></div>
</div></div>

> In general, if we need smooth and accurate derivatives, we need to
> increase the order of the spline. A useful rule to remember is to
> fix the order of the spline basis to be at least two higher than the
> highest order derivative to be used.

For example, in the image above, if we want the first derivative of
the function to be continuous, we need the order to be $1+2=3$; if we
want the first derivative of the function to be smooth(ie. second
derivative to be continuous), we need the order to be $2+2=4$.

> The B-spline basis system has a property that is often useful: the
> sum of the B-spline basis function values at any point $t$ is equal
> to one.

**Example**

<div class="chunk" id="ex5"><div class="rcode"><div class="source"><pre class="knitr r">## sum of basis function values (all equal to 1)
t = seq(0, 10, 0.5)

## calculate the function values corresponding to each basis at each
## time point t. Or use 'p = eval.basis(t, splinebasis)'.
p = predict(splinebasis, t)
rowSums(p)
</pre></div>
</div></div>

> Although spline basis functions are wonderful in many respects, they
> tend to produce rather unstable fits to the data near the beginning
> or the end of the interval over which they are defined. This is
> because in these regions we run out of data to define them, so at
> the boundaries the spline function values are entirely determined by
> a single coefficient. This boundary instability of spline fits
> becomes especially serious for derivative estimation, and the higher
> the order of the derivative, the wilder its behavior tends to be at
> the two boundaries.

## Build functional data objects

In this section, we are going to define a *functional data object* by
combining a set of coefficients $c_k$ with previously defined basis
system.

### Define functions

In R, once we have defined both *functional basis object* and a
*coefficient vector* (matrix or array), we can define the *functional
data object* by using the `fd` 'constructor' function. However, in
practice we seldom need to use the `fd` function directly because other
advanced functions would call it for us. Even though understanding
this function should be helpful.

    fd(coef=NULL, basisobj=NULL, fdnames=NULL)

**Description**

     This is the constructor function for objects of the ‘fd’ class.
     Each function that sets up an object of this class must call this
     function.  This includes functions ‘Data2fd’, ‘smooth.basis’,
     ‘density.fd’, and so forth that estimate functional data objects
     that smooth or otherwise represent data.  Ordinarily, users of the
     functional data analysis software will not need to call this
     function directly, but these notes are valuable to understanding
     the components of a ‘list’ of class ‘fd’.

**Arguments**

- *coef*: a vector, matrix, or three-dimensional array of
  coefficients.  The first dimension (or elements of a vector)
  corresponds to basis functions.  A second dimension corresponds to
  the number of functional observations, curves or replicates.  If
  ‘coef’ is a vector, it represents only a single functional
  observation.  If ‘coef’ is an array, the third dimension corresponds
  to variables for multivariate functional data objects.

- *basisobj*: a functional basis object defining the basis

- *fdnames*: A list of length 3, each member being a string vector
containing labels for the levels of the corresponding dimension of the
discrete data.  The first dimension is for argument values, and is
given the default name "time", the second is for replications, and is
given the default name "reps", and the third is for functions, and is
given the default name "values".

By default, `plot` function uses the `fdnames` values for labels.

**Examples**

<div class="chunk" id="ex6"><div class="rcode"><div class="source"><pre class="knitr r">daybasis65 = create.fourier.basis(c(0,365), 65)
## dummy coefmat
coefmat = matrix(0, 65, 35, dimnames=list(
     daybasis65$names, CanadianWeather$place) )
tempfd = fd(coefmat, daybasis65)
class(tempfd)
str(tempfd)
</pre></div>
</div></div>

### Methods for functional data objects (`fd` class)

Arithmetic methods defined for the functional data objects includes `+`, `-`, `*` and `^`.

Let's see the examples directly.

**Examples**

<div class="chunk" id="ex7"><div class="rcode"><div class="source"><pre class="knitr r">##  Two order 2 splines over unit interval
unitRng = c(0, 1)
bspl2 = create.bspline.basis(unitRng, norder=2)
plot(bspl2, lwd=2)
</pre></div>
<div class="rimage default"><img src="figure/ex7-1.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="source"><pre class="knitr r">f1 = fd(c(-1, 1), bspl2)
f2 = fd(c(1, 3), bspl2)
fsum = f1 + f2
fdif = f2 - f1
fpro = f1*f2
fsqr = f1^2

plot(f1,   lwd=2, xlab="", ylab="Line 1")
</pre></div>
<div class="rimage default"><img src="figure/ex7-2.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
<div class="source"><pre class="knitr r">plot(f2,   lwd=2, xlab="", ylab="Line 2")
</pre></div>
<div class="rimage default"><img src="figure/ex7-3.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
<div class="source"><pre class="knitr r">plot(fsum, lwd=2, xlab="", ylab="Line 1 + Line2")
</pre></div>
<div class="rimage default"><img src="figure/ex7-4.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
<div class="source"><pre class="knitr r">plot(fdif, lwd=2, xlab="", ylab="Line 2 - Line1")
</pre></div>
<div class="rimage default"><img src="figure/ex7-5.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
<div class="source"><pre class="knitr r">plot(fpro, lwd=2, xlab="", ylab="Line 1 * Line2")
</pre></div>
<div class="rimage default"><img src="figure/ex7-6.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
<div class="source"><pre class="knitr r">plot(fsqr, lwd=2, xlab="", ylab="Line 1 ^ 2")
</pre></div>
<div class="rimage default"><img src="figure/ex7-7.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
<div class="source"><pre class="knitr r">## sqare root for f2; can't do it for f1
plot(f2^0.5, lwd=2, xlab="", ylab="Line 2 ^ 0.5")
</pre></div>
<div class="warning"><pre class="knitr r">## Warning in seq.default(rangeval[1], rangeval[2], n = 11): extra argument
## 'n' will be disregarded
</pre></div>
<div class="rimage default"><img src="figure/ex7-8.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
<div class="source"><pre class="knitr r">## reciprocal of a positive function with no values near zero
## can't do it for f1
plot(f2^(-1), lwd=2, xlab="", ylab="Line 2 ^ (-1)")
</pre></div>
<div class="warning"><pre class="knitr r">## Warning in seq.default(rangeval[1], rangeval[2], n = 11): extra argument
## 'n' will be disregarded
</pre></div>
<div class="rimage default"><img src="figure/ex7-9.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
<div class="source"><pre class="knitr r">## use fdnames argument
f3 = fd(c(-1, 1), bspl2, list("aaa", "bbb", "ccc"))
## x axis label should be 'aaa', y axis label should be 'ccc'
plot(f3)
</pre></div>
<div class="rimage default"><img src="figure/ex7-10.png" title="plot of chunk ex7" alt="plot of chunk ex7" class="plot" /></div>
<div class="output"><pre class="knitr r">## [1] &quot;done&quot;
</pre></div>
</div></div>

Other useful methods defined for the `fd` class includes: `coef`, `density`, `deriv`, `[`, `lines`, `mean`, `sum`, `predict` and so on.

<div class="chunk" id="ex8"><div class="rcode"><div class="source"><pre class="knitr r">## mean and sum method
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
</pre></div>
</div></div>

### Smoothing Using Regression Analysis

Recall that $x(t) \approx \sum_{k=1}^K c_k\phi_k(t)$ and the ordinary least squares estimate is $\hat{\mathbf{c}} = (\mathbf{\phi}^T\mathbf{\phi})^{-1}\mathbf{\phi}^T\mathbf{y}$.

Next we use an example to illustrate how to use this method to smooth the functional data

*January Thaw Data* In the following example, we will use 34 years of daily temperature data for Montreal, extracts temperatures for January 16th to February 15th. The data looks like this:

           1961  1962  1963 1964  1965  1966  1967  1968  1969  1970 ...
    jan16 -10.9  -7.5 -16.2 -4.2 -22.5 -15.9 -13.9 -21.4 -15.8 -17.5 ...
	jan17  -6.1 -13.1 -12.0 -8.1 -22.0  -7.8  -6.4 -17.8  -7.3  -8.9 ...
	jan18 -16.1 -15.9 -12.3 -4.2 -16.1  -6.4 -19.8  -8.9   0.8 -17.3 ...
	jan19 -21.7 -13.1  -8.1  0.3 -13.6  -3.6 -17.0   1.1  -2.5 -23.6 ...
	jan20 -19.2 -12.8  -3.4  2.0  -7.2  -1.2  -9.7   0.3  -9.5 -23.9 ...
	...    ...   ...   ...   ...   ...   ...   ...   ...   ...  ...

<div class="chunk" id="ex9"><div class="rcode"><div class="source"><pre class="knitr r">## Smoothing Using Regression Analysis
MtlDaily = as.matrix(MontrealTemp)
thawdata = t(MtlDaily[,16:47])

## only use part of the data
daytime = (16:47)
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
</pre></div>
</div></div>

### The Linear Differential Operator or `Lfd` Class

One of the most distinctive feature of functional data analysis is that we can take advantage of the derivative of a function.

Now we introduce a new concept called **linear differential operators**. When we apply a linear differential operatior to a function, we get a inear combinations of derivatives of the function. A general expression is

$$Lx(t) = \beta_0(t)x(t) + \beta_1(t)Dx(t) + \cdots + \beta_{m-1}(t)D^{m-1}x(t) + D^{m}x(t)$$

where the known *linear differential operator coefficient functions* $\beta_j(t), j=0,..., m-1$ are either constants or functions.

Some special examples are acceleration, $Lx=D^2x$, and harmonic acceleration $L=\omega^2D + D^3$, where $\omega$ is the period. (For now we just focus on how it is defined, we'll go through more detail about harmonic acceleration operator later).

The `Lfd` class is defined by a constructor function `Lfd`:

    Lfd(nderiv=0, bwtlist=vector("list", 0))

**Description**

    A linear differential operator of order m is defined, usually to specify a roughness penalty.

     
**Arguments**

- *nderiv*: a nonnegative integer specifying the order $m$ of the highest
          order derivative in the operator

- *bwtlist*: a list of length $m$.  Each member contains a **functional data
          object** that acts as a weight function for a derivative.  The
          first member weights the function, the second the first
          derivative, and so on up to order $m-1$.

**Examples**

Consider the harmonic acceleration operator where $m = 3$, $\beta_0 = \beta_2 = 0$ and $\beta_1=\omega^2$. In R we could define the harmonic acceleration object `harmaccelLfd0` in this way:

<div class="chunk" id="ex10-1"><div class="rcode"><div class="source"><pre class="knitr r">omega = 2*pi/365
thawconst.basis = create.constant.basis(thawbasis$rangeval)

betalist = vector("list", 3)
betalist[[1]] = fd(0, thawconst.basis)
betalist[[2]] = fd(omega^2, thawconst.basis)
betalist[[3]] = fd(0, thawconst.basis)
harmaccelLfd0 = Lfd(3, betalist)
</pre></div>
</div></div>

However, the code used above is a little bit cumbersome. Most of the time the differential operator we need is just a power of $D$ or all the coefficients $\beta_j$ are constants, and in such situation, we can simply use `int2Lfd` and `vec2Lfd` functions.

<div class="chunk" id="ex10-2"><div class="rcode"><div class="source"><pre class="knitr r">## use 'int2Lfd' to create the acceleration object
accelLfd = int2Lfd(2)

## use 'vec2Lfd' to create the harmonic acceleration object and
## the result is the same with 'harmaccelLfd0'
harmaccelLfd = vec2Lfd(c(0,omega^2,0), thawbasis$rangeval)
all.equal(harmaccelLfd0[-1], harmaccelLfd[-1])

## back to the temperature example introduced above
LmeanTempVec = eval.fd(day.5, meanTempfd, harmaccelLfd)
plot(day.5, LmeanTempVec, type="l", cex=1.2,
     xlab="Day", ylab="Harmonic Acceleration")
abline(h=0)
</pre></div>
</div></div>

> A functional data object for the application of a linear differential operator to an existing functional data object is created by the `deriv.fd` function

<div class="chunk" id="ex10-3"><div class="rcode"><div class="source"><pre class="knitr r">D2tempfd = deriv.fd(temp.fd, 2)
Ltempfd  = deriv.fd(temp.fd, harmaccelLfd)
</pre></div>
</div></div>

