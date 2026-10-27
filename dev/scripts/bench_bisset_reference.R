# Independent one-dimensional quadrature for the original Bisset response rule.
# Parameters are AccumulatR's natural-scale runtime rows GO, S1, IS, S2.
bisset_reference <- local({
  jacobi <- matrix(0,32,32)
  off <- (1:31)/sqrt(4*(1:31)^2-1)
  jacobi[cbind(1:31,2:32)] <- jacobi[cbind(2:32,1:31)] <- off
  rule <- eigen(jacobi,symmetric=TRUE)
  nodes <- (rule$values+1)/2
  weights <- rule$vectors[1,]^2
  exgauss <- function(p) {
    mu <- p['p1']; sigma <- p['p2']; tau <- p['p3']
    weighted_tail <- function(z) {
      w <- z-sigma/tau
      result <- exp((sigma/tau)*(.5*sigma/tau-z)+pnorm(w,log.p=TRUE))
      far <- w < -10
      u <- 1/w[far]^2
      log_mills <- -.5*log(2*pi)-log(-w[far])+
        log(1-u+3*u^2-15*u^3+105*u^4-945*u^5+10395*u^6)
      result[far] <- exp(-z[far]^2/2+log_mills)
      result
    }
    tail <- function(x) weighted_tail((x-mu)/sigma)
    norm <- pnorm(mu/sigma)+tail(0)
    list(pdf=function(t) ifelse(t>0,tail(t)/tau/norm,0),
         cdf=function(t) ifelse(t>0,pmax(0,1-(pnorm((mu-t)/sigma)+tail(t))/norm),0),
         survival=function(t) ifelse(t>0,(pnorm((mu-t)/sigma)+tail(t))/norm,1),
         points=pmax(0,c(mu+c(-8,-3,0,3,8)*sigma,
                         mu+c(1,5,15,30)*tau)),
         expectation=function(fn) {
           scale <- sigma+tau
           integral(function(u) fn(mu+scale*u)*scale/tau*
               weighted_tail(u*(scale/sigma))/norm,
             c(c(-8,-3,0,3,8)*sigma/scale,c(1,5,15,30)*tau/scale),
             lower=-mu/scale)
         })
  }
  integral <- function(fn,points,upper=Inf,lower=0) {
    breaks <- sort(unique(c(lower,points[is.finite(points)&points>lower&points<upper],upper)))
    sum(vapply(seq_len(length(breaks)-1L),function(i)
      integrate(fn,breaks[i],breaks[i+1L],rel.tol=1e-8,abs.tol=1e-10,
                subdivisions=500L)$value,numeric(1)))
  }
  function(parameters,ssd,rt=NULL,s2_without_s1=FALSE) {
    p<-parameters[1,]; v<-p['p1']/p['p4']; d<-(p['p2']+nodes*p['p3'])/p['p4']
    go <- function(t,cdf=FALSE) {
      x <- pmax(t-p['t0'],1e-100)
      z <- outer(x*v,d,'-')/sqrt(x)
      values <- if(cdf) pnorm(z)+exp(outer(rep(1,length(x)),2*v*d)+
          pnorm(-outer(x*v,d,'+')/sqrt(x),log.p=TRUE)) else
        exp(-z*z/2)*rep(d,each=length(x))/sqrt(2*pi*x^3)
      drop(values %*% weights)*(t>p['t0'])
    }
    s1<-exgauss(parameters[2,]); is<-exgauss(parameters[3,]); s2<-exgauss(parameters[4,])
    ps1<-s1$expectation(function(t)1-go(t+ssd,TRUE))
    ps2<-s2$expectation(function(t)is$survival(t))
    q<-parameters[2,'q']
    omission<-(1-q)*ps1*ps2
    if(s2_without_s1) omission <- (1-q)*s2$expectation(function(t) {
      no_early_go <- vapply(t,function(u) {
        end <- u+ssd
        1-go(end,TRUE)+integral(function(g)go(g)*s1$cdf(g-ssd),
          c(ssd+s1$points,p['t0']+d/v),end)
      },numeric(1))
      is$survival(t)*no_early_go
    })
    if(is.null(rt)) return(omission)
    density<-vapply(rt,function(t) {
      x<-t-ssd
      a<-if(x>0) integral(function(u)is$pdf(u)*s2$survival(u),c(is$points,s2$points),x) else 0
      b<-integral(function(u)go(u)*s1$cdf(u-ssd),
                  c(ssd+s1$points,p['t0']+d/v),t)
      q*go(t)+(1-q)*(go(t)*(s1$survival(x)+s1$cdf(x)*a)+is$pdf(x)*s2$survival(x)*b)
    },numeric(1))
    list(omission=omission,density=density)
  }
})
