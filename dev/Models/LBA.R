
rLBA <- function(lR,pars,p_types=c("v","sv","b","A","t0"),posdrift = TRUE,
                 ok=rep(TRUE,length(lR)))
  # lR is an empty latent response factor lR with one level for each accumulator.
  # pars is a matrix of corresponding parameter values named as in p_types
  # pars must be sorted so accumulators and parameter for each trial are in
  # contiguous rows.
{
  bad <- rep(NA, length(lR)/length(levels(lR)))
  out <- data.frame(R = bad, rt = bad)
  nr <- length(levels(lR))
  dt <- matrix(Inf,nrow=nr,ncol=nrow(pars)/nr)
  t0 <- pars[,"t0"]
  pars <- pars[ok,]
  if (!all(p_types %in% dimnames(pars)[[2]]))
    stop("pars must have columns ",paste(p_types,collapse = " "))
  dt[ok] <- (pars[,"b"]-pars[,"A"]*runif(dim(pars)[1]))/
    msm::rtnorm(dim(pars)[1],pars[,"v"],pars[,"sv"],ifelse(posdrift,0,-Inf))
  dt[dt<0] <- Inf
  bad <- apply(dt,2,function(x){all(is.infinite(x))})
  R <- apply(dt,2,which.min)
  pick <- cbind(R,1:dim(dt)[2]) # Matrix to pick winner
  # Any t0 difference with lR due to response production time (no effect on race)
  rt <- matrix(t0,nrow=nr)[pick] + dt[pick]
  R <- factor(levels(lR)[R],levels=levels(lR))
  R[bad] <- NA
  rt[bad] <- Inf
  ok <- matrix(ok,nrow=length(levels(lR)))[1,]
  out$R[ok] <- levels(lR)[R][ok]
  out$R <- factor(out$R,levels=levels(lR))
  out$rt[ok] <- rt[ok]
  out
}
#### Model functions ----