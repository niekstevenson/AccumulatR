# Run from the repository root with the exported, unchanged posterior cases.
# Benchmark the complete likelihood; contexts and parameter matrices are cached.
library(AccumulatR)
source('dev/scripts/bench_cumulative.R')
source('dev/scripts/bench_bisset_reference.R')
folder <- 'dev/scripts/scratch_outputs/bisset'
cases <- readRDS(file.path(folder,'posterior.rds'))
results <- benchmark_cumulative(cases,'866512a',repetitions=25)
saveRDS(results,file.path(folder,'results.rds'))
seconds <- t(vapply(results,function(x)x$seconds[c('baseline','current')],numeric(2)))
group <- sub('/.*','',names(cases))
timing <- aggregate(seconds,list(model=group),sum)
timing$ratio <- timing$current/timing$baseline
write.table(timing,file.path(folder,'timing.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
print(timing,row.names=FALSE)

accuracy <- list()
for(name in names(cases)) {
  case <- cases[[name]]
  data <- case$data
  trials <- data[!duplicated(data$trials),]
  omissions <- which(trials$cond=='ignore' & is.na(trials$R))
  omissions <- omissions[!duplicated(trials$SSD[omissions])]
  alternative <- startsWith(name,'S2_without_S1')
  # Check all saved original-model draws; one final draw per alternative fit.
  draws <- if(alternative) 3L else seq_along(case$parameters)
  if(alternative && length(omissions)) omissions <- omissions[which.max(trials$SSD[omissions])]
  for(j in draws) for(i in omissions) {
    rows <- which(data$trials==trials$trials[i])
    reference <- bisset_reference(case$parameters[[j]][rows,,drop=FALSE],trials$SSD[i],
                                  s2_without_s1=alternative)
    accuracy[[length(accuracy)+1L]] <- data.frame(
      model=name,draw=j,SSD=trials$SSD[i],reference=reference,
      old=exp(results[[name]]$results$baseline$loglik[i,j]),
      current=exp(results[[name]]$results$current$loglik[i,j]))
  }
}
accuracy <- do.call(rbind,accuracy)
saveRDS(accuracy,file.path(folder,'accuracy.rds'))
print(aggregate(cbind(old=abs(accuracy$old-accuracy$reference),
                      current=abs(accuracy$current-accuracy$reference)),
                list(model=sub('/.*','',accuracy$model)),max),row.names=FALSE)
