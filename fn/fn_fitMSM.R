fitMSM <- function(obs_with_state, emI=rbind(c(0, 0), c(0.1, 0)), fp = NULL, cv = NULL){
  QmI <- rbind(c(0, 0.1), c(0.1, 0)) # Initial Transition proba
  # emI - Initial Emission proba - Rows: Hidden States, Cols: Observations
  msmFit <- msm(state~normday, subject=codenum, data=obs_with_state, qmatrix=QmI, ematrix=emI,
                est.initprobs=T, opt.method="bobyqa", covariates = cv, control=list(maxfun=10000, iprint=3),
                fixedpars = fp)
  statetables <- statetable.msm(state, codenum, obs_with_state)
  paths <- viterbi.msm(msmFit)
  return(list(model = msmFit, table = statetables, vPath = paths))
}