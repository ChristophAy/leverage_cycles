simulate_leverage_cycles = function(params, seed, stochastic) {

  T = params[22]
  
  out = matrix(0, nrow = T, ncol = 8)

  run_leverage_cycles_model(out, params, stochastic, seed )

  data = data.frame(out)
  colnames(data) = colnames = c("equity","price","target_leverage", "realized_leverage", "v", "L", "wn", "n")
  data$time = seq(1,T)*0.1

  return(data)
}