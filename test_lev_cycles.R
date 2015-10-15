library(ggplot2)
library(leverageCycles)

get_double_params = function() {

	params = vector()

	params[1] = 1 # time step tau
	params[2] = 1 # tvar
	params[3] = 2.26689 # equity target
	params[4] = -0.5 # b - cyclicality parameter
	params[5] = 1e-06 # sigma0
	params[6] = 0.3 # bank investment weight into risky asset
	params[7] = 0.01 # fund reversion to fundamental
	params[8] = 0.05 # delta - EWMA memory
	params[9] = 0.01#75 # alpha - bank riskiness
	params[10] = 0.95 # theta - portfolio adjustment
	params[11] = 1 # eta - equity redistribution
	params[12] = 25 # fundamental value
	params[13] = 25 # lagged price
	params[14] = 0.5 # initial fund investment weight into risky asset
	params[15] = 0.1 # standard dev. of fund noise process
	params[16] = 0.001 # a0
	params[17] = 0.0158902 # a1
	params[18] = 0.874195 # b1
	params[19] = 1000 # pmax
	params[20] = 1e-10 # eps
	params[21] = 0.05 # rs_quant

	params[22] = 500

	return(params)
}

params = get_double_params()
stochastic = 1
seed = 1

data = simulate_leverage_cycles(params,seed,stochastic)
p = ggplot(data, aes(x = time, y = price)) + geom_line()
print(p)