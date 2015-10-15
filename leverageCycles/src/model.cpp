#include "model.h"

double uniform() {
	/*
	Uniform random number in [0,1]
	*/

	double r;
    double rmax = RAND_MAX;
    r = rand()/rmax;
    return r;
}

double normal(double var) {
	/*
	Normal random number with mu = 0, var =var
	*/

	double x1,x2,R,w,y1,y2;
    double rmax = RAND_MAX;
    x1 = rand()/rmax;
    x2 = rand()/rmax;
    w = 2*M_PI*x1;
    R = sqrt(-2*log(x2));
    y1 = R*cos(w);
    y2 = R*sin(w);
    
    return (var*y1);
}


Model::Model(std::map<std::string,double> double_parameters,int n_time_steps, int n_record, bool dn, bool ins, 
	bool ts_rec,std::string time_stamp, std::map<std::string,double> alpha_control_parameters, bool alpha_control) {
	/*
	Initialized model given standard set of parameters.
	*/

	Dt = double_parameters["Dt"]; // time step tau
	tvar = double_parameters["tvar"]; // VaR horizon
    Et = double_parameters["Et"]; // equity target
    b = double_parameters["b"];  // b - cyclicality parameter
    l0 = double_parameters["l0"]; // sigma0
    w1 = double_parameters["w1"]; // bank investment weight into risky asset
    rho = double_parameters["rho"]*Dt;// fund reversion to fundamental
    d = double_parameters["d"]*Dt;  // delta - EWMA memory
    a = double_parameters["a"];  //alpha - bank riskiness
    theta = double_parameters["theta"]*Dt; //theta - portfolio adjustment
    eta = double_parameters["eta"]*Dt;// eta - equity redistribution
    mu = double_parameters["mu"]; // fundamental value
    po = double_parameters["po"]; // lagged price
    wn = double_parameters["wn"]; // nitial fund investment weight into risky asset

    // Stochastic parameters

	do_noise = dn; // boolean flag for stochastic simulation
	std  = double_parameters["std"]*sqrt(Dt); // standard dev. of fund noise process
	has_generator = false; // set to false on intialization, set random number generator later
    

    // GARCH parameters

    a0 = double_parameters["a0"]; // baseline volatility
    a1 = double_parameters["a1"]; // noise innovation term
    b1 = double_parameters["b1"]; // autoregressive term

    // Further parameters

    T = n_time_steps; // number of time steps
    start_record = n_record; // start recording at time step
    p_max = double_parameters["pmax"]; // price threshold for instability
    eps = double_parameters["eps"]; // offset from fixed point in initialization
    rs_quantile = double_parameters["rs_quant"]; // quantile for computation of realized shortfall
    initialize_stochastic = ins; // flag for stochastic initialization procedure
    time_series_record = ts_rec; // flag if want to record model state time series

    // alpha_control parameters
    rho_alpha = alpha_control_parameters["rho_alpha"]; // adjustment rate of alpha if controlled
    theta_alpha = alpha_control_parameters["theta_alpha"];
    delta_alpha = alpha_control_parameters["delta_alpha"]; //return estimation horizon 
    control_alpha = alpha_control; 


    std::string garch_indicator;
    if (b1 > 0) {
    	garch_indicator = "garch";
    } else {
    	garch_indicator = "no_garch";
    }

    if (time_series_record == true) {
    	std::ostringstream strs;
	    strs << "_ts" << "_Et_" << Et;
	    // const std::string filename = "output/" + time_stamp + "_" + strs.str()+"_.txt";
	    // const std::string filename = "output/" + time_stamp + strs.str()+ garch_indicator + "_.txt";
	    const std::string filename = time_stamp + strs.str()+ garch_indicator + "_.txt";
	    out.open(filename);	
    }

    if (out.good() == false) {
        printf("%s\n","file corrupt" );
    }

    unstable = false;
    min_lev = true;

    negative_equity = 0;

    time_series_out["equity"] = std::vector<double>();
    time_series_out["p"] = std::vector<double>();
    time_series_out["target_leverage"] = std::vector<double>();
    time_series_out["realized_leverage"] = std::vector<double>();
    time_series_out["v"] = std::vector<double>();
    time_series_out["L"] = std::vector<double>();
    time_series_out["wn"] = std::vector<double>();
    time_series_out["n"] = std::vector<double>();

    use_external_random_number = false;

}

void Model::initialize_state() {

	/*
	Initialize model state close to fixed point with offset eps on price / perceived risk
	*/

	double std_i = 0;
    if (do_noise == true && initialize_stochastic == true) {
        std_i = pow(std,2);
    }
	double lev = a*pow((l0+std_i),b);
    double L = (lev - 1)*Et;
    double n = lev*Et*w1/mu;

	state.p = mu+eps;
	state.po = mu;
	state.wn = wn;
	state.n = n;
	state.L = L;
	state.sigma = 0;
	state.error = 0;

	// alpha control set up
	state.alpha_control = a;
	state.av_return = 0;

	update_additional_state_var();

	state.time = 0;
	unstable = false;

	if (do_noise == true && initialize_stochastic == true) {
		state.v = pow(std,2);
	} else {
		state.v = eps;
	}

	negative_equity = 0;

	// time_series_out["equity"].clear()
 //    time_series_out["p"].clear()
 //    time_series_out["target_leverage"].clear()
 //    time_series_out["realized_leverage"].clear()
 //    time_series_out["v"].clear()
 //    time_series_out["L"].clear()
 //    time_series_out["wn"].clear()
 //    time_series_out["n"].clear()
}

double Model::reset_state(Model& reference_model, double d0, bool rel_initialize) {

	double d1 = 0;

	if (rel_initialize == false) {

		// Reset model state relative to reference state

		// Compute difference to core variables

		d1 = d1 + pow(state.p - reference_model.get_state().p,2);
		d1 = d1 + pow(state.po - reference_model.get_state().po,2);
		d1 = d1 + pow(state.wn - reference_model.get_state().wn,2);
		d1 = d1 + pow(state.n - reference_model.get_state().n,2);
		d1 = d1 + pow(state.L - reference_model.get_state().L,2);
		d1 = d1 + pow(state.v - reference_model.get_state().v,2);
		d1 = sqrt(d1);

		// core variables

		state.p = reference_model.get_state().p +  d0/d1* (state.p - reference_model.get_state().p);
		state.po = reference_model.get_state().po +  d0/d1* (state.po - reference_model.get_state().po);
		state.wn = reference_model.get_state().wn + d0/d1* (state.wn - reference_model.get_state().wn);
		state.n = reference_model.get_state().n + d0/d1* (state.n - reference_model.get_state().n);
		state.L = reference_model.get_state().L + d0/d1*(state.L - reference_model.get_state().L);
		state.v = reference_model.get_state().v + d0/d1*(state.v - reference_model.get_state().v);

		// additional variables (non core)
		state.sigma = reference_model.get_state().sigma;
		state.error = reference_model.get_state().error;
		
		state.alpha_control = reference_model.get_state().alpha_control;
		state.av_return = reference_model.get_state().av_return;

	    
	} else {

		// Initialize relative to reference state
		state.p = reference_model.get_state().p + d0; // offset on price
		state.po = reference_model.get_state().po;
		state.wn = reference_model.get_state().wn;
		state.n = reference_model.get_state().n;
		state.L = reference_model.get_state().L;
		state.v = reference_model.get_state().v;

		d1 = d1 + pow(state.p - reference_model.get_state().p,2);
		d1 = d1 + pow(state.po - reference_model.get_state().po,2);
		d1 = d1 + pow(state.wn - reference_model.get_state().wn,2);
		d1 = d1 + pow(state.n - reference_model.get_state().n,2);
		d1 = d1 + pow(state.L - reference_model.get_state().L,2);
		d1 = d1 + pow(state.v - reference_model.get_state().v,2);
		d1 = sqrt(d1);

		// additional variables (non core)
		state.sigma = reference_model.get_state().sigma;
		state.error = reference_model.get_state().error;
		
		state.alpha_control = reference_model.get_state().alpha_control;
		state.av_return = reference_model.get_state().av_return;
	}
	
	update_additional_state_var();

	// state.target_leverage = reference_model.get_state().target_leverage;
	// state.old_equity = reference_model.get_state().old_equity;
	// state.equity = reference_model.get_state().equity;
	// state.realized_leverage = reference_model.get_state().realized_leverage;
	// state.bank_fund_rel_size = reference_model.get_state().bank_fund_rel_size;
    return d1;
}

void Model::map() {
	/*
	Updates the state of the model from t to t+1
	*/

	double p_t = state.p;
	double po_t = state.po;
	double wn_t = state.wn;
	double n_t = state.n;
	double L_t = state.L;
	double v_t = state.v;
	double sigma_t = state.sigma;
	double error_t = state.error;

	double A = n_t*p_t/w1; // corresponding assets
    double E = A-L_t; // corresponding equity
    double dE = (Et - E)*eta; // equity adjustment

    double lev_target;// target leverage

    if (control_alpha == true) {
    	// double  m = .0005;
    	// double da = 0;
    	// if (state.av_return > m) {
    	// 	da = rho_alpha*(a - state.alpha_control) + theta_alpha*(state.av_return-m);	
    	// } else if (state.av_return < m) {
    	// 	da = rho_alpha*(a - state.alpha_control) + theta_alpha*(state.av_return+m);	
    	// }
    	// double da = rho_alpha*((mu + .1) - p_t);
    	double da = rho_alpha*(a - state.alpha_control) + theta_alpha*state.av_return;
    	state.alpha_control = state.alpha_control + da;
    	state.av_return = (1-delta_alpha)*state.av_return + delta_alpha*log(p_t/po_t);
    	lev_target = state.alpha_control*pow((v_t + l0),b);
    } else {
    	lev_target = a*pow((v_t + l0),b);
    }
    double c1 = (1 - w1)*n_t*p_t/w1 + dE; // risk free asset of the bank after equity adjustment
    double c2 = (1 - wn_t)*(1 - n_t)*p_t/wn_t - dE; //risk free asset of the fund after equity adjustment

    double dA;
    if (min_lev == true) {
    	dA = std::max((lev_target*E - A)*theta,-L_t); // balance sheet adjustment of bank to match target leverage	
    } else {
    	dA = (lev_target*E - A)*theta; // balance sheet adjustment of bank to match target leverage
    }
    

    double dw; // change of fund risky asset investment weight

    if (do_noise == true) {
    	state.sigma = a0 + a1*pow(error_t,2) + b1*sigma_t;
    	// state.error = sqrt(state.sigma)*normal(1);
    	state.error = sqrt(state.sigma)*normal_random();
    	dw = wn_t/p_t *( (mu-p_t)*rho + state.error );
    } else {
    	dw = wn_t/p_t *( (mu-p_t)*rho );
    }

    state.wn = std::min(std::max(wn_t + dw,0.1),0.9); // fund investment weight
    state.p = pow(1 - w1*n_t - (1 - n_t)*state.wn,-1)*(w1*(c1 + dA) + state.wn *c2);// price
    state.n = w1*(n_t*state.p + c1 + dA)/state.p; // ownership
    state.po = p_t; // lagged price
    state.v = (1-d)*v_t + d*pow(log(p_t/po_t)*tvar/Dt,2); // perceived risk
    state.L = L_t + dA; // liabilities

    update_additional_state_var();
    state.time = state.time +1;

}

void Model::record() {
	/*
	Records model state time series and further derived quantities.
	*/

	// if (time_series_record == true) {
	// 	out << state.p << "," << state.v << "," << state.L << ", "<<  state.wn  << ", " << state.n << ", " << 
	// 	state.target_leverage << ", " << state.realized_leverage << ", " << state.equity << "\n";

	// }

	time_series_out["equity"].push_back(state.equity);
    time_series_out["p"].push_back(state.p);
    time_series_out["target_leverage"].push_back(state.target_leverage);
    time_series_out["realized_leverage"].push_back(state.realized_leverage);
    time_series_out["v"].push_back(state.v);
    time_series_out["L"].push_back(state.L);
    time_series_out["wn"].push_back(state.wn);
    time_series_out["n"].push_back(state.n);
	
	if (state.time > start_record) {
		// record negative asset returns
		if (state.p < state.po) {
			n_asset_returns.push_back(log(state.p/state.po));
			// check if equity > 0
			if (state.equity > 0 && state.old_equity > 0) {
				n_equity_returns.push_back(log(state.equity/state.old_equity));	
			} else{
				// do nothing for now. 
				// TODO: Implement flag for when equity drops below zero.
			}
			
		}

		if (state.equity < 0) {
			negative_equity += state.equity/T;
		}

		target_leverage_t.push_back(state.target_leverage);
		realized_leverage_t.push_back(state.realized_leverage);
		relative_size_t.push_back(state.bank_fund_rel_size);
		asset_returns.push_back(log(state.p/state.po));
		price_t.push_back(state.p);
	}
	

}

void Model::print_state() {
	std::cout << state.p << ", " << state.wn << ", " << state.n << ", " << state.v << ", " << state.L << "\n";
}

void Model::update_additional_state_var() {
	/*
	Computes quantities derived from model state
	*/
	if (control_alpha == true) {
		state.target_leverage = state.alpha_control*pow((state.v + l0),b);
	} else {
		state.target_leverage = a*pow((state.v + l0),b);
	}
	state.old_equity = state.equity;
	state.equity = state.n*state.p/w1 - state.L;
	state.realized_leverage = state.n*state.p/w1/state.equity;
	state.bank_fund_rel_size = state.n*state.p/w1/((1-state.n)*state.p/state.wn);

}

void Model::comput_realized_shortfall() {

	/*
	Compute realized shortfall from rs_quantile worst negative returns on assets and equity
	*/

	std::sort(n_asset_returns.begin(),n_asset_returns.end());
	std::sort(n_equity_returns.begin(),n_equity_returns.end());
	asset_RS = 0;
	equity_RS = 0;
	int count = 0;

	// First only look at assets. Need to split loops since vectors may not be the same size.
	if (n_asset_returns.size() > 0) {
		for (int i = 0; i < n_asset_returns.size()*rs_quantile + 0.1; i ++) {
			asset_RS += n_asset_returns[i];
			count++;
		}
		asset_RS = asset_RS/count;
	} else {
		// flag that no negative asset returns found in this run.
		asset_RS = 1.;
	}

	// Then look at equity returns.
	count = 0;
	if (n_equity_returns.size() > 0) {
		for (int i = 0; i < n_equity_returns.size()*rs_quantile + 0.1; i ++) {
			equity_RS += n_equity_returns[i];
			count++;
		}
		equity_RS = equity_RS/count;
	} else {
		// flag that no negative asset returns found in this run.
		equity_RS = 1.;
	}
}

void Model::compute_timeseries_averages() {
	double sum;
	sum = std::accumulate(realized_leverage_t.begin(), realized_leverage_t.end(), 0.0);
	average_realized_leverage = sum / realized_leverage_t.size();
	sum = std::accumulate(target_leverage_t.begin(), target_leverage_t.end(), 0.0);
	average_target_leverage = sum / target_leverage_t.size();
	sum = std::accumulate(relative_size_t.begin(), relative_size_t.end(), 0.0);
	average_relative_size = sum / relative_size_t.size();

	sum = 0;
	double sum2 = 0;
	int count = 0;

	for (int i = 0; i < asset_returns.size(); i++) {
		sum+= asset_returns[i];
		sum2+= pow(asset_returns[i],2);
		count++;
	}

	asset_volatility =  0;

	if (count > 0) {
		double Er = sum/count;
		double Er2 = sum2/count;
		asset_volatility = Er2 - pow(Er,2);
	}
}

void Model::run() {

	/*
	Produce one time series of model.
	*/

	initialize_state();
	bool do_run = true;
	while (state.time < T && do_run == true) {
		map();
		record();
		if (state.p > p_max || state.p < 0) {
			do_run = false;
			unstable = true;
			if (time_series_record == true) printf("Unstable: breaking here. \n");
		}
	}
	compute_timeseries_averages();
	comput_realized_shortfall();
	
	if (time_series_record == true) {
		printf("average realized_leverage %f , average target leverage %f, average relative size %f , RS  %f, vol %f \n", 
		average_realized_leverage, average_target_leverage, average_relative_size, asset_RS, asset_volatility);
		out.close();	
	}

	wipe_all_timeseries_records();
	
}

void Model::set_alpha(double alpha) {
	a = alpha;
}

void Model::set_b(double bee) {
	b = bee;
}

bool Model::get_unstable(){
	return unstable;
}

double Model::get_asset_RS() {
	return asset_RS;
}

double Model::get_equity_RS() {
	return equity_RS;
}

double Model::get_average_target_leverage() {
	return average_target_leverage;
}

double Model::get_average_realized_leverage() {
	return average_realized_leverage;
}

double Model::get_average_relative_size() {
	return average_relative_size;
}

double Model::get_volatility() {
	return asset_volatility;
}

void Model::wipe_all_timeseries_records() {

	n_asset_returns.clear(); // negative asset returns
    n_equity_returns.clear(); // negative equity returns
    target_leverage_t.clear(); // target leverage time series
    realized_leverage_t.clear(); // realized leverage time series
    relative_size_t.clear();
    asset_returns.clear();
    price_t.clear();
}

bool Model::is_GARCH() {
	/*
	Returns true if model is running GARCH
	*/
	bool isg = false;
	if (b1 != 0 && a1 != 0) {
		isg = true;
	}
	return isg;
}

void Model::set_alpha_theta(double tee) {
	theta_alpha = tee;
}

void Model::set_alpha_delta(double dee) {
	delta_alpha = dee;
}

void Model::set_d(double dee) {
	d = dee;
}

// void Model::set_random_generator(std::default_random_engine gen) {
// 	generator = gen;
// 	has_generator = true;
// }

double Model::normal_random() {
	double num;
	double mean = 0;
	double s_dev = 1;
	if (use_external_random_number == false) {
		// if (has_generator == true) {
		// std::normal_distribution<double> distribution(mean,s_dev);
		// num = distribution(generator);
		// } else {
		// 	num = normal(1);
		// }
		num = normal(1);
	} else {
		num = external_random_number;
	}

	
	return num;
}

double Model::uniform_random() {
	double num;
	// if (has_generator == true) {
 //  		std::uniform_real_distribution<double> distribution(0.0,1.0);
 //  		num = distribution(generator);
	// } else {
	// 	num = uniform();
	// }
	num = uniform();

	return num;
}

void Model::set_Et(double eetee) {
	Et = eetee;
}

void Model::set_GARCH_params(double aa0, double aa1, double bb1) {
	a0 = aa0;
	a1 = aa1;
	b1 = bb1;
}

void Model::set_all_double_parameters(std::map<std::string,double> double_parameters) {
	
	Dt = double_parameters["Dt"]; // time step tau
	tvar = double_parameters["tvar"]; // VaR horizon
    Et = double_parameters["Et"]; // equity target
    b = double_parameters["b"];  // b - cyclicality parameter
    l0 = double_parameters["l0"]; // sigma0
    w1 = double_parameters["w1"]; // bank investment weight into risky asset
    rho = double_parameters["rho"]*Dt;// fund reversion to fundamental
    d = double_parameters["d"]*Dt;  // delta - EWMA memory
    a = double_parameters["a"];  //alpha - bank riskiness
    theta = double_parameters["theta"]*Dt; //theta - portfolio adjustment
    eta = double_parameters["eta"]*Dt;// eta - equity redistribution
    mu = double_parameters["mu"]; // fundamental value
    po = double_parameters["po"]; // lagged price
    wn = double_parameters["wn"]; // nitial fund investment weight into risky asset

    // Stochastic parameters
	std  = double_parameters["std"]*sqrt(Dt); // standard dev. of fund noise process
    
    // GARCH parameters

    a0 = double_parameters["a0"]; // baseline volatility
    a1 = double_parameters["a1"]; // noise innovation term
    b1 = double_parameters["b1"]; // autoregressive term
}

State Model::get_state() {
	return state;
}

void Model::set_use_external_random_number(bool use) {
	use_external_random_number = use;
}

void Model::set_external_random_number(double num) {
	external_random_number = num;
}

void Model::set_min_lev(bool ml) {
	min_lev = ml;
}

double Model::get_n_equity() {
	return negative_equity;
}

std::map<std::string,std::vector<double> > Model::get_time_series_out() {
	return time_series_out;
}
