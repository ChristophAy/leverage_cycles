#include "model.h"
#include <Rcpp.h>

using namespace Rcpp;

void parse_R_parameters(std::map<std::string,double>& double_parameters, int& n_time_steps, bool& do_noise, NumericVector double_params, int stochastic) {

    double_parameters["Dt"] = double_params[0]; // time step tau
    double_parameters["tvar"] = double_params[1]; // VaR horizon
    double_parameters["Et"] = double_params[2]; // equity target
    double_parameters["b"] = double_params[3];  // b - cyclicality parameter
    double_parameters["l0"] = double_params[4]; // sigma0
    double_parameters["w1"] = double_params[5]; // bank investment weight into risky asset
    double_parameters["rho"] = double_params[6];// fund reversion to fundamental
    double_parameters["d"] = double_params[7];  // delta - EWMA memory
    double_parameters["a"] = double_params[8];  //alpha - bank riskiness
    double_parameters["theta"] = double_params[9]; //theta - portfolio adjustment
    double_parameters["eta"] = double_params[10];// eta - equity redistribution
    double_parameters["mu"] = double_params[11]; // fundamental value
    double_parameters["po"] = double_params[12]; // lagged price
    double_parameters["wn"] = double_params[13]; // initial fund investment weight into risky asset
    double_parameters["std"] = double_params[14]; // standard dev. of fund noise process
    double_parameters["a0"] =   double_params[15];//pow(double_parameters["std"],2); 
    double_parameters["a1"] = double_params[16];//0.07;
    double_parameters["b1"] = double_params[17];//0.9;
    double_parameters["pmax"] = double_params[18];
    double_parameters["eps"] = double_params[19];
    double_parameters["rs_quant"] = double_params[20];

    n_time_steps = int(double_params[21]);

    if (stochastic == 1) {
        do_noise = true;
    } else {
        do_noise = false;
    }

}


// [[Rcpp::export]]
NumericMatrix run_leverage_cycles_model(NumericMatrix out, 
    NumericVector double_params, int stochastic, int seed) {

    srand (seed);

	std::map<std::string,double> double_parameters;
	std::map<std::string,double> search_parameters;
	std::map<std::string,double> sweep_parameters;
	std::map<std::string,double> alpha_control_parameters;
	std::map<std::string,double> delta_sweep_parameters;

	int n_time_steps;
    int n_record;
    bool do_noise;
    bool initialize_stochastic;
    bool ts_rec;

    // Default parameters

    double_parameters["Dt"] = 1.; // time step tau
	double_parameters["tvar"] = 1.; // VaR horizon
    double_parameters["Et"] = 1; // equity target
    double_parameters["b"] = -0.5;  // b - cyclicality parameter
    double_parameters["l0"] = pow(10,-6); // sigma0
    double_parameters["w1"] = 0.3; // bank investment weight into risky asset
    double_parameters["rho"] = 0.01;// fund reversion to fundamental
    double_parameters["d"] = 0.05;  // delta - EWMA memory
    double_parameters["a"] = 0.02;  //alpha - bank riskiness
    double_parameters["theta"] = 0.95; //theta - portfolio adjustment
    double_parameters["eta"] = 1.;// eta - equity redistribution
    double_parameters["mu"] = 25; // fundamental value
    double_parameters["po"] = 25; // lagged price
    double_parameters["wn"] = 0.5; // initial fund investment weight into risky asset
    double_parameters["std"] = 0.1; // standard dev. of fund noise process
    double_parameters["a0"] = 	0.1;//pow(double_parameters["std"],2); 
    double_parameters["a1"] = 0;//0.07;
    double_parameters["b1"] = 0;//0.9;
    double_parameters["pmax"] = pow(10,3);
    double_parameters["eps"] = pow(10,-10);
    double_parameters["rs_quant"] = 0.05;
    
    n_time_steps = 5000;
    n_record = 100;
    do_noise = true;
    initialize_stochastic = true;
    ts_rec = true;

    
    search_parameters["a_min"] = pow(10,-2);
    search_parameters["a_max"] = pow(10,3);
    search_parameters["rel_tol"] = pow(10,-3);
    search_parameters["Nmax"] = 5000;
    search_parameters["Nruns_search"] = 10;

    
    sweep_parameters["Na"] = 100;
    sweep_parameters["b_min"] = -0.5;
    sweep_parameters["b_max"] = 0.5;
    sweep_parameters["Nb"] = 100;
    sweep_parameters["Nruns_sweep"] = 30;

    // alpha control parameters

    bool alpha_control = false;
    alpha_control_parameters["rho_alpha"] = 0.5;
    alpha_control_parameters["delta_alpha"] = 0.2;
    alpha_control_parameters["delta_alpha_min"] = pow(10,-2);
    alpha_control_parameters["delta_alpha_max"] = 0.7;
    alpha_control_parameters["theta_alpha"] = 1.0;
    alpha_control_parameters["theta_alpha_min"] = pow(10,-2);
    alpha_control_parameters["theta_alpha_max"] = 7.0;
    alpha_control_parameters["N_theta"] = 100;
    alpha_control_parameters["N_delta"] = 100;
    alpha_control_parameters["Nruns_alpha_sweep"] = 30;

    // delta sweep parameters
    delta_sweep_parameters["Na_delta"] = 100;
    delta_sweep_parameters["delta_min"] = 0.001;
    delta_sweep_parameters["delta_max"] = 0.5;
    delta_sweep_parameters["N_delta_sweep"] = 100;
    delta_sweep_parameters["Nruns_delta_sweep"] = 30;


    std::string time_stamp = "";
    ts_rec = false;

    parse_R_parameters(double_parameters, n_time_steps, do_noise, double_params, stochastic);
    
	Model model (double_parameters,n_time_steps,n_record, do_noise,initialize_stochastic,
		ts_rec,time_stamp, alpha_control_parameters, alpha_control);

    model.set_min_lev(false);

    model.run();

    std::map<std::string,std::vector<double> > time_series = model.get_time_series_out();
    std::vector<std::string> output_order;
    output_order.push_back("equity");
    output_order.push_back("p");
    output_order.push_back("target_leverage");
    output_order.push_back("realized_leverage");
    output_order.push_back("v");
    output_order.push_back("L");
    output_order.push_back("wn");
    output_order.push_back("n");

    for (int i = 0; i < out.nrow(); i++) {
        for (int j = 0; j < out.ncol(); j++) {
            out(i,j) = time_series[output_order[j]][i];
        }
    }

    return out;
}