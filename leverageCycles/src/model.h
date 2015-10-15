#ifndef __MODEL_H_INCLUDED__
#define __MODEL_H_INCLUDED__

#include <vector>
#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <string>
#include <algorithm>
#include <math.h>
#include <sstream>
#include <set>
#include <numeric>
#include "model.h"
#include <time.h>
// #include <random>

struct State {
	/*
	Holds all relevant state variable for the dynamical system and some additional derived quantities we want to keep track of.
	*/

	// Variables of dynamical system

	double p; // price of risky asset
	double po; // last price of risky asset
	double v; // perceived risk of risky asset
	double n; // bank ownership of risky asset
	double wn; // fund risky asset investment weight
	double L; // bank liabilities
	double sigma; // GARCH volatility
	double error; // GARCH error

	// additional things to keep track off (not part of dynamical system)

	double target_leverage;
	double realized_leverage;
	double equity;
	double old_equity;
	double bank_fund_rel_size;

	double alpha_control; // current value of alpha if controlled by policy maker
	double av_return; // return estimate

	int time; // model time step
};

class Model {
	// Holds all parameters of simulation and method to run dynamical system.

	// Model parameters (deterministic)

	double Dt; // time step tau
	double tvar; // VaR horizon
    double Et; // equity target
    double b;  // b - cyclicality parameter
    double l0; // sigma0
    double w1; // bank investment weight into risky asset
    double rho;// fund reversion to fundamental
    double d;  // delta - EWMA memory
    double a;  //alpha - bank riskiness
    double theta; //theta - portfolio adjustment
    double eta;// eta - equity redistribution
    double mu; // fundamental value
    double po; // lagged price
    double wn; // nitial fund investment weight into risky asset

    // Parameters for control of alpha 

    double rho_alpha; // adjustment rate of alpha if controlled
    double delta_alpha; //return estimation horizon 
    double theta_alpha; // alpha adjustment rate
    bool control_alpha; // flag if want to control alpha to stabilize prices


    // Stochastic parameters

	bool do_noise; // boolean flag for stochastic simulation
	double std; // standard dev. of fund noise process

    // Random number generator
    // std::default_random_engine generator;
    bool has_generator;
    
    // GARCH parameters

    double a0; // baseline volatility
    double a1; // error innovation term
    double b1; // autoregressive term

    // Further parameters

    int T; // number of time steps
    int start_record; // time for transients to die out
    double p_max; // price threshold for instability
    double eps; // offset from fixed point in initialization
    bool initialize_stochastic; // flag for stochastic initialization procedure
    double rs_quantile; // quantile for computation of realized shortfall
    bool time_series_record; // flag if want to record time series of model state

    // State of the model

    State state;

    // Recording
    std::ofstream out; // file for output of model state time series

    std::vector<double> n_asset_returns; // negative asset returns
    std::vector<double> n_equity_returns; // negative equity returns
    std::vector<double> target_leverage_t; // target leverage time series
    std::vector<double> realized_leverage_t; // realized leverage time series
    std::vector<double> relative_size_t; // realitive size bank to fund time series
    std::vector<double> asset_returns; // asset returns
    std::vector<double> price_t; // price time series

    double asset_RS; // realized shortfall for risky asset
    double equity_RS; // realized shortfall for bank equity
    bool unstable; // flag if model is in unstable regime
    double average_target_leverage;
    double average_realized_leverage;
    double average_relative_size;
    double asset_volatility; // volatility of risky asset.
    double negative_equity; // how far does equity drop below zero

    bool use_external_random_number; // indicates if external stream of random numbers should be used
    double external_random_number; // external random number to be supplied
    bool min_lev;

    std::map<std::string,std::vector<double> > time_series_out; // holds time series output of model



    // Constructors and functions

public:
	Model(std::map<std::string,double>,int, int, bool, bool,bool,std::string, std::map<std::string,double>, bool ); // initializes the model
	void run(); // run model for T time steps
	void wipe_all_timeseries_records(); // clears all storage vectors
	void set_alpha(double); // set alpha
	void set_b(double); // set b
	bool get_unstable(); // check if model is unstable
	double get_asset_RS();
	double get_equity_RS();
	double get_average_target_leverage();
	double get_average_realized_leverage();
	double get_average_relative_size();
	double get_volatility();
	bool is_GARCH();
    void set_alpha_theta(double); // set adjustment rate for alpha control
    void set_alpha_delta(double); // set estimation horizon for alpha control
    void set_d(double); // set delta (risk estimation horizon)
    // void set_random_generator(std::default_random_engine); // set random number generator
    double normal_random(); // generate normall random number
    double uniform_random(); // generate uniform random number
    void set_GARCH_params(double,double,double); // set parameters of GARCH noise
    void set_Et(double); // set target equity
    void set_all_double_parameters(std::map<std::string,double>); // resets all double parameters
    State get_state(); // returns state of the model
    double reset_state(Model&, double, bool); // returns state of the model
    void map(); // update state of model
    void initialize_state(); // initialize state close to fixed point
    void record(); // record state of the model
    void comput_realized_shortfall(); // compute realized shortfall of asset and equity returns
    void compute_timeseries_averages(); // compute time series averages
    
    void set_use_external_random_number(bool); // sets external random number flag
    void set_external_random_number(double); // provides external random number
    void set_min_lev(bool);
    double get_n_equity();
    std::map<std::string,std::vector<double> > get_time_series_out();


private:
	// void map(); // update state of model
	// void initialize_state(); // initialize state close to fixed point
	// void record(); // record state of the model
	void print_state(); // print selection of state variables
	void update_additional_state_var(); // update derived state variables
	// void comput_realized_shortfall(); // compute realized shortfall of asset and equity returns
	// void compute_timeseries_averages(); // compute time series averages


};


#endif // __MODEL_H_INCLUDED__ 