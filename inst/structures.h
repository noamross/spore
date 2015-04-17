#ifndef STRUCTURES
#define STRUCTURES

#include <RcppArmadillo.h>

struct process_ {
  int max_i;
	double lambda;
	double lambda_ex;
	double alpha;
	double mu;
	double r;
	double d;
	double K;
	};

struct init_ {
	double N;
	double P;
	double A;
}


struct sim_ {
  double macro_step_size;
  int n_relaxation_steps;
  int n_simulation_steps;
  int n_bursts;
  double time_max
)
};

struct parms_ {
	process_ process;
	sim_ sim;
};


struct burst_ {
   vec of micro_state
   
   method: clear
   method: repopulate from macro_state
};

struct macro_state_ {
	double N;
	double P;
	double A;
	double dN;
	double dP;
	double dA;
	double ddN;
	double ddP;
	double ddA;
	double time;
	
	method: write to file
	method: clear and repopulate from burst_ restriction
	
};

struct micro_state_ {
	vec pop
	double time
	datalog_ log
};


struct datalog_ {
};

#endif /* STRUCTURES */