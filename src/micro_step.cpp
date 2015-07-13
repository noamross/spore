
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <boost/unordered_map.hpp>

using namespace Rcpp;



// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
List micro_state_c_step(const NumericVector micro_state, const List parms, const double control, const double time) {
  NumericVector micro_state_out = clone(micro_state);
  //int max_i = parms["max_i"];
  NumericVector is = parms["is"];
  int infections = sum(is * micro_state_out);
  int N = sum(micro_state_out);
  if (N == 0) {
    return List::create(_["micro_state"] = micro_state_out, _["time_next"] = R_PosInf);
  }
  double lambda = parms["lambda"];
  double mu = parms["mu"];
  double lambda_ex = parms["lambda_ex"];
  double alpha = parms["alpha"];
  //double alpha_power = parms["alpha_power"];
  double r = parms["r"];
  double K = parms["K"];
  double d = parms["d"];
  int i;
  NumericVector irate;
  IntegerVector samp;
  NumericVector rates = NumericVector::create(
    lambda * infections + lambda_ex * exp(-control) * N,
    mu * infections,
    alpha * infections + (d * N),
    r * N * (1 - N/K)
  );

  if(rates[3] < 0) {
    rates[3] = 0;
  }

  double rate = sum(rates);
  double time_next = time + Rf_rexp(1 / rate);
  IntegerVector options = IntegerVector::create(0, 1, 2, 3);
  IntegerVector event = RcppArmadillo::sample(options, 1, false, rates);
  //Rcout << time << ", " << rates[0] << ", " << rates[1] << ", " << rates[2] << ", " << rates[3] << ", " << event[0] << " ," << N << ", " << infections << std::endl;


  switch(event[0]) {
  case 0:
    irate = micro_state_out;
    samp = RcppArmadillo::sample(is, 1, false, irate);
    i = samp[0];
    micro_state_out[i] = micro_state_out[i] - 1;
    if (i != micro_state_out.size() + 1) {
      micro_state_out[i + 1] = micro_state_out[i + 1] + 1;
    }
    break;
  case 1:
    irate =  mu * is * micro_state_out;
    //     irate = irate / std::accumulate(irate.begin(),
    //                                     irate.end(), 0.0);
    samp = RcppArmadillo::sample(is, 1, false, irate);
    i = samp[0];
    micro_state_out[i] = micro_state_out[i] - 1;
    micro_state_out[i - 1] = micro_state_out[i - 1] + 1;
    break;
  case 2:
    irate = alpha * is * micro_state_out + d*micro_state_out;
    //     irate = irate / std::accumulate(irate.begin(),
    //                                     irate.end(), 0.0);
    samp = RcppArmadillo::sample(is, 1, false, irate);
    i = samp[0];
    micro_state_out[i] = micro_state_out[i] - 1;
    break;
  case 3:
    micro_state_out[0] = micro_state_out[0] + 1;
    break;
  }
  //  return List::create(_["samp"] = samp, _["i"] = i, _["event"] = event[0],_["rates"] =  rates, _["micro_state"] = micro_state, _["micro_state_out"] = micro_state_out, _["time_next"] = time_next);
//  Rcout <<
  return List::create(_["micro_state"] = micro_state_out, _["time_next"] = time_next);
}

//' @export
// [[Rcpp::export]]
NumericVector micro_state_c_stepto(const NumericVector micro_state, const List parms, const double control, const double time, const double timeto, const int record = 0, const int run = 1) {

  double time0 = time;
  // RECORD = !is.null(record) && ("connection" %in% class(parms$micro_record))

  // 	if(RECORD) {
  // 		micro_state_c.record(micro_state, control, time, connection = record, run = run, start = time0)
  // 	}

  List out;
  double time_next;
  NumericVector micro_state_out = clone(micro_state);
  while(time0 < timeto) {
    out = micro_state_c_step(micro_state_out, parms, control, time0);
    time_next = out["time_next"];
    if(time_next < timeto) {
      micro_state_out = out["micro_state"];
    }
    time0 = time_next;

      // 		if(RECORD) {
      // 			micro_state_c.record(micro_state, control, time, connection = record, run = run, start = time0)
      // 		}
  }
//  Rcout << micro_state_out << std::endl;
  return micro_state_out;
}


// No check for missing values
// [[Rcpp::export]]
NumericVector tabulate1(NumericVector x, const int max) {
  NumericVector counts(max + 1);

  int n = x.size();
  for (int i = 0; i < n; i++) {
    int pos = x[i];
    if (pos < max && pos >= 0) counts[pos]++;
  }

  return counts;
}

//' @export
// [[Rcpp::export]]
NumericVector lift_macro_state(const NumericVector macro_state, const List parms) {
  NumericVector vals;
  vals = rpois(macro_state[0], macro_state[1]/macro_state[0]);
  NumericVector micro_state = tabulate1(vals, parms["max_i"]);
	return micro_state;
}



