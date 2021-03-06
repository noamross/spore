
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <functional>
#include <boost/unordered_map.hpp>

using namespace Rcpp;



// [[Rcpp::depends(RcppArmadillo)]]
//' @export
// [[Rcpp::export]]
List micro_state_c_step(const IntegerVector micro_state, const List parms, const double control, const double time) {
  NumericVector micro_state_out = as<NumericVector>(wrap(micro_state));
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
    (lambda * infections + lambda_ex * exp(-control)) * N,
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
  //Rcout << infections << std::endl;
  //Rcout << as<arma::vec>(rates) << std::endl;
  //Rcout << rate << ", " << rates[0] << ", " << rates[1] << ", " << rates[2] << ", " << rates[3] << ", " << event[0] << std::endl;



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
  return List::create(_["micro_state"] = as<IntegerVector>(wrap(micro_state_out)), _["time_next"] = time_next);
}

//' @export
// [[Rcpp::export]]
IntegerVector micro_state_c_stepto(const IntegerVector micro_state, const List parms, const double control, const double time, const double timeto, const int record = 0, const int run = 1) {

  double time0 = time;
  // RECORD = !is.null(record) && ("connection" %in% class(parms$micro_record))

  // 	if(RECORD) {
  // 		micro_state_c.record(micro_state, control, time, connection = record, run = run, start = time0)
  // 	}

  List out;
  double time_next;
  IntegerVector micro_state_out = clone(micro_state);
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


// // No check for missing values
// // [[Rcpp::export]]
// NumericVector tabulate1(NumericVector x, const int max) {
//   NumericVector counts(max + 1);
//
//   int n = x.size();
//   for (int i = 0; i < n; i++) {
//     int pos = x[i];
//     if (pos < max && pos >= 0) counts[pos]++;
//   }
//
//   return counts;
// }
//
// //' @export
// // [[Rcpp::export]]
// NumericVector lift_macro_state(const NumericVector macro_state, const List parms) {
//   NumericVector vals;
//   vals = rpois(macro_state[0], macro_state[1]/macro_state[0]);
//   NumericVector micro_state = tabulate1(vals, parms["max_i"]);
// 	return micro_state;
// }


//' @export
// [[Rcpp::export]]
IntegerVector tabulate2(IntegerVector x, const int max) {
  IntegerVector counts(max + 1);

  int n = x.size();
  for (int i = 0; i < n; i++) {
    int pos = x[i];
    if (pos < max && pos >= 0) counts[pos]++;
  }

  return counts;
}

//' @export
// [[Rcpp::export]]
IntegerVector lift_macro_state(const IntegerVector macro_state, const List parms) {
  IntegerVector N = seq_len(macro_state[0]) - 1;
  int P = macro_state[1];
  int max_i = parms["max_i"];
  NumericVector prob = NumericVector::create();
  IntegerVector buckets = RcppArmadillo::sample(N, P, true, prob);
  IntegerVector bucket_counts = tabulate2(buckets, macro_state[0]);
  IntegerVector out = tabulate2(bucket_counts, max_i);
  out[0] = out[0] - 1;
  //if(out[1] == 287) {
    return out;
//  } else{
 //   return buckets;
//  }
}

//' @export
// [[Rcpp::export]]
IntegerVector restrict_micro_state(IntegerVector micro_state) {
  IntegerVector macro_state(2);
  for(int i = 0; i < micro_state.size(); i++) {
    macro_state[0] += micro_state[i];
    macro_state[1] += micro_state[i] * i;
  }
  return macro_state;
}

//' @export
// [[Rcpp::export]]
List macro_state_c_step(const IntegerVector macro_state, const List parms, const double control, const double time) {
   IntegerVector micro_state = lift_macro_state(macro_state, parms);
   List micro_state_out = micro_state_c_step(micro_state, parms, control, time);
   IntegerVector macro_state_out = restrict_micro_state(micro_state_out["micro_state"]);
   return List::create(_["macro_state"] = macro_state_out, _["time_next"] = micro_state_out["time_next"]);
}

//' @export
// [[Rcpp::export]]
IntegerVector macro_state_c_stepto(const IntegerVector macro_state, const List parms, const double control, const double time, const double timeto) {
   IntegerVector micro_state = lift_macro_state(macro_state, parms);
   IntegerVector micro_state_out = micro_state_c_stepto(micro_state, parms, control, time, timeto);
   IntegerVector macro_state_out = restrict_micro_state(micro_state_out);
   return macro_state_out;
}

//' @export
// [[Rcpp::export]]
List macro_state_c_step_diff(const IntegerVector macro_state, const List parms, const double control, const double time) {
  IntegerVector micro_state = lift_macro_state(macro_state, parms);
  List micro_state_out = micro_state_c_step(micro_state, parms, control, time);
  IntegerVector macro_state_out_diff = restrict_micro_state(micro_state_out["micro_state"]) - macro_state;
  double dtime = micro_state_out["time_next"];
  dtime -= time;
  return List::create(_["macro_state_out_diff"] = macro_state_out_diff, _["dtime"] = dtime);
}

//' @export
// [[Rcpp::export]]
NumericVector macro_state_c_step_aves(const IntegerVector macro_state, const List parms, const double control, const double time) {
  IntegerVector diff = IntegerVector::create(0, 0);
  IntegerVector diffs = IntegerVector::create(0, 0);
  double times = 0;
  NumericVector dtime = NumericVector::create(0);
  List runout;
  int reps = parms["n_sims"];
  for (int i = 0; i < reps; i++) {
    runout = macro_state_c_step_diff(macro_state, parms, control, time);
    diff = runout["macro_state_out_diff"];
    diffs = diffs + diff;
    dtime = runout["dtime"];
    times = times + dtime[0];
  }
  NumericVector out2(2);
  out2[0] = double (diffs[0]) /times;
  out2[1] = double (diffs[1]) /times;
  times = 0; diffs[0] = 0; diffs[1] = 0; diff[0] = 0; diff[1] = 0; dtime[0] = 0;
  return out2;
}



