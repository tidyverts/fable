#include <Rcpp.h>
#include <vector>

// Filters a linear innovations ("single source of error") state space model:
//   y[t] = w'x[t-1] + e[t]
//   x[t] = F x[t-1] + g e[t]
// This one recursion underlies fitting, forecasting and simulation for
// TBATS - all of the model-specific structure lives in how F, g, w and x0
// are constructed in R (see tbats_ssm() in R/tbatsmodel.R). Unlike a general
// Kalman filter, g does not vary over time, so no matrix inversion is needed
// at each step, only a matrix-vector product - but with observations
// numbering in the thousands, doing that step-by-step loop in compiled code
// rather than R avoids the overhead of `n` separate R-level matrix multiplies.
// [[Rcpp::export]]
Rcpp::List tbats_filter(const Rcpp::NumericVector &y,
                        const Rcpp::NumericMatrix &F,
                        const Rcpp::NumericVector &g,
                        const Rcpp::NumericVector &w,
                        const Rcpp::NumericVector &x0) {
  const int n = y.size();
  const int m = x0.size();

  Rcpp::NumericVector fitted(n);
  Rcpp::NumericVector e(n);
  Rcpp::NumericMatrix states(n, m);

  std::vector<double> prev(x0.begin(), x0.end());
  std::vector<double> cur(m);

  for (int t = 0; t < n; t++) {
    double yhat = 0.0;
    for (int i = 0; i < m; i++) yhat += w[i] * prev[i];
    const double err = y[t] - yhat;

    for (int i = 0; i < m; i++) {
      double fx = 0.0;
      for (int j = 0; j < m; j++) fx += F(i, j) * prev[j];
      cur[i] = fx + g[i] * err;
    }

    fitted[t] = yhat;
    e[t] = err;
    for (int i = 0; i < m; i++) states(t, i) = cur[i];
    prev.swap(cur);
  }

  return Rcpp::List::create(
    Rcpp::Named("fitted") = fitted,
    Rcpp::Named("e") = e,
    Rcpp::Named("states") = states
  );
}
