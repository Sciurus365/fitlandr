#include <Rcpp.h>
using namespace Rcpp;
using namespace std;

//' @export
// [[Rcpp::export]]
double fast_bilinear_cpp(NumericVector x, NumericVector y, NumericMatrix z,
                         double x0, double y0){
	int x_ind = floor((x0 - x(0)) / (x(1) - x(0)));
	int y_ind = floor((y0 - y(0)) / (y(1) - y(0)));
	if (x_ind < 0 || x_ind > x.length() - 2 || y_ind < 0 || y_ind > y.length() - 2) {
		return 0;
	}

	double z11 = z(x_ind, y_ind);
	double z12 = z(x_ind, y_ind + 1);
	double z21 = z(x_ind + 1, y_ind);
	double z22 = z(x_ind + 1, y_ind + 1);
	double ex = (x0 - x(x_ind)) / (x(x_ind + 1) - x(x_ind));
	double ey = (y0 - y(y_ind)) / (y(y_ind + 1) - y(y_ind));
	double result = (1 - ex) * (1 - ey) * z11 + (1 - ex) * (ey) * z12 + (ex) * (1 - ey) * z21 + (ex) * (ey) * z22;

	return result;
}
