#include <cmath>

double norm2(const double x[], int d) {
	double result = 0;
	for (int i = 0; i < d; i++) result += x[i] * x[i];
	return sqrt(result);
}

double InnerProduct(const double x[], const double y[], int d) {
	double sum = 0;
	for (int i = 0; i < d; i++) sum += x[i] * y[i];
	return sum;
}

void Normalize(double x[], int d) {
	double Norm = 0;
	for (int i = 0; i < d; i++) Norm += x[i] * x[i];
	Norm = sqrt(Norm);
	for (int i = 0; i < d; i++) x[i] /= Norm;
}
