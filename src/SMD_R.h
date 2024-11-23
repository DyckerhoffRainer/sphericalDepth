#pragma once
#ifndef __SignedMeasureDepthR__
#define __SignedMeasureDepthR__

namespace DataDepth {

	/*******************************************************************************************************************************/
	enum Algorithm { standard, adaptive, single, multiple, nGP, GP };

	// signed measure depth in dimension d = 1, i.e., on the line
	template<typename val_t>
	void smD_d1s(const double x[], const val_t valX[], const int& n, const double z[], const int valZ[], const int& m,
		const int ind[], const int& l, val_t result[]);
	// signed measure depth in dimension d = 2, i.e., in the plane
	template<typename val_t>
	void smD_d2(const double x[], const val_t valX[], const int& n, const double z[], const int valZ[], const int& m,
		const int ind[], const int& l, val_t result[]);

	// angular halfspace depth in dimension d = 1, i.e., on the 0-sphere { -1, 1 }
	//void aHD_d1(const double x[], int& n, const double z[], int& m, const int ind[], int& l, int result[]);
	// angular halfspace depth in dimension d = 2, i.e., on the 1-sphere
	//void aHD_d2(const double x[], int& n, const double z[], int& m, const int ind[], int& l, int result[], const Algorithm& Alg = adaptive);
	// angular halfspace depth in dimension d = 3, i.e., on the 2-sphere
	//void aHD_d3(const double x[], int& n, const double z[], int& m, const int ind[], int& l, int result[], const Algorithm& Alg = nGP2);

	// angular halfspace depth in arbitrary dimension, recursive algorithm
	//void aHD_Rec(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, int result[], const Algorithm& Alg = standard);
	// angular halfspace depth in arbitrary dimension, combinatorial algorithm, reduction to d == target
	//void aHD_Comb(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, int result[], const Algorithm& Alg = standard);

	// angular halfspace depth in arbitrary dimension, recursive algorithm
	//extern "C" void aHD_R(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, int result[], const Algorithm& Alg = standard);

	template<typename val_t>
	void aHD_R(const double x[], const val_t val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, val_t result[], const Algorithm& Alg = standard, const int& nThreads = 0);
	template<typename val_t>
	void aHD_R(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, val_t result[], const Algorithm& Alg = standard, const int& nThreads = 0);

	extern "C" void aHD_R_int_val(const double x[], const int val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, int result[], const Algorithm& alg, const int& nThreads) {
		aHD_R(x, val, n, d, z, m, ind, l, target, result, alg, nThreads);
	};
	extern "C" void aHD_R_double_val(const double x[], const double val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, double result[], const Algorithm& alg, const int& nThreads) {
		aHD_R(x, val, n, d, z, m, ind, l, target, result, alg, nThreads);
	};

	extern "C" void aHD_R_int(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, int result[], const Algorithm& alg, const int& nThreads) {
		aHD_R(x, n, d, z, m, ind, l, target, result, alg, nThreads);
	};
	extern "C" void aHD_R_double(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, double result[], const Algorithm& alg, const int& nThreads) {
		aHD_R(x, n, d, z, m, ind, l, target, result, alg, nThreads);
	};



	template<typename val_t>
	void aHD_C(const double x[], const val_t val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, val_t result[], const Algorithm& Alg = standard, const int& nThreads = 0);
	template<typename val_t>
	void aHD_C(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, val_t result[], const Algorithm& Alg = standard, const int& nThreads = 0);

	extern "C" void aHD_C_int_val(const double x[], const int val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, int result[], const Algorithm& alg, const int& nThreads) {
		aHD_C(x, val, n, d, z, m, ind, l, target, result, alg, nThreads);
	};
	extern "C" void aHD_C_double_val(const double x[], const double val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, double result[], const Algorithm& alg, const int& nThreads) {
		aHD_C(x, val, n, d, z, m, ind, l, target, result, alg, nThreads);
	};

	extern "C" void aHD_C_int(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, int result[], const Algorithm& alg, const int& nThreads) {
		aHD_C(x, n, d, z, m, ind, l, target, result, alg, nThreads);
	};
	extern "C" void aHD_C_double(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, const int& target, double result[], const Algorithm& alg, const int& nThreads) {
		aHD_C(x, n, d, z, m, ind, l, target, result, alg, nThreads);
	};


}

#endif
