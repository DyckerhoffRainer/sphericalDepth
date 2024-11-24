#pragma once
#ifndef __SignedMeasureDepthR__
#define __SignedMeasureDepthR__

namespace DataDepth {

	/*******************************************************************************************************************************/
	enum AlgType { comb, rec };
	enum AlgSubtype { standard, adaptive, single, multiple, nGP, GP };

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
	void aHD(const double x[], const val_t val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, val_t result[], const AlgType& alg = comb, const int& target = 2, const AlgSubtype& algSub = standard, const int& nThreads = 0);
	template<typename val_t>
	void aHD(const double x[],                    int& n, const int& d, const double z[], int& m, const int ind[], int& l, val_t result[], const AlgType& alg = comb, const int& target = 2, const AlgSubtype& algSub = standard, const int& nThreads = 0);

	extern "C" void aHD_int_val(const double x[], const int val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, int result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads) {
		aHD(x, val, n, d, z, m, ind, l, result, alg, target, algSub, nThreads);
	};
	extern "C" void aHD_double_val(const double x[], const double val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, double result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads) {
		aHD(x, val, n, d, z, m, ind, l, result, alg, target, algSub, nThreads);
	};

	extern "C" void aHD_int(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, int result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads) {
		aHD(x, n, d, z, m, ind, l, result, alg, target, algSub, nThreads);
	};
	extern "C" void aHD_double(const double x[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, double result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads) {
		aHD(x, n, d, z, m, ind, l, result, alg, target, algSub, nThreads);
	};


}

#endif
