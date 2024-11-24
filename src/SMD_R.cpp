#include <memory>
#include <cmath>
#include <limits>
#include <algorithm>
#include <set>
#include <functional>
#include <iostream>
#include <random>
#include <omp.h>
#include <execution>
#include "auxLinAlg.h"
#include "Matrix.h"
#include "HD.h"
#include "SMD_R.h"

using namespace std;
using namespace placeholders;
using namespace dyMatrixClass;


namespace DataDepth {

	int NumberOfCoresToUse{};

	const double eps = 1e-13; // 13
	const double eps_angle = 1e-13;  // 13
	const double eps_point = 1e-13;  // 13
	const double eps_point_atan = 1e-12;  // 12
	const double eps_HH = 1e-13;  // 13
	const int PrecBinary = 39; // 1.82 * 10^(-12)
	constexpr double pi = 3.14159265358979323846;

	const int limitSingle = 30;

	template<typename val_t>
	struct pointInfo {
		const double* x;
		val_t val;
	};

	template<typename val_t>
	struct pointInfoP {
		const double* x;
		val_t val;
		int index;
	};

	struct angleInfo {
		double val;
		int p1;
		int p2;
	};

	class intSet {
	protected:
		int n;
		unique_ptr<int[]> elems;
		unique_ptr<bool[]> isElem;
	public:
		intSet(int size) : n{ 0 }, elems{ new int[size] {} }, isElem{ new bool[size] {} } {};
		int operator[](int i) const { return elems[i]; }
		int size() const { return n; }
		void insert(int elem) {
			if (!isElem[elem]) {
				isElem[elem] = true;
				elems[n++] = elem;
			}
		}
		bool isElement(int elem) { return isElem[elem]; }
		void clear() {
			for (int i = 0; i < n; i++) isElem[elems[i]] = false;
			n = 0;
		}
	};

	int nCr(int n, int k) {
		int result{ n };
		for (int i = 2; i <= k; i++) result = result * (--n) / i;
		return result;
	}

	int FindX(int r, int m) {
		const long long fak[13] = { 1,1,2,6,24,120,720,5040,40320,362880,3628800,39916800, 479001600 };
		constexpr double pi = 3.14159265358979323846;
		if (r == 0) return m;
		if (m == 1) return r + 1;
		if (m == 2) return (int)ceil((1 + sqrt(1 + 8 * r)) / 2);
		if (m <= 12) return (int)ceil(pow(r * fak[m], 1.0 / m) + (m - 1) / 2.0);
		return (int)ceil(m * exp(log(r) / m + log(2 * pi * m) / (2 * m) + 1.0 / (12 * m * m) - 1.0 / (360 * m * m * m * m) - 1) + (m - 1) / 2.0);
	}

	void unRank(int r, int k, int ans[]) {
		for (int i = k; i > 0; --i) {
			int x = FindX(r, i);
			while (r >= nCr(x, i)) x++;
			ans[i - 1] = x - 1;
			r -= nCr(x - 1, i);
		}
	} // Element()

	template <typename T> int sgn(T val) {
		return (T(0) < val) - (val < T(0));
	}

	inline double sqr(double x) { return x * x; }

	inline double RoundBinary(double x) { return ldexp(round(ldexp(x, PrecBinary)), -PrecBinary); }

	// Several routines to sort angles or points

	bool cmpAngles(const angleInfo& a, const angleInfo& b) {
		if (a.val < b.val - eps_angle) return true;
		if (a.val < b.val + eps_angle) {
			if (a.p1 < b.p1) return true;
			if ((a.p1 == b.p1) && (a.p2 < b.p2)) return true;
		}
		return false;
	}

	bool compAngles(const angleInfo& a, const angleInfo& b, const int ranks[]) {
		if (a.val < b.val - eps_angle) return true;
		if (a.val < b.val + eps_angle) {
			if (ranks[a.p1] < ranks[b.p1]) return true;
			if ((ranks[a.p1] == ranks[b.p1]) && (ranks[a.p2] < ranks[b.p2])) return true;
		}
		return false;
	}

	bool compPAngles(const angleInfo* a, const angleInfo* b, const int ranks[]) {
		if (a->val < b->val - eps_angle) return true;
		if (a->val < b->val + eps_angle) {
			if (ranks[a->p1] < ranks[b->p1]) return true;
			if ((ranks[a->p1] == ranks[b->p1]) && (ranks[a->p2] < ranks[b->p2])) return true;
		}
		return false;
	}

	template<typename val_t>
	bool cmpPointsEpsAtan(const pointInfo<val_t>& a, const pointInfo<val_t>& b) {
		if (atan(a.x[0]) < atan(b.x[0]) - eps_point) return true;
		if ((atan(a.x[0]) <= atan(b.x[0]) + eps_point) && (atan(a.x[1]) > atan(b.x[1]) + eps_point)) return true;
		return false;
	}

	template<typename val_t>
	bool cmpPPointsEpsAtan(const pointInfo<val_t>* a, const pointInfo<val_t>* b) {
		if (atan(a->x[0]) < atan(b->x[0]) - eps_point_atan) return true;
		if ((atan(a->x[0]) <= atan(b->x[0]) + eps_point_atan) && (atan(a->x[1]) > atan(b->x[1]) + eps_point_atan)) return true;
		return false;
	}

	template<typename val_t>
	bool ltPointsRounded(const pointInfoP<val_t>& a, const pointInfoP<val_t>& b, int d) {
		bool result{ false };
		for (int i = 0; i < d; i++) {
			double x{ RoundBinary(a.x[i]) }, y{ RoundBinary(b.x[i]) };
			if (x < y) return true;
			if (x > y) break;
		}
		return false;
	}

	template<typename val_t>
	bool eqPointsRounded(const pointInfoP<val_t>& a, const pointInfoP<val_t>& b, int d) {
		for (int i = 0; i < d; i++)
			if (RoundBinary(a.x[i]) != RoundBinary(b.x[i])) return false;
		return true;
	}

	bool eqPointsOnSphere(const double a[], const double b[]) {
		double x[3], y[3];
		x[2] = 1.0 / sqrt(1.0 + sqr(a[0]) + sqr(a[1]));
		x[1] = a[1] * x[2];
		x[0] = a[0] * x[2];
		y[2] = 1.0 / sqrt(1.0 + sqr(b[0]) + sqr(b[1]));
		y[1] = b[1] * y[2];
		y[0] = b[0] * y[2];
		return sqr(x[0] - y[0]) + sqr(x[1] - y[1]) + sqr(x[2] - y[2]) < sqr(eps_point);
	}

	bool eqPoints(const double a[], const double b[]) {
		return (fabs(a[0] - b[0]) < eps_point) && (fabs(a[1] - b[1]) < eps_point);
	}

	template<typename val_t>
	bool cmpRealsEps(const pointInfo<val_t>& a, const pointInfo<val_t>& b) {
		return (a.x[0] < b.x[0] - eps_point);
	}

	/**************************************************************************************************************/

	// Several auxiliary routines

	bool pointsOnEquator(const double x[], int n) {
		double eps_equator = min(1e-2, 1 - pow(0.5, 1.0 / n));
		// guarantees that under a uniform distribution on the sphere the probability of getting no points
		// near the equator is at least 1/2
		for (int i = 0; i < n; i++)
			if (!isnan(x[3 * i]) && (fabs(x[3 * i + 2]) < eps_equator)) return true;
		return false;
	}

	void Householder(double z[], const double u[]) {
		// Hz is set to H * z
		// H is the Householder matrix that maps e1 to u (and vice versa)
		// z and u are assumed to have norm 1
		if (!isnan(z[0]) && !isnan(u[0])) {
			int d{ 3 };
			double lambda{ 0.0 }, ip;
			ip = InnerProduct(z, u, d);
			if (fabs(u[0]) < 1 - eps_HH)
				lambda = (u[0] < 0) ? (ip - z[0]) / (1 - u[0]) : (ip + z[0]) / (1 + u[0]);
			for (int i = 0; i < d; i++) z[i] -= lambda * u[i];
			z[0] += (u[0] < 0) ? lambda : -lambda;
		}
		else z[0] = nan("");
	}

	void HouseholderProjection(const double z[], const double u[], int d, double Hz[], bool& pos) {
		// Hz is set to proj(H * z)
		// H is the Householder matrix that maps e1 to u (and vice versa)
		// proj is the projection on the last d-1 coordinates
		// z and u are assumed to have norm 1
		// Finally Hz is normalized to have norm 1. If || Hz ||^2 <= 2*eps, Hz is set to 0.
		if (!isnan(z[0]) && !isnan(u[0])) {
			double lambda{ 0.0 }, ip;
			ip = InnerProduct(z, u, d);
			if (fabs(u[0]) < 1 - eps_HH)
				lambda = (u[0] < 0) ? (ip - z[0]) / (1 - u[0]) : (ip + z[0]) / (1 + u[0]);
			double Norm{ 0.0 };
			for (int i = 1; i < d; i++) {
				Hz[i - 1] = z[i] - lambda * u[i];
				Norm += Hz[i - 1] * Hz[i - 1];
			}
			if (Norm > 2 * eps_HH) {
				Norm = sqrt(Norm);
				for (int i = 0; i < d - 1; i++) Hz[i] /= Norm;
			}
			else {
				Hz[0] = nan("");
				pos = (ip > 0);
			}

		}
		else Hz[0] = nan("");
	}

	void project(const cMatrix& Q, const double x[], int n, double px[]) {
		int d{ Q.Rows() }, pd{ Q.Cols() };
		double norm;

		for (int i = 0; i < n; i++) {
			norm = 0;
			for (int j = 0; j < pd; j++) {
				px[i * pd + j] = 0;
				for (int k = 0; k < d; k++) px[i * pd + j] += x[i * d + k] * Q(k, j);
				norm += px[i * pd + j] * px[i * pd + j];
			}
			norm = sqrt(norm);
			if (norm > 2 * eps_HH) {
				for (int j = 0; j < pd; j++) px[i * pd + j] /= norm;
			}
			else px[i * pd] = nan("");
		}
	}

	template<typename val_t>
	void project2(const cMatrix& Q, const double x[], const val_t val[], int n, int dComp, double px[], val_t& HD) {
		int d{ Q.Rows() }, dSpan{ d - dComp };
		unique_ptr<int[]> ind{ new int[n] };
		int ny{ 0 };
		double norm;

		for (int i = 0; i < n; i++) {
			norm = 0;
			for (int j = 0; j < dComp; j++) {
				px[i * dComp + j] = 0;
				for (int k = 0; k < d; k++) px[i * dComp + j] += x[i * d + k] * Q(k, dSpan + j);
				norm += px[i * dComp + j] * px[i * dComp + j];
			}
			norm = sqrt(norm);
			if (norm > 2 * eps_HH) {
				for (int j = 0; j < dComp; j++) px[i * dComp + j] /= norm;
			}
			else {
				px[i * dComp] = nan("");
				ind[ny++] = i;
			}
		}
		HD = 0;
		if (ny > dSpan) {
			unique_ptr<double[]> y{ new double[ny * dSpan]{} };
			unique_ptr<val_t[]> valY{ new val_t[ny] };
			for (int i = 0; i < ny; i++) {
				for (int j = 0; j < dSpan; j++) {
					for (int k = 0; k < d; k++) y[i * dSpan + j] += x[ind[i] * d + k] * Q(k, j);
				}
				valY[i] = val[ind[i]];
			}
			HD = nHD_Comb2(y.get(), valY.get(), ny, dSpan);
		}
	}

	/***********************************************************************************************/

	// The following routines cope with the case that points are not in general position
	// The routines 'adjustForMultiplePoints' and "adjustForMultipleReals' cope with points
	//   appearing multiple times in the data points
	// The routines assume that the array p[] is sorted in ascending order

	template<typename val_t>
	void adjustForMultiplePointsAtan(pointInfo<val_t> p[], int& n) {
		int i = 0, j = 0, k = 0;
		while (k < n) {
			k = j + 1;
			while ((k < n) && eqPointsOnSphere(p[j].x, p[k].x)) k++;
			p[i] = p[j];
			for (j++; j < k; j++) p[i].val += p[j].val;
			i++;
		}
		n = i;
	}

	template<typename val_t>
	void adjustForMultipleReals(pointInfo<val_t> p[], int& n) {
		int i = 0, j = 0, k = 0;
		while (k < n) {
			k = j + 1;
			while ((k < n) && p[k].x[0] <= p[j].x[0] + eps_point) k++;
			p[i] = p[j];
			for (j++; j < k; j++) p[i].val += p[j].val;
			i++;
		}
		n = i;
	}

	template<typename val_t>
	void adjustForMultipleReals(pointInfo<val_t> p[], int& n, const double x[], int rank[]) {
		int i = 0, j = 0, k = 0;
		while (k < n) {
			k = j + 1;
			while ((k < n) && p[k].x[0] <= p[j].x[0] + eps_point) k++;
			p[i] = p[j];
			rank[distance(x, p[j].x)] = i;
			for (j++; j < k; j++) {
				p[i].val += p[j].val;
				rank[distance(x, p[j].x)] = i;
			}
			i++;
		}
		n = i;
	}

	/***************************************************************************************************************************/

	// Signed halfspace depth of z_1,...z_m and x_i1,...x_il w.r.t. x_1,...x_n in R^1, i.e., on the real line

	template<typename val_t>
	void smD_d1m(const double x[], const val_t valX[], const int& nn,
		const double z[], const int valZ[], const int& m,
		const int ind[], const int& l, val_t result[]) {

		int n{ 0 }, lz, uz;
		val_t sumV{ 0 };

		unique_ptr<pointInfo<val_t>[]> p{ new pointInfo<val_t>[nn] };
		for (int i = 0; i < nn; i++) {
			if (!isnan(x[i])) {
				p[n].x = x + i;
				p[n].val = valX[i];
				n++;
			}
		}

		unique_ptr<int[]> rank;
		sort(p.get(), p.get() + n, [](pointInfo<val_t>& a, pointInfo<val_t>& b) { return (a.x[0] < b.x[0]); });
		/*
		if ((NumberOfCoresToUse <= 1) || omp_in_parallel())
			sort(execution::seq, p.get(), p.get() + n, [](pointInfo<val_t>& a, pointInfo<val_t>& b) { return (a.x[0] < b.x[0]); } );
		else
			sort(execution::par, p.get(), p.get() + n, [](pointInfo<val_t>& a, pointInfo<val_t>& b) { return (a.x[0] < b.x[0]); });
		*/
		if (l == 0)
			adjustForMultipleReals(p.get(), n);
		else {
			rank.reset(new int[nn]);
			adjustForMultipleReals(p.get(), n, x, rank.get());
		}

		unique_ptr<val_t[]> cumsumU{ new val_t[n + 1]{0} };
		for (int i = 0; i < n; i++) cumsumU[i + 1] = sumV += p[i].val;

		unique_ptr<val_t[]> minU{ new val_t[n + 1] };
		unique_ptr<val_t[]> maxU{ new val_t[n + 1] };
		unique_ptr<val_t[]> minL{ new val_t[n + 1] };
		unique_ptr<val_t[]> maxL{ new val_t[n + 1] };
		maxL[0] = minL[0] = 0;
		minU[n] = maxU[n] = sumV;
		for (int i = 1;   i <= n; i++) maxL[i] = max(maxL[i - 1], cumsumU[i]);
		for (int i = n-1; i >= 0; i--) minU[i] = min(minU[i + 1], cumsumU[i]);
		for (int i = 1;   i <= n; i++) minL[i] = min(minL[i - 1], cumsumU[i]);
		for (int i = n-1; i >= 0; i--) maxU[i] = max(maxU[i + 1], cumsumU[i]);

		for (int k = 0; k < m; k++) {
			if (!isnan(z[k])) {
				uz = distance(p.get(), lower_bound(p.get(), p.get() + n, pointInfo<val_t>{z + k, 0}, cmpRealsEps<val_t>));
				if ((uz >= n) || (p[uz].x[0] > z[k] + eps)) lz = uz; else lz = uz + 1;
				result[k] = (valZ[k] == 1) ? min(minU[lz], sumV - maxL[uz]) : min(minL[uz] - sumV, -maxU[lz]);
			}
			else result[k] = numeric_limits<val_t>::max();
		}
		if (l > 0) {
			for (int k = 0; k < l; k++) {
				if (!isnan(x[ind[k]])) {
					lz = 1 + (uz = rank[ind[k]]);
					result[m + k] = (valX[ind[k]] > 0) ? min(minU[lz], sumV - maxL[uz]) : min(minL[uz] - sumV, -maxU[lz]); // replaced valZ[ind[k]] == 1 by > 0
				}
				else result[m + k] = numeric_limits<val_t>::max();
			}
		}
	}

	// Signed halfspace depth of z_1,...z_m and x_i1,...x_il w.r.t. x_1,...x_n in R^1, i.e., on the real line

	template<typename val_t>
	void smD_d1s(const double x[], const val_t valX[], const int& nn, const double z[], const int valZ[], const int& m,
		const int ind[], const int& l, val_t result[]) {

		int n{ 0 }, lz, rz;
		val_t sumV{ 0 };

		unique_ptr<pointInfo<val_t>[]> p{ new pointInfo<val_t>[nn] };
		for (int i = 0; i < nn; i++) {
			if (!isnan(x[i])) {
				p[n].x = x + i;
				p[n].val = valX[i];
				n++;
			}
		}
		sort(p.get(), p.get() + n, [](pointInfo<val_t>& a, pointInfo<val_t>& b) { return (a.x[0] < b.x[0]); });
		/*
		if ((NumberOfCoresToUse <= 1) || omp_in_parallel())
			sort(execution::seq, p.get(), p.get() + n, [](pointInfo<val_t>& a, pointInfo<val_t>& b) { return (a.x[0] < b.x[0]); });
		else
			sort(execution::par, p.get(), p.get() + n, [](pointInfo<val_t>& a, pointInfo<val_t>& b) { return (a.x[0] < b.x[0]); });
		*/
		adjustForMultipleReals(p.get(), n);

		unique_ptr<val_t[]> cumsumU{ new val_t[n + 1]{0} };
		for (int i = 0; i < n; i++) cumsumU[i + 1] = sumV += p[i].val;

		for (int k = 0; k < m; k++) {
			if (!isnan(z[k])) {
				lz = 0;
				while ((lz < n) && (p[lz].x[0] < z[k] - eps)) lz++;
				if ((lz < n) && (p[lz].x[0] > z[k] + eps)) rz = lz; else rz = lz + 1;

				if (valZ[k] == 1) {
					val_t minU = sumV, maxU = 0;
					for (int i = 1; i < lz + 1; i++)
						if (cumsumU[i] > maxU) maxU = cumsumU[i];
					for (int i = rz; i <= n; i++)
						if (cumsumU[i] < minU) minU = cumsumU[i];
					result[k] = min(minU, sumV - maxU);
				}
				else {
					val_t minU = 0, maxU = sumV;
					for (int i = 1; i < lz + 1; i++)
						if (cumsumU[i] < minU) minU = cumsumU[i];
					for (int i = rz; i <= n; i++)
						if (cumsumU[i] > maxU) maxU = cumsumU[i];
					result[k] = min(minU - sumV, -maxU);
				}
			}
			else result[k] = numeric_limits<val_t>::max();
		}
		for (int k = 0; k < l; k++) {
			if (!isnan(x[ind[k]])) {
				lz = 0;
				while ((lz < n) && p[lz].x[0] < x[ind[k]] - eps) lz++;
				lz++;

				if (valX[ind[k]] > 0) {
					val_t minU = sumV, maxU = 0;
					for (int i = 1; i < lz; i++)
						if (cumsumU[i] > maxU) maxU = cumsumU[i];
					for (int i = lz; i <= n; i++)
						if (cumsumU[i] < minU) minU = cumsumU[i];
					result[m + k] = min(minU, sumV - maxU);
				}
				else {
					val_t minU = 0, maxU = sumV;
					for (int i = 1; i < lz; i++) //minU = min(minU, cumsumU[i]);
						if (cumsumU[i] < minU) minU = cumsumU[i];
					for (int i = lz; i <= n; i++) //maxU = max(maxU, cumsumU[i]);
						if (cumsumU[i] > maxU) maxU = cumsumU[i];
					result[m + k] = min(minU - sumV, -maxU);
				}
			}
			else result[m + k] = numeric_limits<val_t>::max();
		}
	}

	// Signed halfspace depth of z_1,...z_m and x_i1,...x_il w.r.t. x_1,...x_n in R^2.
	// This routine assumes that the points are in general position
	template<typename val_t>
	void smD_GP(const double x[], const val_t valX[], const int& nn, const double z[], const int valZ[], const int& m,
		const int ind[], const int& l, val_t result[]) {
		int n{ 0 }, rz, iz;
		int curVal = +1;

		pointInfo<val_t>* p = new pointInfo<val_t>[nn + 1];
		for (int i = 0; i < nn; i++) {
			if (!isnan(x[2 * i])) {
				p[n].x = x + 2 * i;
				p[n].val = valX[i];
				n++;
			}
		}
		sort(p, p + n, [](pointInfo<val_t>& a, pointInfo<val_t>& b) { return (a.x[0] < b.x[0]); });


		int N = n + 1;
		angleInfo* angles = new angleInfo[N * (N - 1) / 2];
		int nAngle = 0;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				angles[nAngle].val = atan2(p[j].x[0] - p[i].x[0], p[i].x[1] - p[j].x[1]);
				angles[nAngle].p1 = i;
				angles[nAngle].p2 = j;
				nAngle++;
			}
		}
		sort(angles, angles + n * (n - 1) / 2, [](angleInfo& a, angleInfo& b) { return (a.val < b.val); });

		pointInfo<val_t>** op = new pointInfo<val_t> *[N];
		angleInfo** oa = new angleInfo * [N * (N - 1) / 2];

		val_t sumV;
		val_t* cumsumU = new val_t[N];
		int* ranks = new int[N];

		for (int k = 0; k < m; k++) {
			if (!isnan(z[2 * k])) {
				for (int i = 0; i < n; i++) op[i] = &p[i];
				for (int i = 0; i < n * (n - 1) / 2; i++) oa[i] = &angles[i];
				for (int i = 0; i < n; i++) ranks[i] = i;

				iz = n;
				for (int i = 0; i < n; i++) {
					if (eqPoints(z + 2 * k, p[i].x)) iz = i;
				}
				p[iz].x = z + 2 * k;

				if (iz == n) {
					N = n + 1;
					p[n].val = 0;
					op[n] = &p[n];
					inplace_merge(op, op + n, op + N, [](pointInfo<val_t>* a, pointInfo<val_t>* b) { return (a->x[0] < b->x[0]); });
					rz = 0;
					while (op[rz]->x != p[n].x) rz++;
					for (int i = rz; i < n; i++) ranks[i]++;
					ranks[n] = rz;

					angleInfo* newAngles = angles + n * (n - 1) / 2;
					for (int i = 0; i < n; i++) {
						if (i < rz) {
							newAngles[i].val = atan2(p[n].x[0] - p[i].x[0], p[i].x[1] - p[n].x[1]);
							newAngles[i].p1 = i;
							newAngles[i].p2 = n;
						}
						else {
							newAngles[i].val = atan2(p[i].x[0] - p[n].x[0], p[n].x[1] - p[i].x[1]);
							newAngles[i].p1 = n;
							newAngles[i].p2 = i;
						}
					}
					sort(newAngles, newAngles + n, [](angleInfo& a, angleInfo& b) { return (a.val < b.val); });
					for (int i = n * (n - 1) / 2; i < N * (N - 1) / 2; i++) oa[i] = &angles[i];
					inplace_merge(oa, oa + n * (n - 1) / 2, oa + N * (N - 1) / 2, [](angleInfo* a, angleInfo* b) { return (a->val < b->val); });
				}
				else N = n;

				if (valZ[k] != curVal) {
					curVal = -curVal;
					for (int i = 0; i < n; i++) p[i].val = -p[i].val;
				}

				cumsumU[0] = op[0]->val;
				for (int i = 1; i < N; i++) cumsumU[i] = cumsumU[i - 1] + op[i]->val;
				sumV = cumsumU[N - 1];

				int p1, p2, r1;
				val_t minU = sumV, maxU = 0;
				for (int i = 0; i < ranks[iz]; i++)
					if (cumsumU[i] > maxU) maxU = cumsumU[i];
				for (int i = ranks[iz]; i < N; i++)
					if (cumsumU[i] < minU) minU = cumsumU[i];

				for (int k = 0; k < N * (N - 1) / 2; k++) {
					p1 = oa[k]->p1;
					p2 = oa[k]->p2;
					r1 = ranks[p1];
					swap(ranks[p1], ranks[p2]);
					cumsumU[r1] += (p[p2].val - p[p1].val);

					if (ranks[iz] <= r1) {
						if (cumsumU[r1] < minU) minU = cumsumU[r1];
					}
					else {
						if (cumsumU[r1] > maxU) maxU = cumsumU[r1];
					}
				}
				result[k] = min(minU, sumV - maxU);
			}
			else result[k] = numeric_limits<val_t>::max();
		}
		for (int k = 0; k < l; k++) {
			if (!isnan(x[2 * ind[k]])) {
				iz = 0;
				while (!eqPoints(p[iz].x, x + 2 * ind[k])) iz++;

				if (sgn(valX[ind[k]]) != curVal) {
					curVal = -curVal;
					for (int i = 0; i < n; i++) p[i].val = -p[i].val;
				}

				for (int i = 0; i < n; i++) ranks[i] = i;

				cumsumU[0] = p[0].val;
				for (int i = 1; i < n; i++) cumsumU[i] = cumsumU[i - 1] + p[i].val;
				sumV = cumsumU[n - 1];

				int p1, p2, r1;
				val_t minU = sumV, maxU = 0;
				for (int i = 0; i < ranks[iz]; i++)
					if (cumsumU[i] > maxU) maxU = cumsumU[i];
				for (int i = ranks[iz]; i < n; i++)
					if (cumsumU[i] < minU) minU = cumsumU[i];

				for (int k = 0; k < nAngle; k++) {
					p1 = angles[k].p1;
					p2 = angles[k].p2;
					r1 = ranks[p1];
					swap(ranks[p1], ranks[p2]);
					cumsumU[r1] += (p[p2].val - p[p1].val);

					if (ranks[iz] <= r1) {
						if (cumsumU[r1] < minU) minU = cumsumU[r1];
					}
					else {
						if (cumsumU[r1] > maxU) maxU = cumsumU[r1];
					}
				}
				result[m + k] = min(minU, sumV - maxU);
			}
			else result[m + k] = numeric_limits<val_t>::max();
		}

		delete[] ranks;
		delete[] cumsumU;
		delete[] p;
		delete[] angles;
		delete[] op;
		delete[] oa;
	}

	// Signed halfspace depth for some of the data points x_i, the indices of these points are passed in 'ind', points are in R^3
	// Points need not be in general position
	// Smart pointers and a class for the set of positions is used

	template<typename val_t>
	void smD_nGP(const double x[], const val_t valX[], const int& nn,
		const double z[], const int valZ[], const int& m,
		const int ind[], const int& l, val_t result[]) {
		int n{ 0 }, rz, iz;
		int curVal = +1;

		unique_ptr<pointInfo<val_t>[]> p{ new pointInfo<val_t>[nn + 1] };
		for (int i = 0; i < nn; i++) {
			if (!isnan(x[2 * i])) {
				p[n].x = x + 2 * i;
				p[n].val = valX[i];
				n++;
			}
		}
		sort(p.get(), p.get() + n, cmpPointsEpsAtan<val_t>);
		adjustForMultiplePointsAtan(p.get(), n);
		int N = n + 1;

		unique_ptr<angleInfo[]> angles{ new angleInfo[N * (N - 1) / 2] };
		int nAngle = 0;
		for (int i = 0; i < n; i++) {
			for (int j = i + 1; j < n; j++) {
				if (fabs(atan(p[i].x[0]) - atan(p[j].x[0])) < eps_point_atan) angles[nAngle].val = 0;
				else
					angles[nAngle].val = atan2(p[j].x[0] - p[i].x[0], p[i].x[1] - p[j].x[1]);
				angles[nAngle].p1 = i;
				angles[nAngle].p2 = j;
				nAngle++;
			}
		}
		sort(angles.get(), angles.get() + n * (n - 1) / 2, cmpAngles);

		unique_ptr<pointInfo<val_t>* []> op{ new pointInfo<val_t> *[N] };
		unique_ptr<angleInfo* []> oa{ new angleInfo * [N * (N - 1) / 2] };

		val_t sumV;
		unique_ptr<val_t[]> cumsumU{ new val_t[N] };
		unique_ptr<int[]> ranks{ new int[N] };
		intSet position(N);

		for (int k = 0; k < m; k++) {
			if (!isnan(z[2 * k])) {
				for (int i = 0; i < n; i++) op[i] = &p[i];
				for (int i = 0; i < n * (n - 1) / 2; i++) oa[i] = &angles[i];
				for (int i = 0; i < n; i++) ranks[i] = i;

				iz = n;
				for (int i = 0; i < n; i++) {
					//if (hypot(atan(z[2 * k]) - atan(p[i].x[0]), atan(z[2 * k + 1]) - atan(p[i].x[1])) < eps) iz = i;
					if (eqPointsOnSphere(z + 2 * k, p[i].x)) iz = i;
				}
				p[iz].x = z + 2 * k;

				if (iz == n) {
					N = n + 1;
					p[n].val = 0;
					op[n] = &p[n];
					inplace_merge(op.get(), op.get() + n, op.get() + N, cmpPPointsEpsAtan<val_t>);
					rz = 0;
					while (op[rz]->x != p[n].x) rz++;
					for (int i = rz; i < n; i++) ranks[i]++;
					ranks[n] = rz;

					angleInfo* newAngles = angles.get() + n * (n - 1) / 2;
					for (int i = 0; i < n; i++) {
						if (i < rz) {
							if (fabs(atan(p[i].x[0]) - atan(p[n].x[0])) < eps_point_atan) newAngles[i].val = 0;
							else
								newAngles[i].val = atan2(p[n].x[0] - p[i].x[0], p[i].x[1] - p[n].x[1]);
							newAngles[i].p1 = i;
							newAngles[i].p2 = n;
						}
						else {
							if (fabs(atan(p[i].x[0]) - atan(p[n].x[0])) < eps_point_atan) newAngles[i].val = 0;
							else
								newAngles[i].val = atan2(p[i].x[0] - p[n].x[0], p[n].x[1] - p[i].x[1]);
							newAngles[i].p1 = n;
							newAngles[i].p2 = i;
						}
					}
					sort(newAngles, newAngles + n, bind(compAngles, _1, _2, ranks.get()));
					for (int i = n * (n - 1) / 2; i < N * (N - 1) / 2; i++) oa[i] = &angles[i];
					// Here the assertion that both sequences are sorted is violated in rare cases
					inplace_merge(oa.get(), oa.get() + n * (n - 1) / 2, oa.get() + N * (N - 1) / 2, bind(compPAngles, _1, _2, ranks.get()));
				}
				else N = n;

				if (valZ[k] != curVal) {
					curVal = -curVal;
					for (int i = 0; i < n; i++) p[i].val = -p[i].val;
				}

				cumsumU[0] = op[0]->val;
				for (int i = 1; i < N; i++) cumsumU[i] = cumsumU[i - 1] + op[i]->val;
				sumV = cumsumU[N - 1];

				int p1, p2, r1;
				val_t v1, v2, minU = sumV, maxU = 0;
				for (int i = 0; i < ranks[iz]; i++) //maxU = max(maxU, cumsumU[i]);
					if (cumsumU[i] > maxU) maxU = cumsumU[i];
				for (int i = ranks[iz]; i < N; i++) //minU = min(minU, cumsumU[i]);
					if (cumsumU[i] < minU) minU = cumsumU[i];

				for (int k = 0; k < N * (N - 1) / 2; k++) {
					p1 = oa[k]->p1;
					p2 = oa[k]->p2;
					r1 = ranks[p1];
					v1 = p[p1].val;
					v2 = p[p2].val;
					swap(ranks[p1], ranks[p2]);
					cumsumU[r1] += (v2 - v1);

					position.insert(r1);
					if ((k == N * (N - 1) / 2 - 1) || (oa[k]->val < oa[k + 1]->val - eps_angle)) {
						for (int p = 0; p < position.size(); p++) {
							if (ranks[iz] <= position[p])
								if (cumsumU[position[p]] < minU) minU = cumsumU[position[p]];
							if (ranks[iz] > position[p])
								if (cumsumU[position[p]] > maxU) maxU = cumsumU[position[p]];
						}
						position.clear();
					}
				}
				result[k] = min(minU, sumV - maxU);
			}
			else result[k] = numeric_limits<val_t>::max();
		}

		for (int k = 0; k < l; k++) {
			if (!isnan(x[2 * ind[k]])) {
				iz = 0;
				//while ((hypot(atan(p[iz].x[0]) - atan(x[2 * ind[k] + 0]), atan(p[iz].x[1]) - atan(x[2 * ind[k] + 1])) >= eps)) iz++;
				while (!eqPointsOnSphere(p[iz].x, x + 2 * ind[k])) iz++;
				//if (p[iz].val == -1)
				//if (valX[ind[k]] == -1)
				//    for (int i = 0; i < n; i++) p[i].val = -p[i].val;
				if (sgn(valX[ind[k]]) != curVal) {
					curVal = -curVal;
					for (int i = 0; i < n; i++) p[i].val = -p[i].val;
				}

				for (int i = 0; i < n; i++) ranks[i] = i;

				cumsumU[0] = p[0].val;
				for (int i = 1; i < n; i++) cumsumU[i] = cumsumU[i - 1] + p[i].val;
				sumV = cumsumU[n - 1];

				int p1, p2, r1;
				val_t v1, v2, minU = sumV, maxU = 0;
				for (int i = 0; i < ranks[iz]; i++) //maxU = max(maxU, cumsumU[i]);
					if (cumsumU[i] > maxU) maxU = cumsumU[i];
				for (int i = ranks[iz]; i < n; i++) //minU = min(minU, cumsumU[i]);
					if (cumsumU[i] < minU) minU = cumsumU[i];

				for (int k = 0; k < nAngle; k++) {
					p1 = angles[k].p1;
					p2 = angles[k].p2;
					r1 = ranks[p1];
					v1 = p[p1].val;
					v2 = p[p2].val;
					swap(ranks[p1], ranks[p2]);
					cumsumU[r1] += (v2 - v1);

					position.insert(r1);
					if ((k == nAngle - 1) || (angles[k].val < angles[k + 1].val - eps_angle)) {
						for (int p = 0; p < position.size(); p++) {
							if (ranks[iz] <= position[p])
								if (cumsumU[position[p]] < minU) minU = cumsumU[position[p]];
							if (ranks[iz] > position[p])
								if (cumsumU[position[p]] > maxU) maxU = cumsumU[position[p]];
						}
						position.clear();
					}
				}
				result[m + k] = min(minU, sumV - maxU);
			}
			else result[m + k] = numeric_limits<val_t>::max();
		}
	}

	/***************************************************************************************************************************/

	// Angular halfspace depth in d = 1, i.e., on { +1, -1 }

	template<typename val_t>
	void aHD_d1(const double x[], const val_t val[], int& n, const double z[], int& m, const int ind[], int& l, val_t result[]) {
		val_t pos{ 0 }, neg{ 0 };

		for (int i = 0; i < n; i++) {
			if (x[i] > eps) pos += val[i];
			if (x[i] < -eps) neg += val[i];
		}
		for (int i = 0; i < m + l; i++) result[i] = numeric_limits<val_t>::max();

		for (int i = 0; i < m; i++) {
			if (z[i] > eps) result[i] = min(pos, result[i]);
			if (z[i] < -eps) result[i] = min(neg, result[i]);
		}
		for (int i = 0; i < l; i++) {
			if (x[ind[i]] > eps) result[m + i] = min(pos, result[m + i]);
			if (x[ind[i]] < -eps) result[m + i] = min(neg, result[m + i]);
		}
	}

	// Angular halfspace depth in d = 2, i.e., on the circle

	template<typename val_t>
	void aHD_d2(const double x[], const val_t val[], int& n, const double z[], int& m, const int ind[], int& l, val_t result[], const AlgSubtype& alg) {
		unique_ptr<double[]> xx{ new double[n] };
		unique_ptr<double[]> zz{ new double[m] };
		unique_ptr<val_t[]> valX{ new val_t[n] };
		unique_ptr<int[]> valZ{ new int[m] };
		val_t totalMass{ }, negMass{ };
		AlgSubtype localAlg{ alg };
		if ((alg != single) && (alg != multiple)) localAlg = (m + l > limitSingle) ? multiple : single;
		// Parallelization of these loops was tested but is not advisable
		//#pragma omp parallel for default(none), shared(x,val,xx,valX,totalMass,negMass)
		for (int i = 0; i < n; i++) {
			if (!isnan(x[2 * i])) {
				if (fabs(x[2 * i + 1]) >= eps) {
					xx[i] = atan(x[2 * i] / x[2 * i + 1]);
					valX[i] = (x[2 * i + 1] > 0) ? val[i] : -val[i];
				}
				else {
					xx[i] = pi / 2; // numeric_limits<double>::infinity();
					valX[i] = (x[2 * i] > 0) ? val[i] : -val[i];
				}
				//#pragma omp critical
				{
					totalMass += val[i];
					if (valX[i] < 0) negMass += val[i];
				}
			}
			else xx[i] = nan("");
		}
		//#pragma omp parallel for default(none), shared(z,zz,valZ)
		for (int i = 0; i < m; i++) {
			if (!isnan(z[2 * i])) {
				if (fabs(z[2 * i + 1]) >= eps) {
					zz[i] = atan(z[2 * i] / z[2 * i + 1]);
					valZ[i] = (z[2 * i + 1] > 0) ? 1 : -1;
				}
				else {
					zz[i] = pi / 2; // numeric_limits<double>::infinity();
					valZ[i] = (z[2 * i] > 0) ? 1 : -1;
				}
			}
			else zz[i] = nan("");
		}
		switch (localAlg) {
		case single:   smD_d1s(xx.get(), valX.get(), n, zz.get(), valZ.get(), m, ind, l, result); break;
		case multiple: smD_d1m(xx.get(), valX.get(), n, zz.get(), valZ.get(), m, ind, l, result); break;
		default:       smD_d1m(xx.get(), valX.get(), n, zz.get(), valZ.get(), m, ind, l, result);
		}
		for (int i = 0; i < m; i++)
			if (!isnan(zz[i]))
				result[i] += (valZ[i] == 1) ? negMass : totalMass - negMass;
			else
				result[i] = numeric_limits<val_t>::max();
		for (int i = 0; i < l; i++)
			if (!isnan(xx[ind[i]]))
				result[m + i] += (valX[ind[i]] > 0) ? negMass : totalMass - negMass;  // here should be ind[i] instzead of simply i
			else
				result[m + i] = numeric_limits<val_t>::max();
	}

	// Angular halfspace depth in d = 3, i.e., on the 2-sphere

	template<typename val_t>
	void aHD_d3(const double xx[], const val_t val[], int& n, const double zz[], int& m, const int ind[], int& l, val_t result[], const AlgSubtype& alg) {
		unique_ptr<double[]> x{ new double[3 * n] };
		unique_ptr<double[]> z{ new double[3 * m] };
		unique_ptr<double[]> px{ new double[2 * n] };
		unique_ptr<double[]> pz{ new double[2 * m] };
		unique_ptr<val_t[]> valX{ new val_t[n] };
		unique_ptr<int[]> valZ{ new int[m] };
		val_t totalMass{ }, negMass{ };

		AlgSubtype localAlg{ (alg == standard) ? nGP : alg };

		static mt19937 gen(1234);
		normal_distribution<double> normal(0, 1);

		uninitialized_copy(xx, xx + 3 * n, x.get());
		uninitialized_copy(zz, zz + 3 * m, z.get());
		while (pointsOnEquator(x.get(), n) || pointsOnEquator(z.get(), m)) {
			double u[3];
			double norm{ 0.0 };
			for (int i = 0; i < 3; i++) {
				u[i] = normal(gen);
				norm += u[i] * u[i];
			}
			norm = sqrt(norm);
			for (int i = 0; i < 3; i++) u[i] /= norm;
			for (int i = 0; i < n; i++) Householder(&x[3 * i], u);
			for (int i = 0; i < m; i++) Householder(&z[3 * i], u);
		}
		for (int i = 0; i < n; i++) {
			if (!isnan(x[3 * i])) {
				px[2 * i + 0] = x[3 * i + 0] / x[3 * i + 2];
				px[2 * i + 1] = x[3 * i + 1] / x[3 * i + 2];
				totalMass += val[i];
				valX[i] = (x[3 * i + 2] > 0) ? val[i] : -val[i];
				if (valX[i] < 0) negMass += val[i];
			}
			else px[2 * i] = nan("");
		}
		for (int i = 0; i < m; i++) {
			if (!isnan(z[3 * i])) {
				pz[2 * i + 0] = z[3 * i + 0] / z[3 * i + 2];
				pz[2 * i + 1] = z[3 * i + 1] / z[3 * i + 2];
				valZ[i] = (z[3 * i + 2] > 0) ? 1 : -1;
			}
			else pz[2 * i] = nan("");
		}
		switch (localAlg) {
		case nGP: smD_nGP(px.get(), valX.get(), n, pz.get(), valZ.get(), m, ind, l, result); break;
		case GP:  smD_GP (px.get(), valX.get(), n, pz.get(), valZ.get(), m, ind, l, result); break;
		default:  smD_nGP (px.get(), valX.get(), n, pz.get(), valZ.get(), m, ind, l, result);
		}

		for (int i = 0; i < m; i++)
			if (!isnan(pz[2 * i]))
				result[i] += (z[3 * i + 2] > 0) ? negMass : totalMass - negMass;
			else
				result[i] = numeric_limits<val_t>::max();
		for (int i = 0; i < l; i++)
			if (!isnan(px[2 * ind[i]]))
				result[m + i] += (x[3 * ind[i] + 2] > 0) ? negMass : totalMass - negMass;
			else
				result[m + i] = numeric_limits<val_t>::max();
	}

	/***************************************************************************************************************************/

	// Angular halfspace depth in arbitrary dimension, recursive algorithm

	template<typename val_t>
	void aHD_Rec(const double x[], const val_t val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, val_t result[], const AlgType& alg, const int& target, const AlgSubtype& algSub) {
		if (d > target) {
			int pd{ d - 1 };

			for (int i = 0; i < m + l; i++) result[i] = numeric_limits<val_t>::max();
			#pragma omp parallel for default(none), shared(n,d,pd,x,z,m,val,target,ind,l,alg,algSub,result)
			for (int i = 0; i < n; i++) {
				if (!isnan(x[i * d])) {
					unique_ptr<double[]> pz{ new double[m * (d - 1)] {} };
					unique_ptr<double[]> px{ new double[n * (d - 1)] {} };
					unique_ptr<val_t[]> res{ new val_t[m + l] };
					// project data points
					bool pos;
					val_t nPos{ 0 }, nNeg{ 0 };
					for (int j = 0; j < m; j++)
						HouseholderProjection(z + j * d, x + i * d, d, pz.get() + j * (d - 1), pos);
					for (int j = 0; j < n; j++) {
						HouseholderProjection(x + j * d, x + i * d, d, px.get() + j * (d - 1), pos);
						if (pos) nPos += val[j]; else nNeg += val[j];
					}
					val_t minP = min(nPos, nNeg);
					aHD_Rec(px.get(), val, n, pd, pz.get(), m, ind, l, res.get(), alg, target, algSub);
					#pragma omp critical
					{
						for (int i = 0; i < m + l; i++)
							if (res[i] < numeric_limits<val_t>::max())
								result[i] = min(result[i], res[i] + minP);
					}
				}
			}
		}
		else {
			switch (d) {
			case 1: aHD_d1(x, val, n, z, m, ind, l, result); return;
			case 2: aHD_d2(x, val, n, z, m, ind, l, result, algSub); return;
			case 3: aHD_d3(x, val, n, z, m, ind, l, result, algSub); return;
			}
		}
	}

	// Angular halfspace depth in arbitrary dimension, combinatorial algorithm

	template<typename val_t>
	void aHD_Comb(const double x[], const val_t val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, val_t result[], const AlgType& alg, const int& target, const AlgSubtype& algSub) {
		const int pd{ d - target };

		if (d > target) {
			for (int i = 0; i < m + l; i++) result[i] = numeric_limits<val_t>::max();
			#pragma omp parallel for default(none), shared(n,d,pd,x,z,m,val,target,ind,l,alg,algSub,result)
			for (int r = 0; r < nCr(n, pd); r++) {
				unique_ptr<double[]> pz{ new double[m * target] };
				unique_ptr<double[]> px{ new double[n * target] };
				unique_ptr<val_t[]> res{ new val_t[m + l] };
				unique_ptr<int[]> comb{ new int[pd] };
				unRank(r, pd, comb.get());
				val_t HD;
				cMatrix A(d, pd);
				// copy vectors to matrix A
				for (int i = 0; i < d; i++)
					for (int j = 0; j < pd; j++) A(i, j) = x[comb[j] * d + i];

				cMatrix Q;
				unique_ptr<double[]> D{ new double[pd] };
				unique_ptr<int[]> p{ new int[pd] };
				int rank{};
				A.QR_Factorization(A, D.get(), p.get(), false, &rank);
				if (rank == pd) {
					// project data points
					A.AccumulateHH(D.get(), 0, pd, Q, false, pd, d);
					project(Q, z, m, pz.get());
					A.AccumulateHH(D.get(), 0, pd, Q, false, 0, d);
					project2(Q, x, val, n, target, px.get(), HD);

					switch (target) {
					case 1: aHD_d1(px.get(), val, n, pz.get(), m, ind, l, res.get()); break;
					case 2: aHD_d2(px.get(), val, n, pz.get(), m, ind, l, res.get(), algSub); break;
					case 3: aHD_d3(px.get(), val, n, pz.get(), m, ind, l, res.get(), algSub); break;
					}

					#pragma omp critical
					{
						for (int i = 0; i < m + l; i++)
							if (res[i] < numeric_limits<val_t>::max()) {
								result[i] = min(result[i], res[i] + HD);
							}
					}
				}
			}
		}
		else {
			switch (d) {
			case 1: aHD_d1(x, val, n, z, m, ind, l, result); return;
			case 2: aHD_d2(x, val, n, z, m, ind, l, result, algSub); return;
			case 3: aHD_d3(x, val, n, z, m, ind, l, result, algSub); return;
			}
		}
	}

	/***************************************************************************************************************************/

	// Wrapper for aHD_Comb and aHD_Rec that does the following preprocessing.
	// Checks whether there are multiple data points.
	// Checks whether the points are contained in some lower dimensional subspace. If that is the case, the following
	// computations are done in that lower dimensional subspace.

	template<typename val_t>
	void aHD(const double xx[], const val_t val[], int& nn, const int& d, const double z[], int& m, const int ind[], int& l, val_t result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads) {

		int nThreadsSav{ omp_get_max_threads() };
		int NestedSav{ omp_get_nested() };
		omp_set_num_threads(NumberOfCoresToUse = (nThreads == 0) ? omp_get_num_procs() : min(nThreads, omp_get_num_procs()));
		omp_set_nested(0);

		unique_ptr<int[]> newInd{ new int[l] };
		unique_ptr<int[]> perm{ new int[nn] };

		int n{ nn };
		unique_ptr<pointInfoP<val_t>[]> xp{ new pointInfoP<val_t>[nn] };
		for (int i = 0; i < nn; i++) {
			xp[i].x = &xx[i * d];
			xp[i].val = val[i];
			xp[i].index = i;
			perm[i] = i;
		}
		if (d > 2) {
			sort(xp.get(), xp.get() + nn, bind(ltPointsRounded<val_t>, _1, _2, d));
			int i{ 0 }, j, k{ 0 };
			while (i < nn) {
				j = i + 1;
				while ((j < nn) && (eqPointsRounded(xp[i], xp[j], d))) j++;
				xp[k].x = xp[i].x;
				xp[k].val = xp[i].val;
				perm[xp[i].index] = k;
				for (int l = i + 1; l < j; l++) {
					xp[k].val += xp[l].val;
					perm[xp[l].index] = k;
				}
				i = j;
				k++;
			}
			n = k;
		}
		unique_ptr<double[]> x{ new double[n * d] };
		unique_ptr<val_t[]> valX{ new val_t[n] };
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < d; j++) x[i * d + j] = xp[i].x[j];
			valX[i] = xp[i].val;
		}
		for (int i = 0; i < l; i++) newInd[i] = perm[ind[i]];

		cMatrix A(d, n);
		// copy vectors to matrix A
		for (int i = 0; i < d; i++)
			for (int j = 0; j < n; j++) A(i, j) = x[j * d + i];
		cMatrix Q;
		unique_ptr<double[]> D{ new double[min(n,d)] };
		unique_ptr<int[]> p{ new int[n] };
		int rank{};
		A.QR_Factorization(A, D.get(), p.get(), false, &rank);
		unique_ptr<int[]> pInv{ new int[n] };
		unique_ptr<int[]> pind{ new int[l] };
		unique_ptr<double[]> zz{ new double[m * d] };
		unique_ptr<double[]> px{ new double[n * rank] };
		unique_ptr<double[]> pz{ new double[m * rank] };
		unique_ptr<bool[]> aHD{ new bool[m] {} };
		bool compHD{ false };
		if (rank < d) {
			// Rank deficiency!
			A.AccumulateHH(D.get(), 0, rank, Q, false, 0, rank);
			project(Q, x.get(), n, px.get());
			// project data points z_1,...,z_m
			A.AccumulateHH(D.get(), 0, rank, Q, false, 0, d);
			project(Q, z, m, zz.get());
			// for each point z_i check whether it is in sp(x_1,\dots,x_n)
			for (int i = 0; i < m; i++) {
				double norm{ 0.0 };
				for (int j = rank; j < d; j++) norm += zz[i * d + j] * zz[i * d + j];
				if (sqrt(norm) < eps) { // z is in sp(x_1,...,x_n)
					for (int j = 0; j < rank; j++) pz[i * rank + j] = zz[i * d + j];
					aHD[i] = true;
				}
				else {
					for (int j = 0; j < rank; j++) pz[i * rank + j] = nan("");
					aHD[i] = false;
					compHD = true;
				}
			}
			switch (alg) {
			case comb: aHD_Comb(px.get(), valX.get(), n, rank, pz.get(), m, newInd.get(), l, result, alg, target, algSub); break;
			case rec:  aHD_Rec (px.get(), valX.get(), n, rank, pz.get(), m, newInd.get(), l, result, alg, target, algSub); break;
			}
			if (compHD) { // compute the halfspace depth of 0 w.r.t. the points px_1,...px_n
				val_t res = DataDepth::nHD_Comb2(px.get(), valX.get(), n, rank);
				for (int i = 0; i < m; i++)
					if (!aHD[i]) result[i] = res;
			}
		}
		else switch (alg) {
		case comb: aHD_Comb(x.get(), valX.get(), n, d, z, m, newInd.get(), l, result, alg, target, algSub); break;
		case rec:  aHD_Rec (x.get(), valX.get(), n, d, z, m, newInd.get(), l, result, alg, target, algSub); break;
		}
		omp_set_num_threads(nThreadsSav);
		omp_set_nested(NestedSav);
	}

	template<typename val_t>
	void aHD(const double xx[], nullptr_t, int& nn, const int& d, const double z[], int& m, const int ind[], int& l, val_t result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads) {
		unique_ptr<val_t[]> val{ new val_t[nn] };
		if (is_integral<val_t>::value)
			for (int i = 0; i < nn; i++) val[i] = 1;
		else
			for (int i = 0; i < nn; i++) val[i] = 1.0 / nn;
		aHD(xx, val.get(), nn, d, z, m, ind, l, result, alg, target, algSub, nThreads);
	}

	template<typename val_t>
	void aHD(const double xx[], int& nn, const int& d, const double z[], int& m, const int ind[], int& l, val_t result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads) {
		aHD(xx, nullptr, nn, d, z, m, ind, l, result, alg, target, algSub, nThreads);
	}

	// Explicit instantiation
	template void aHD<int>   (const double x[], const int    val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, int    result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads);
	template void aHD<double>(const double x[], const double val[], int& n, const int& d, const double z[], int& m, const int ind[], int& l, double result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads);

	template void aHD<int>   (const double x[],                     int& n, const int& d, const double z[], int& m, const int ind[], int& l, int    result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads);
	template void aHD<double>(const double x[],                     int& n, const int& d, const double z[], int& m, const int ind[], int& l, double result[], const AlgType& alg, const int& target, const AlgSubtype& algSub, const int& nThreads);


}

