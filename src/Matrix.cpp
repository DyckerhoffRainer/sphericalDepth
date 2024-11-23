#include <memory>
#include <numeric>
#include "Matrix.h"

namespace dyMatrixClass {


	using namespace std;

	//********** Matrix format flag handling *****************************************************

	enum class matrix_format_flags { fShort, fLong };

	long& matrix_format_flag(ios_base& s) {
		static int my_index = ios_base::xalloc();
		return s.iword(my_index);
	}

	long get_matrix_format_flag(ios_base& s) {
		return matrix_format_flag(s);
	}

	void set_matrix_format_flag(ios_base& s, long n) {
		matrix_format_flag(s) = n;
	}

	static void set_matrix_format(ios_base& s, matrix_format_flags mat_fmt)	{
		matrix_format_flag(s) = (long)mat_fmt;
	}

	ostream& shortMat(ostream & os) {
		set_matrix_format(os, matrix_format_flags::fShort);
		return os;
	}

	ostream& longMat(ostream & os) {
		set_matrix_format(os, matrix_format_flags::fLong);
		return os;
	}


	//********** Matrix class ********************************************************************

#undef _MESSAGES_

	cMatrix::cMatrix(cMatrix&& A)											// Move constructor
		: m{ A.m }, n{ A.n }, elems{ move(A.elems) } {
		//cout << "Matrix move constructor" << endl;
	};

	cMatrix& cMatrix::operator=(cMatrix&& A) {								// Move assignment
		m = A.m;
		n = A.n;
		elems = move(A.elems);
		//cout << "Matrix move assignment" << endl;
		return *this;
	};

	cMatrix::cMatrix(const cMatrix& A)										// Copy constructor
		: m{ A.m }, n{ A.n }, elems{ new double[A.m*A.n] } {
		std::uninitialized_copy_n(A.elems.get(), m*n, elems.get());
		//cout << "Matrix copy constructor" << endl;
	}

	cMatrix& cMatrix::operator=(const cMatrix& A) {							// Copy assignment
		if (this != &A) {
			m = A.m;
			n = A.n;
			elems = unique_ptr<double[]>(new double[m*n]);
			std::uninitialized_copy_n(A.elems.get(), m*n, elems.get());
		}
		//cout << "Matrix copy assignment" << endl;
		return *this;
	}

	cMatrix::cMatrix(std::initializer_list<std::initializer_list<double>> args) {					// List initializer constructor
		m = args.size();
		n = 0;
		for (int i = 0; i < args.size(); i++) n = std::max(n, (int)args.begin()[i].size());
		elems = unique_ptr<double[]>(new double[m*n]{});
		for (int i = 0; i < args.size(); i++) std::uninitialized_copy(args.begin()[i].begin(), args.begin()[i].end(), elems.get() + i*n);
		//cout << "Matrix 2-dim list constructor" << endl;
	}

	cMatrix::cMatrix(std::initializer_list<double> args) {										// List initializer constructor (column vector)
		m = args.size();
		n = 1;
		elems = unique_ptr<double[]>(new double[m]{});
		std::uninitialized_copy(args.begin(), args.end(), elems.get());
		//cout << "Matrix 1-dim list constructor" << endl;
	}

	void cMatrix::SetSize(int _m, int _n) {
		m = _m;
		n = _n;
		elems = unique_ptr<double[]>(new double[m*n]{});
	};

	void cMatrix::Resize(int _m, int _n, bool keepData) {
		if (keepData) m = _m;
	}

	cMatrix cMatrix::Row(int i) const {
		cMatrix x(1, n);
		uninitialized_copy_n(elems.get() + i*n, n, x.elems.get());
		return x;
	}

	cMatrix cMatrix::Col(int j) const {
		cMatrix x(m, 1);
		for (int i = 0; i < m; i++) x(i, 0) = elems[i*n + j];
		return x;
	}

	cMatrix cMatrix::Diag() const {
		cMatrix x(min(m,n), 1);
		for (int i = 0; i < min(m,n); i++) x(i, 0) = elems[i*n + i];
		return x;
	}

	cMatrix cMatrix::SubMatrix(int i1, int i2, int j1, int j2) const {
		cMatrix A(i2 - i1 + 1, j2 - j1 + 1);
		for (int i = i1; i <= i2; i++)
			for (int j = j1; j <= j2; j++)
				A(i - i1, j - j1) = elems[i*n + j];
		return A;
	}


	cMatrix& cMatrix::apply(std::function<double(double)> f) {
		for (int i = 0; i < m*n; i++) elems[i] = f(elems[i]);
		return *this;
	}

	cMatrix& cMatrix::operator*=(double c) {
		for (int i = 0; i < m*n; i++) elems[i] *= c;
		return *this;
	}

	cMatrix& cMatrix::operator/=(double c) {
		for (int i = 0; i < m*n; i++) elems[i] /= c;
		return *this;
	}

	cMatrix& cMatrix::operator+=(const cMatrix& B) {
		assert((m == B.m) && (n == B.n));
		for (int i = 0; i < m*n; i++) elems[i] += B.elems[i];
		return *this;
	}

	cMatrix& cMatrix::operator-=(const cMatrix& B) {
		assert((m == B.m) && (n == B.n));
		for (int i = 0; i < m*n; i++) elems[i] -= B.elems[i];
		return *this;
	}

	cMatrix& cMatrix::operator*=(const cMatrix& B) {
		assert((n == B.m));
		cMatrix& A{ *this };
		cMatrix tmp(m, B.n);
		for (int i = 0; i < tmp.m; i++) {
			for (int j = 0; j < tmp.n; j++) {
				for (int k = 0; k < n; k++) tmp(i, j) += A(i, k) * B(k, j);
			}
		}
		return *this = std::move(tmp);
	}

	cMatrix cMatrix::operator-() {
		cMatrix B(m, n);
		for (int i = 0; i < m*n; i++) B.elems[i] = -elems[i];
		return B;
	}

	cMatrix operator*(const cMatrix& A, double c) {
		cMatrix res(A.m, A.n);
		for (int i = 0; i < res.m * res.n; i++) res.elems[i] = c * A.elems[i];
		return res;
	}

	cMatrix operator*(double c, const cMatrix& A) {
		cMatrix res(A.m, A.n);
		for (int i = 0; i < res.m * res.n; i++) res.elems[i] = c*A.elems[i];
		return res;
	}

	cMatrix operator/(const cMatrix& A, double c) {
		cMatrix res(A.m, A.n);
		for (int i = 0; i < res.m * res.n; i++) res.elems[i] = A.elems[i] / c;
		return res;
	}

	cMatrix operator+(const cMatrix& A, const cMatrix& B) {
		assert((A.m == B.m) && (A.n == B.n));
		cMatrix res(A.m, A.n);
		for (int i = 0; i < res.m * res.n; i++) res.elems[i] = A.elems[i] + B.elems[i];
		return res;
	}

	cMatrix operator-(const cMatrix& A, const cMatrix& B) {
		assert((A.m == B.m) && (A.n == B.n));
		cMatrix res(A.m, A.n);
		for (int i = 0; i < res.m * res.n; i++) res.elems[i] = A.elems[i] - B.elems[i];
		return res;
	}

	cMatrix operator*(const cMatrix& A, const cMatrix& B) {
		assert((A.n == B.m));
		cMatrix res(A.m, B.n);
		for (int i = 0; i < res.m; i++)
			for (int j = 0; j < res.n; j++)
				for (int k = 0; k < A.n; k++)
					res(i, j) += A(i, k) * B(k, j);
		return res;
	}

	bool operator==(const cMatrix& A, const cMatrix& B) {
		if ((A.m != B.m) || (A.n != B.n)) return false;
		for (int i = 0; i < A.m * A.n; i++)
			if (A.elems[i] != B.elems[i]) return false;
		return true;
	}

	bool operator!=(const cMatrix& A, const cMatrix& B) {
		if ((A.m != B.m) || (A.n != B.n)) return true;
		for (int i = 0; i < A.m * A.n; i++)
			if (A.elems[i] != B.elems[i]) return true;
		return false;
	}

	bool operator<=(const cMatrix& A, const cMatrix& B) {
		assert((A.m == B.m) && (A.n == B.n));
		for (int i = 0; i < A.m * A.n; i++)
			if (A.elems[i] > B.elems[i]) return false;
		return true;
	}

	bool operator>=(const cMatrix& A, const cMatrix& B) {
		assert((A.m == B.m) && (A.n == B.n));
		for (int i = 0; i < A.m * A.n; i++)
			if (A.elems[i] < B.elems[i]) return false;
		return true;
	}

	bool operator<(const cMatrix& A, const cMatrix& B) {
		assert((A.m == B.m) && (A.n == B.n));
		for (int i = 0; i < A.m * A.n; i++)
			if (A.elems[i] >= B.elems[i]) return false;
		return true;
	}

	bool operator>(const cMatrix& A, const cMatrix& B) {
		assert((A.m == B.m) && (A.n == B.n));
		for (int i = 0; i < A.m * A.n; i++)
			if (A.elems[i] <= B.elems[i]) return false;
		return true;
	}

	cMatrix::operator double() const {
		assert((m == 1) && (n == 1));
		return elems[0];
	}

	cMatrix Trans(const cMatrix& A) {
		cMatrix res(A.n, A.m);
		for (int i = 0; i < res.m; i++)
			for (int j = 0; j < res.n; j++) res(i, j) = A(j, i);
		return res;
	}

	std::ostream& operator<<(std::ostream& os, const cMatrix& A) {
		if (get_matrix_format_flag(os) == (long)matrix_format_flags::fShort) {
			os << "{";
			for (int i = 0; i < A.m; i++) {
				os << "{";
				for (int j = 0; j < A.n; j++) {
					os << ((abs(A(i,j)) < 1e-14) ? 0 : A(i, j));
					if (j < A.n - 1) os << ",";
					else os << "}";
				}
				if (i < A.m - 1) os << ",";
			}
			os << "}";
		}
		else {
			for (int i = 0; i < A.m; i++) {
				os << "[";
				for (int j = 0; j < A.n; j++) {
					os << std::setw(10) << ((abs(A(i, j)) < 1e-14) ? 0 : A(i, j));
					if (j < A.n - 1) os << ",";
				}
				os << "]" << std::endl;
			}
		}
		return os;
	}

	std::istream& operator>>(std::istream& os, const cMatrix& A) {
		for (int i = 0; i < A.m; i++) {
			for (int j = 0; j < A.n; j++) os >> A(i, j);
			}
		return os;
	}

	//***************** Vector *******************************************************

	cVector& cVector::operator=(cMatrix&& A) {								// Move assignment
		cMatrix::operator=(move(A));
		return *this;
	};

	cVector& cVector::operator=(const cMatrix& A) {							// Copy assignment
		if (this != &A) {
			cMatrix::operator=(A);
		}
		return *this;
	}

	cVector::cVector(std::initializer_list<double> args) {					// List initializer constructor (column vector)
		m = args.size();
		n = 1;
		elems = unique_ptr<double[]>(new double[m]{});
		std::uninitialized_copy(args.begin(), args.end(), elems.get());
		//cout << "Vector list constructor" << endl;
	}

	double InnerProduct(const cVector& A, const cVector& B) {
		return inner_product(A.elems.get(), A.elems.get() + A.m, B.elems.get(), 0.0);
	}


	//***************** Matrixnormen *******************************************************

	double FrobeniusNorm(const cMatrix& A) {
		return sqrt(inner_product(A.elems.get(), A.elems.get() + A.m * A.n, A.elems.get(), 0.0));
	}

	double RowSumNorm(const cMatrix& A) {
		double norm{};
		for (int i = 0; i < A.m; i++) {
			double sum{};
			for (int j = 0; j < A.n; j++) sum += abs(A(i, j));
			norm = max(sum, norm);
		}
		return norm;
	}

	double ColSumNorm(const cMatrix& A) {
		double norm{};
		for (int j = 0; j < A.n; j++) {
			double sum{};
			for (int i = 0; i < A.m; i++) sum += abs(A(i, j));
			norm = max(sum, norm);
		}
		return norm;
	}

	cMatrix& cMatrix::SetEpsToZero(double eps) {
		for (int i = 0; i < m*n; i++)
			if (abs(elems[i]) < eps) elems[i] = 0;
		return *this;
	}

	double Sum(const cMatrix& A) {
		double sum{ 0.0 };
		for (int i = 0; i < A.m*A.n; i++) sum += A.elems[i];
		return sum;
	}


/************************************************************************************************

	Cholesky decomposition of a positive semidefinite matrix, A = GG^t.
	Only the lower trinagular part of A{ *this } is used.
	G is a lower triangular matrix, G may coincide with A in which case A is overwritten by G.

	This routine also handles the case where A is only positive semidefinite (but not positive
	definite).

	Parameters:
		A{*this} matrix to be factored.
		G       on exit lower triangular matrix s.t. GG^T = A.
		rank	rank of A (or equivalent G)

************************************************************************************************/

	void cMatrix::Cholesky(cMatrix& G, int* rank) const {
		const double eps = 1e-12;
		if (rank != nullptr) *rank = 0;

		G = *this;

		for (int j = 0; j < n; j++) {
			double eps1{ eps * abs(G(j, j)) };
			for (int k = 0; k < j; k++) G(j, j) -= G(j, k) * G(j, k);
			if (abs(G(j, j)) <= eps1)
				G(j, j) = 0;
			else if (G(j, j) > 0) {
				G(j, j) = sqrt(G(j, j));
				if (rank != nullptr) (*rank)++;
			}
			else throw definiteness_error("Matrix not positive semidefinite!");
			for (int i = j + 1; i < n; i++) {
				double eps2{ eps * abs(G(i, i)) };
				for (int k = 0; k < j; k++) G(i, j) -= G(i, k) * G(j, k);
				if (G(j, j) > 0)
					G(i, j) /= G(j, j);
				else if (abs(G(i, j)) <= eps1 * eps2)
					G(i, j) = 0;
				else throw definiteness_error("Matrix not positive semidefinite!");
			}
		}
		// Elemente über der Diagonalen gleich null setzen
		for (int i = 0; i < n; i++)
			for (int j = i + 1; j < n; j++) G(i, j) = 0;
	};

/************************************************************************************************

	Computes the LU-factorization of the given sqare matrix: PA = LU
	L is a lower triangular matrix with ones on the diagonal.
	U is an upper triangular matrix
	L and U are stored in LU, the ones on the diagonal of L are not stored.
	P describes the row permutation resulting from pivoting.

	Parameters:
		A{*this} matrix to be factored. A is assumed to be a square matrix.
		LU       stores matrices L and U as described above, LU may coincide with A
		p        stores the row permutation
		det      determinant of A

************************************************************************************************/

	void cMatrix::LU_Factorization(cMatrix& LU, int p[], double& det) const {
		const double eps{ 1e-15 };
		unique_ptr<double[]> scale{ new double[n]{} };
		// Copy A to LU
		LU = *this;
		// Compute scaling factors for the rows
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) scale[i] += LU(i, j) * LU(i, j);
			scale[i] = 1 / sqrt(scale[i]);
		}
		// Initialize det
		det = 1;
		// Outer loop
		for (int k = 0; k < n; k++) {
			int imax{ k };
			double max{ 0 };
			// Calculate the elements L(k, k), ..., L(n, k) }
			for (int i = k; i < n; i++) {
				for (int j = 0; j < k; j++) LU(i, k) -= LU(i, j) * LU(j, k);
				double tmp{ abs(LU(i, k) * scale[i]) };
				// Determine maximum of the scaled absolute values of L(k, k), ..., L(n, k)
				if (tmp > max) {
					max = tmp;
					imax = i;
				}
			}
			// If imax != k exchange rows k and imax
			if (imax != k) {
				det = -det;
				LU.XChangeRows(k, imax);
				scale[imax] = scale[k];
			}
			// store number of pivot row
			p[k] = imax;
			// Update determinant
			det *= LU(k, k);
			// continue if matrix is not singular
			if (max < 8 * eps) throw runtime_error("Matrix singular!");
			// Calculate the elements U(k, k + 1), ..., U(k, n)
			for (int j = k + 1; j < n; j++)
				for (int i = 0; i < k; i++) LU(k, j) -= LU(k, i) * LU(i, j);
			for (int i = k + 1; i < n; i++) LU(i, k) /= LU(k, k);
		}
	}


/************************************************************************************************

	Householder transformation of the given matrix A{ *this }.
	The Hoiuseholder matrix H is determined s.t. in the case 'Transpose = false'

		H * | A[first, pos]	| = | *	| ,
			|		.		|	| 0 |
			|		.		|	| . |
			|		.		|	| . |
			| A[last, pos]	|	| 0 |

	and in the case 'Transpose = true'

		| A[pos, first] | ^T * H =	| *	|^T.
		|		.		|			| 0 |
		|		.		|			| . |
		|		.		|			| . |
		| A[pos, last]	|			| 0 |

	The Householder vector v is stored at the index positionbs
		A[first, pos], ..., A[last, pos]	(Transpose = false)
	or
		A[pos, first], ..., A[pos, last]	(Transpose = true).

	The matrix A is then tranformed by
		A = H * A	(Transpose = false)
	or
		A = A * H	(Transpose = true).
	Here it is assumed that
		A[i, j] = 0, i = first, ..., last, j = 0, ..., pos - 1  (Transpose = false)
	or
		A[i, j] = 0, i = 0, ..., pos - 1, j = first, ..., last	(Transpose = true).
	The new value * is returned in 'New'.

************************************************************************************************/

	void cMatrix::Householder1(double& New, int first, int last, int pos, bool transpose) {
		cMatrix& A{ *this };
		double scale{ 0 };
		double g{ 0 };
		if (!transpose) {
			for (int i = first; i <= last; i++) scale += abs(A(i, pos));
			if (scale != 0) {
				double r{ 0 };
				for (int i = first; i <= last; i++) {
					A(i, pos) /= scale;
					r += A(i, pos) * A(i, pos);
				}
				double f{ A(first, pos) };
				g = -copysign(sqrt(r), f);
				double h{ f * g - r };
				A(first, pos) = f - g;
				for (int j = pos + 1; j < n; j++) {
					double s{ 0 };
					for (int i = first; i <= last; i++) s += A(i, pos) * A(i, j);
					f = s / h;
					for (int i = first; i <= last; i++) A(i, j) += f * A(i, pos);
				}
				for (int i = first; i <= last; i++) A(i, pos) *= scale;
			}
		}
		else {
			for (int j = first; j <= last; j++) scale += abs(A(pos, j));
			if (scale != 0) {
				double r{ 0 };
				for (int j = first; j <= last; j++) {
					A(pos, j) /= scale;
					r += A(pos, j) * A(pos, j);
				}
				double f{ A(pos, first) };
				g = -copysign(sqrt(r), f);
				double h{ f * g - r };
				A(pos, first) = f - g;
				for (int i = pos + 1; i < m; i++) {
					double s{ 0 };
					for (int j = first; j <= last; j++) s += A(pos, j) * A(i, j);
					f = s / h;
					for (int j = first; j <= last; j++) A(i, j) += f * A(pos, j);
				}
				for (int j = first; j <= last; j++) A(pos, j) *= scale;
			}
		}
		New = g * scale;
	}

/************************************************************************************************

	Householder transformation of the given matrix A{ *this }.
	The Hoiuseholder matrix H is determined s.t.

		H * | A[first, pos]	| = | *	|.
			|		.		|	| 0 |
			|		.		|	| . |
			|		.		|	| . |
			| A[last, pos]	|	| 0 |


	The Householder vector v is stored at the index positionbs
		A[first, pos], ..., A[last, pos].
	The matrix A is then transformed by A = H * A * H.
	It is assumed that A[i, j] = 0, i = first, ..., last, j = 0, ..., pos - 1
	The new value * is returned in 'New'.

	The parameter 'Symm' is 'true' if A is symmetric and 'false' otherwise.

************************************************************************************************/

	void cMatrix::Householder2(double& New, int first, int last, int pos, bool Symm) {
		cMatrix& A{ *this };
		double g{ 0 };
		// Skalierung mit der 1 - Norm des Spaltenvektors
		double scale{ 0 };
		for (int i = first; i <= last; i++) scale += abs(A(i, pos));
		if (scale != 0) {
			// Householder - Transformation I - \beta v v^t berechnen.Der Vektor v
			// wird in A(first, pos], ..., A(last, pos] gespeichert.
			double r{ 0 };
			for (int i = first; i <= last; i++) {
				A(i, pos) /= scale;
				r += A(i, pos) * A(i, pos);
			}
			double f{ A(first, pos) };
			g = -copysign(sqrt(r), f);          // g = -sign(A(first, pos)) * || u ||
			double h{ f * g - r };              // h = v^t v / 2 = -1 / \beta
			A(first, pos) = f - g;
			if (!Symm) {
				// HH - Matrix von links multiplizieren
				for (int j = pos + 1; j < n; j++) {
					double s{ 0 };
					for (int i = first; i <= last; i++) s += A(i, pos) * A(i, j);
					f = s / h;
					for (int i = first; i <= last; i++) A(i, j) += f * A(i, pos);
				}
				// HH - Matrix von rechts multiplizieren
				for (int i = 0; i < n; i++) {
					double s{ 0 };
					for (int j = first; j <= last; j++) s += A(j, pos) * A(i, j);
					f = s / h;
					for (int j = first; j <= last; j++) A(i, j) += f * A(j, pos);
				}
			}
			else {
				// Berechne p = \beta D v.Hierbei ist D = A(i, j), i, j = first..last
				// p wird in u[k, k + 1], ..., u[k, n - 1] gespeichert.
				for (int i = first; i <= last; i++) {
					A(pos, i) = 0;
					for (int j = first; j <= last; j++) A(pos, i) += A(i, j) * A(j, pos);
					A(pos, i) = -A(pos, i) / h;
				}
				// Berechne w = p - (p^t v) / (v^t v) v
				// w überschreibt p.
				double s{ 0 };
				for (int i = first; i <= last; i++) s += A(pos, i) * A(i, pos);
				s /= 2;
				for (int i = first; i <= last; i++) A(pos, i) += A(i, pos) * s / h;
				// Update der Matrix D : D = D - vw^t - wv^t.Hierbei ist wieder
				// D = A(i, j), i, j = k + 1, ..., n - 1.
				for (int i = first; i <= last; i++)
					for (int j = first; j <= last; j++)
						A(i, j) -= (A(i, pos) * A(pos, j) + A(j, pos) * A(pos, i));
			}
			// Rückskalieren
			for (int i = first; i <= last; i++) A(i, pos) *= scale;
		}
		New = g * scale;
	}

/************************************************************************************************

	The Householder transformations that are stored in factored form in Q are applied to the
	matrix A{ *this }, which yields matrix B.
	A und B may coincide. In that case A is overwritten by B.
	It is assumed that Q.m = A.m (Transpose = false) resp. Q.n = A.m (Transpose = true).
	The routine computes
	  B = H_k H_{k-1} *...* H_2 * H_1 * A
	If 'Transpose=false' then the Householder vectors are stored as a lower triangular matrix in Q.
	If 'Transpose=true' then the Householder vectors are stored as an upper triangular matrix in Q.

	Parameters:
		A{*this} matrix to which the Householder tansformations are applied to.
		w        vector of diagonal values of R, where R is the upper triangular matrix of a
				 QR-Factorization; w is returned from function 'QR_Factporization'
		Q        matrix in which the Householder vectors are stored
		B		 matrix B = H_k H_{k-1} *...* H_2 * H_1 * A

*************************************************************************************************/

	void cMatrix::ApplyHH(double w[], const cMatrix& Q, cMatrix& B, bool Transpose) const {
		B = *this;
		int mn{ min(Q.m, Q.n) };
		if (!Transpose)
			// apply HH - transformations to the RHS b
			for (int k = 0; k < mn; k++) {
				if (w[k] != 0) {
					for (int j = 0; j < n; j++) {
						double s{ 0 };
						for (int i = k; i < m; i++) s += Q(i, k) * B(i, j);
						double f{ (s / Q(k, k)) / w[k] };
						for (int i = k; i < m; i++) B(i, j) += f * Q(i, k);
					}
				}
			}
		else
			// apply HH - transformations to the RHS b
			for (int k = 0; k < mn; k++){
				if (w[k] != 0){
					for (int j = 0; j < n; j++) {
						double s{ 0 };
						for (int i = k; i < m; i++) s += Q(k, i) * B(i, j);
						double f{ (s / Q(k, k)) / w[k] };
						for (int i = k; i < m; i++) B(i, j) += f * Q(k, i);
					}
				}
			}
	}

/************************************************************************************************/

	void ModifiedHouseholder(const double* x, int n, double* w, double* y, double& New) {
		double r{ 0 };
		New = 0;
		for (int i = 0; i < n; i++) r += x[i] * x[i];
		if (r != 0) {
			New = -copysign(sqrt(r), x[0]);     // g = -sign(A(first, pos)) * || a ||
			w[0] = x[0] - New;
			for (int i = 1; i < n; i++) {
				w[i] = x[i] / New;
				y[i] = x[i] / w[0];
			}
			w[0] /= New;
			y[0] = 1.0;
		}
	}

/************************************************************************************************/

	void cMatrix::ApplyModifiedHouseholder(const double w[], const double y[], int first, int last, int pos, matrix_direction Direction) {
		cMatrix& A{ *this };
		if (((int)Direction & (int)matrix_direction::mLeft) != 0)
			// HH - Matrix von links multiplizieren
			for (int j = pos; j < n; j++) {
				double s{ A(first, j) };
				for (int i = first + 1; i <= last; i++) s += y[i - first] * A(i, j);
				for (int i = first; i <= last; i++) A(i, j) += s * w[i - first];
			}
		if (((int)Direction & (int)matrix_direction::mRight) != 0)
			// HH - Matrix von rechts multiplizieren
			for (int i = 0; i < n; i++) {
				double s{ A(i, first) };
				for (int j = first + 1; j <= last; j++) s += y[j - first] * A(i, j);
				for (int j = first; j <= last; j++) A(i, j) += s * w[j - first];
			}
	}

/************************************************************************************************

	Calculates the values c, s of a Givens rotation s.t.

		| c   s | *| x | = | * |
		| -s  c |  | z |   | 0 |

	The value * is returned in r.

************************************************************************************************/

	void GivensRotation(double& c, double& s, double x, double z, double& r) {
		if (z == 0) {
			c = 1;
			s = 0;
			r = x;
		}
		else if (x == 0) {
			c = 0;
			s = 1;
			r = z;
		}
		else {
			r = hypot(x, z);
			c = x / r;
			s = z / r;
			if ((abs(x) > abs(z)) && (c < 0)) {
				c = -c;
				s = -s;
				r = -r;
			}
		}
	}

/************************************************************************************************

	The Givens rotation that is given by c and s is applied to the given matrix A{ *this }.
	The rows of A that are affected by the rotation are k1 and k2.
	the flags 'Direction' control if the Givens rotation is pre- or post-multiplied.

		mLeft: |  c  s | * A   mRight : A * | c -s |
			   | -s  c |					| s  c |

	Note, that both flags may be set.

************************************************************************************************/

	void cMatrix::ApplyGivens(double c, double s, int k1, int k2, matrix_direction Direction) {
		cMatrix& A{ *this };
		if (((int)Direction & (int)matrix_direction::mLeft) != 0)
			for (int j = 0; j < n; j++) {
				double x{ A(k1, j) };
				double y{ A(k2, j) };
				A(k1, j) =  c * x + s * y;
				A(k2, j) = -s * x + c * y;
			}
		if (((int)Direction & (int)matrix_direction::mRight) != 0)
			for (int i = 0; i < n; i++) {
				double x{ A(i, k1) };
				double y{ A(i, k2) };
				A(i, k1) =  c * x + s * y;
				A(i, k2) = -s * x + c * y;
			}
	}

/************************************************************************************************

	The Householder transformations that are stored in factored form in U{ *this} are
	accumulated in the matrix Q. Thus, Q is an orthogonal matrix that is equivalent to the
	Householder transformations stored in factored form in U{ *this}.

	In fact, only columns [c1,c2) of Q are computed.

	U und Q may coincide if 'transpose=false' and 'c1 = 0' and 'c2 = Q.n' and U.m >= U.n.
	In that case U is overwritten by Q.

	The routine computes columns j, c1 <= j < c2, of
	  Q = H_k H_{k-1} *...* H_2 * H_1 * I
	If 'Transpose=false' then the Householder vectors are stored as a lower triangular matrix in U.
	If 'Transpose=true' then the Householder vectors are stored as an upper triangular matrix in U.

	Parameters:
		U{*this} matrix in which the Householder vectors are stored
		D        vector of diagonal values of R, where R is the upper triangular matrix of a
				 QR-Factorization; D is returned from function 'QR_Factporization'
		ofs		 column (or row) of U that contains the first Householder vector
		nh		 number of Householder vectors stored in U
		Q		 matrix Q = H_k H_{k-1} *...* H_2 * H_1 * I
		c1		 first column of Q to be computed
		c2		 last column of Q to be computed

*************************************************************************************************/

	void cMatrix::AccumulateHH(double D[], int Ofs, int nh, cMatrix& Q, bool Transpose, int c1, int c2) const {
		if ((&Q == this) && !Transpose && (c1 == 0) && (this->m >= this->n) && (c2 == this->n)) {
			AccumulateHH(D, Ofs, nh);
			return;
		}
		const cMatrix& U{ *this };
		if (!Transpose) {
			if (c2 == -1) c2 = m;
			Q.SetSize(m, c2-c1);
			for (int i = c1; i < c2; i++) Q(i, i-c1) = 1;
			for (int k = nh - 1; k >= 0; k--)
				if (D[k + Ofs] != 0) {
					// Spalten 2 bis n berechnen
					for (int j = max((Ofs + k) + 1,c1); j < min(m,c2); j++) {
						double s{ 0 };
						for (int i = (Ofs + k) + 1; i < m; i++) s += U(i, k) * Q(i, j-c1);
						double f{ (s / U(Ofs + k, k)) / D[k + Ofs] };
						for (int i = (Ofs + k); i < m; i++) Q(i, j-c1) += f * U(i, k);
					}
					// Spalte 1 berechnen
					if ((Ofs + k >= c1) && (Ofs + k < c2)) {
						for (int i = Ofs + k; i < m; i++) Q(i, Ofs + k - c1) = U(i, k) / D[k + Ofs];
						// Und noch eine Eins auf der Diagonalen addieren
						Q(Ofs + k, Ofs + k - c1) += 1;
					}
  			}
		}
		else {
			if (c2 == -1) c2 = n;
			Q.SetSize(n, c2 - c1);
			for (int i = c1; i < c2; i++) Q(i, i-c1) = 1;
			for (int k = nh - 1; k >= 0; k--)
				if (D[k + Ofs] != 0) {
					// Spalten 2 bis n berechnen
					for (int j = max((Ofs + k) + 1, c1); j < min(n, c2); j++) {
						double s{ 0 };
						for (int i = (Ofs + k) + 1; i < n; i++) s += U(k, i) * Q(i, j-c1);
						double f{ (s / U(k, Ofs + k)) / D[k + Ofs] };
						for (int i = Ofs + k; i < n; i++) Q(i, j-c1) += f * U(k, i);
					}
					// Spalte 1 berechnen
					if ((Ofs + k >= c1) && (Ofs + k < c2)) {
						for (int i = Ofs + k; i < n; i++) Q(i, Ofs + k - c1) = U(k, i) / D[Ofs + k];
						// Und noch eine Eins auf der Diagonalen addieren
						Q(Ofs + k, Ofs + k - c1) += 1;
					}
				}
		}
	}

/************************************************************************************************

	The Householder transformations that are stored in factored form in U{ *this} are
	accumulated in the matrix Q. Thus, Q is an orthogonal matrix that is equivalent to the
	Householder tarnsformations stored in in factored form in U{ *this}.
	U is overwritten by Q.

	Parameters:
		U{*this} matrix in which the Householder vectors are initialyy stored.
				 U is overwritten by H_k H_{k-1} *...* H_2 * H_1 * I
		D        vector of diagonal values of R, where R is the upper triangular matrix of a
				 QR-Factorization; D is returned from function 'QR_Factporization'
		ofs		 column (or row) of U that contains the first Householder vector
		nh		 number of Householder vectors stored in U

*************************************************************************************************/

	void cMatrix::AccumulateHH(double D[], int Ofs, int nh) const {
		const cMatrix& U{ *this };
		for (int i = Ofs + nh; i < m; i++)
			for (int j = Ofs + nh; j < n; j++) U(i, j) = (i == j) ? 1 : 0;
		for (int k = nh - 1; k >= 0; k--) {
			// Erste Zeile (bis auf erste Spalte) gleich Null setzen
			for (int j = (Ofs + k) + 1; j < n; j++) U(Ofs + k, j) = 0;
			if (D[k + Ofs] != 0) {
				// Spalten 2 bis n berechnen
				for (int j = (Ofs + k) + 1; j < n; j++) {
					double s{ 0 };
					for (int i = (Ofs + k) + 1; i < m; i++) s += U(i, k) * U(i, j);
					double f{ (s / U(Ofs + k, k)) / D[k + Ofs] };
					for (int i = (Ofs + k); i < m; i++) U(i, j) += f * U(i, k);
				}
				// Spalte 1 berechnen
				for (int i = Ofs + k; i < m; i++) U(i, Ofs + k) = U(i, k) / D[k + Ofs];
			}
			else {
				// Spalte 1 berechnen
				for (int i = Ofs + k; i < m; i++) U(i, Ofs + k) = 0;
			}
			// Und noch eine Eins auf der Diagonalen addieren
			U(Ofs + k, Ofs + k) += 1;
		}
		for (int i = 0; i < Ofs; i++) {
			U(i, i) = 1;
			for (int j = i + 1; j < m; j++) U(i, j) = 0;
			for (int j = i + 1; j < m; j++) U(j, i) = 0;
		}
	}

/************************************************************************************************

	Computes the QR-Factorization with column pivoting of the given matrix:
		AP    = QR (transpose = false)
		A^T P = QR (transpose = true)
	Here,
		Q is an orthogonal matrix.
		R is an upper triangular matrix
		P describes the column permutation resulting from pivoting.

	The orthogonal matrix Q is stored in factored form in the
		lower triangular part of QR (trasnpose = false) or in the
		upper triangular part of QR (trasnpose = true).
	The matrix Q can be explicitly constructed by calling the routine 'AccumulateHH'.

	If 'transpose = false', R is stored in the upper triangular part of QR.
	If 'transpose = true', R^T is stored in the lower triangular part of QR.
	In any case the main diagonal of R is not stored in QR but in the vector D.

	Parameters:
		A{*this} matrix to be factored.
		QR       stores matrices Q and R as described above, QR may coincide with A, in which
				 case A is overwritten by QR.
		D		 stores the main diagonal of R.
		p        stores the column permutation.
		rank	 rank of A

************************************************************************************************/

	void cMatrix::QR_Factorization(cMatrix& QR, double D[], int p[], bool Transpose, int* rank) const {
		const double eps2{ 1e-14 };
		QR = { *this };
		if (rank) *rank = 0;
		int mn{ min(m, n) };
		double fn0{ 0 };
		if (!Transpose)  {
			unique_ptr<double[]> gamma{ new double[n]{} };
			for (int j = 0; j < n; j++) p[j] = j;
			// Compute euclidean norm of columns
			for (int j = 0; j < n; j++) {
				for (int i = 0; i < m; i++) gamma[j] += QR(i, j) * QR(i, j);
				fn0 += gamma[j];
			}
			for (int k = 0; k < mn; k++) {
				// Find column of maximum norm
				double Max{};
				double fn{ 0 };
				int l{ k };
				for (int j = k; j < n; j++) {
					fn += gamma[j];
					if (gamma[j] > Max)  {
						Max = gamma[j];
						l = j;
					}
				}
				if (fn < eps2 * fn0) break;
				if (l != k)  {
					// If necessary, exchange columns.
					swap(p[k], p[l]);
					swap(gamma[k], gamma[l]);
					QR.XChangeCols(k, l);
				}
				if (rank) (*rank)++;
				// Householder transformation
				QR.Householder1(D[k], k, m - 1, k);
				// Update column norms
				for (int j = k + 1; j < n; j++) gamma[j] -= QR(k, j) * QR(k, j);
			}
		}
		else {
			unique_ptr<double[]> gamma{ new double[m]{} };
			for (int i = 0; i < m; i++) p[i] = i;
			// Compute euclidean norm of rows
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) gamma[i] += QR(i, j) * QR(i, j);
				fn0 += gamma[i];
			}
			for (int k = 0; k < mn; k++) {
				// Find row of maximum norm
				double Max{};
				double fn{ 0 };
				int l{ k };
				for (int i = k; i < m; i++) {
					fn += gamma[i];
					if (gamma[i] > Max) {
						Max = gamma[i];
						l = i;
					}
				}
				if (fn < eps2 * fn0) break;
				if (l != k)  {
					// If necessary, exchange rows.
					swap(p[k], p[l]);
					swap(gamma[k], gamma[l]);
					QR.XChangeRows(k, l);
				}
				if (rank) (*rank)++;
				// Householder transformation
				QR.Householder1(D[k], k, n - 1, k, true);
				// Update row norms
				for (int i = k + 1; i < m; i++) gamma[i] -= QR(i, k) * QR(i, k);
			}
		}
	}

/************************************************************************************************

	Bidiagonalization via Householder transformations of the given matrix A{ *this }.
	A is an m x n - matrix. The routine determines orthogonal matrices U and V s.t.
		U^T * A   * V = BD, if m >= n or
		U^T * A^T * V = BD, if m < n.
	The orthogonal matrices U and V are computed only if U (V) != nullptr.
	The main diagonal of BD is stored in HD, the upper diagonal is stored in ND.

	Parameters:
		A{*this} matrix to be bidiagonalized.
		HD		main diagonal of BD. The elements are stored at the index positions 0..n-1.
		ND		upper diagonal of BD. The elements are stored at the index positions 1..n-1.
		U		pointer to matrix U. If U != nullptr matrix U is returned in *U.
				+U may coincide with A.
		V		pointer to matrix V. If V != nullptr matrix V is returned in *V.

************************************************************************************************/

	void cMatrix::BiDiag(double HD[], double ND[], cMatrix* U, cMatrix* V) const {
		// Bidiagonalizaton via Householder transformations
		cMatrix A{ *this };
		ND[0] = 0;
		int nn{ min(m, n) };
		bool mgen{ m >= n }; // true, if m >= n
		for (int k = 0; k < nn; k++) {
			A.Householder1(HD[k], k, max(m, n) - 1, k, !mgen);
			if (k < nn - 1) A.Householder1(ND[k + 1], k + 1, nn - 1, k, mgen);
		}
		// If necessary, compute U and V
		if (V) A.AccumulateHH(ND, 1, nn - 1, *V, mgen);
		if (U) A.AccumulateHH(HD, 0, nn, *U, !mgen);
	}

/************************************************************************************************

	Tridiagonalization via Householder transformations of the given symmetric matrix A{ *this }.
	A is a symmetric square n x n - matrix. The routine determines an orthogonal matrix U s.t.
		U^T * A * U = TD,
	where TD is a tridiagonal matrix.
	The orthogonal matrix U is computed only if U != nullptr.
	The main diagonal of TD is stored in HD, the sub- (or super-) diagonal is stored in ND.

	Parameters:
		A{*this} matrix to be tridiagonalized.
		HD		main diagonal of TD. The elements are stored at the index positions 0..n-1.
		ND		subdiagonal of TD. The elements are stored at the index positions 1..n-1.
		U		pointer to matrix U, if U != nullptr matrix U is returned in *U.
				*U may coincide with A.

************************************************************************************************/

	void cMatrix::TriDiag(double HD[], double ND[], cMatrix* U) const {
		// Tridiagonalizatopn via Householder transformations
		cMatrix* V{ (U != nullptr) ? U : new cMatrix };
		cMatrix& A{ *V };
		A = *this;

		ND[0] = 0;
		for (int k = 0; k < n - 2; k++) {
			// Hauptdiagonalelemente können direkt übernommen werden.
			HD[k] = A(k, k);
			A.Householder2(ND[k + 1], k + 1, n - 1, k, true);
		}
		// Restliche Diagonalelemente auffüllen
		HD[n - 2] = A(n - 2, n - 2);
		ND[n - 1] = A(n - 1, n - 2);
		HD[n - 1] = A(n - 1, n - 1);
		// If necessary, compute U
		if (U) A.AccumulateHH(ND, 1, n - 2); else delete V;
	}

/************************************************************************************************

	Transforms the square matrix A{ *this } to Hessenberg form using Householder transformations.
	A is a square n x n - matrix. The routine determines an orthogonal matrix U s.t.
		U^T * A * U = H,
	where H is an upper Hessenberg matrix.
	The orthogonal matrix U is computed only if U != nullptr.

	Parameters:
		A{*this} matrix to be transformed to Hessenberg form.
		H		stores the Hessenberg matrix H. H may coincide with A, in which
				case A is overwritten by H.
		U		pointer to matrix U. If U != nullptr matrix U is returned in *U.
				*U may coincide with A.

************************************************************************************************/

	void cMatrix::Hessenberg(cMatrix& H, cMatrix* U) const {
		H = *this;
		unique_ptr<double[]> SD{ new double[n] };
		SD[0] = 0;
		for (int k = 0; k < n - 1; k++) H.Householder2(SD[k + 1], k + 1, n - 1, k, false);
		if (U) H.AccumulateHH(SD.get(), 1, n - 1, *U);
		for (int j = 0; j < n - 1; j++) {
			H(j + 1, j) = SD[j + 1];
			for (int i = j + 2; i < n; i++) H(i, j) = 0;
		}
	}

/************************************************************************************************

	Solves the linear equation L*U*x = P*B.
	L,U and P represent the LU factorization of a square matrix, that is computed by a call
	of the routine 'LU_Factorization'.
	A is a sqaure n n x n matrix, B an n x r matrix.

	Parameters:
		LU{ *this } is the LU factorization of a square matrix A.
		P		describes the row permutations made during the Lu-factorization
		B		matrix of right hand sides
		X		solution of the system A * X = B. Xmay coincide with B in which case B is
				overwritten by X.

************************************************************************************************/

	void cMatrix::LGS_Solve_LU(int p[], const cMatrix& B, cMatrix& X) {
		cMatrix& LU{ *this };
		int r{ B.n };
		// Copy B to X
		X = B;
		// Permutation of the rows of X
		for (int i = 0; i < n; i++)
			if (p[i] != i) X.XChangeRows(i, p[i]);
		// Loop through all r right hand sides
		for (int k = 0; k < r; k++) {
			// Solve L*y = b
			for (int i = 0; i < n; i++) {
				for (int j = 0; j < i; j++) X(i, k) += LU(i, j) * X(j, k);
				X(i, k) = -X(i, k);
			}
			// Solve U*x = y
			for (int i = n - 1; i >= 0; i--) {
				for (int j = i + 1; j < n; j++) X(i, k) += LU(i, j) * X(j, k);
				X(i, k) /= -LU(i, i);
			}
		}
	}

/************************************************************************************************

	Solves the linear equation A * x = B via LR factorization.
	A is a sqaure n n x n matrix, B an n x r matrix.

	Parameters:
		A		square n x n matrix
		B		matrix of right hand sides
	Returns:
		X		solution of the system A * X = B

************************************************************************************************/

	cMatrix LGS_Solve_LU(const cMatrix& A, const cMatrix& B) {
		unique_ptr<int[]> p{ new int[A.n] };
		cMatrix LU(A.n, A.n);
		cMatrix X(B.m, B.n);
		double det{};
		A.LU_Factorization(LU, p.get(), det);
		LU.LGS_Solve_LU(p.get(), B, X);
		return X;
	}

/************************************************************************************************

	Solves the linear equation Q * R * P^T * X = B.
	Q,R and P are the result of a QR-factorization of a matrix A (transpose=false)
	                                                        or A^T (transpose=true).
	In general, the QR-factorization yields
	  AP = Q [ R11 R12 ] ,
	         [  0   0  ]
	where R11 is a regular upper triangular matrix.
	Let PX = [ X1 ] and Q^T * B = [ C1 ] .
	         [ X2 ]               [ C2 ]
	The routine returns X s.t. X1 = R11^{-1} * C1, X2 = 0.

	If m > n and rank = n, this is the unique LS solution to the overdetermined system A x = B.
	If m > n and rank < n, there is no unique LS solution. Then, the returned X is the so-called
		basic solution to the LS problem (cf. Golub, Reinsch, Matrix Computations)
	If m = n and rank = m, this is the unique solution of the linear equation A x = B.
	If m < n and rank = m, there is an infinite number of solutions of A x = b. Then, the
		returned X is a so-called basic solution.
	If m <= n and rank < m, then the system may have no solution. In that case a basic
		LS solution is returned.

	If 'Transpose = false', then QR is the QR-factorization of A, A P = Q R. In this case, R is
	stored in factored form in the lower triangular part of QR and R in the upper traingular part.
	If 'Transpose = true', then QR is the QR-factorization of A^T, A^T P = Q R. In this case, R is
	stored in factored form in the upper triangular part of QR and R^T in the lower traingular part.

	Parameters:
		QR{ *this } is the QR factorization of a matrix A (or A^T).
		D		vector of diagonal values of R, where R is the upper triangular matrix of a
				QR-Factorization; D is returned from function 'QR_Factporization'
		P		describes the column permutation resulting from pivoting.
		rank	Rank of R (or, equivalently, A).
		B		matrix of right hand sides
		X		solution of the system Q * R * P^T * X = B.
		transpose If Transpoe = false

************************************************************************************************/

	void cMatrix::LGS_Solve_QR(double D[], int p[], int rank, const cMatrix& B, cMatrix& X, bool transpose) {
		cMatrix& QR{ *this };
		int r{ B.n };
		if (!transpose) {
			cMatrix QtB(m, r);
			X.SetSize(n, r);
			// Rechte Seite von links mit Q multiplizieren
			B.ApplyHH(D, QR, QtB, false);
			// gestaffeltes Gleichungssystem lösen
			for (int j = 0; j < r; j++) {
				for (int k = rank - 1; k >= 0; k--) {
					X(p[k], j) = QtB(k, j);
					for (int i = k + 1; i < rank; i++) X(p[k], j) -= X(p[i], j) * QR(k, i);
					X(p[k], j) /= D[k];
				}
			}
		}
		else {
			cMatrix QtB(n, r);
			X.SetSize(m, r);
			// Rechte Seite von links mit Q multiplizieren
			B.ApplyHH(D, QR, QtB, true);
			// gestaffeltes Gleichungssystem lösen
			for (int j = 0; j < r; j++) {
				for (int k = rank - 1; k >= 0; k--) {
					X(p[k], j) = QtB(k, j);
					for (int i = k + 1; i < rank; i++) X(p[k], j) -= X(p[i], j) * QR(i, k);
					X(p[k], j) /= D[k];
				}
			}
		}
	}

/************************************************************************************************

	Solves the linear equation A * x = B via QR factorization.
	A is a square n x n matrix, B an n x r matrix.
	For a detailed description please refer to the previous overload of LGS_Solve_QR above.

	Parameters:
		A	m x n matrix
		B	m x r matrix of right hand sides
	Returns:
		X	solution of the system A * X = B

************************************************************************************************/

	cMatrix LGS_Solve_QR(const cMatrix& A, const cMatrix& B) {
		cMatrix QR(A.m, A.n);
		cMatrix X(B.m, B.n);
		int rank{};
		unique_ptr<double[]> D{ new double[A.n] };
		unique_ptr<int[]> p{ new int[A.n] };
		A.QR_Factorization(QR, D.get(), p.get(), false, &rank);
		QR.LGS_Solve_QR(D.get(), p.get(), rank, B, X);
		return X;
	}

/************************************************************************************************

	Computes a  basis of the kernel (nullspace) of a given m x n matrix A{ *this } via
	QR-factorization.

	Parameters:
		A{ *this } m x n matrix of which the nullspace has to be computed
	Returns:
		X	basis of the kernel of A.

************************************************************************************************/

	cMatrix cMatrix::Nullspace_QR() {
		cMatrix QR, Q;
		unique_ptr<double[]> D{ new double[n] };
		unique_ptr<int[]> p{ new int[n] };
		int rank{};
		QR_Factorization(QR, D.get(), p.get(), true, &rank);
		QR.AccumulateHH(D.get(), 0, rank, Q, true, rank, n);
		return Q;
	}

/************************************************************************************************

	Computes a  basis of the orthogonal complement of the subspace spanned by the columns of
	a given m x n matrix A{ *this } via QR-factorization.

	Parameters:
		A{ *this } m x n matrix.
	Returns:
		X	basis of the orthogonal complement of the subspace spanned by the columns of A.

************************************************************************************************/

	cMatrix cMatrix::OrthogonalComplement() {
		cMatrix QR, Q;
		unique_ptr<double[]> D{ new double[n] };
		unique_ptr<int[]> p{ new int[n] };
		int rank{};
		QR_Factorization(QR, D.get(), p.get(), false, &rank);
		QR.AccumulateHH(D.get(), 0, rank, Q, false, rank, m);
		return Q;
	}

/************************************************************************************************

	Solves the linear equation A x = b where A is a symmetric positive definite matrix.
	In a first step the Cholesky decomposition A = L L^T is computed. In the second step
	the linear equation L L^T x = b is solved.

	Parameters:
		A{ *this } symmetric positive definite n x n matrix.
		b		matrix of right hand sides
		pRank	if pRank != nullptr, the rank of A is returned in *pRank.
	Returns:
		X		solution of the system A * X = B

************************************************************************************************/

	cMatrix LGS_Solve_PosDef(const cMatrix& A, const cMatrix& B, int* pRank) {
		cMatrix G{ A };
		cMatrix X(B.m, B.n);
		int tmp{};
		int& rank{ (pRank != nullptr) ? *pRank : tmp };
		G.Cholesky(G, &rank);
		if (rank < A.n) throw runtime_error("Matrix singular!");
		X = B;
		// Lösen des LGS G*x1 = X^ty
		for (int k = 0; k < X.n; k++) {
			for (int i = 0; i < A.n; i++) {
				for (int j = 0; j < i; j++) X(i, k) -= G(i, j) * X(j, k);
				X(i, k) /= G(i, i);
			};
			// Lösen des LGS G^t*x2 = x1
			for (int i = A.n - 1; i >= 0; i--) {
				for (int j = i + 1; j < A.n; j++) X(i, k) -= G(j, i) * X(j, k);
				X(i, k) /= G(i, i);
			}
		};
		return X;
	};

/************************************************************************************************

	Computes the (pseudo)inverse of a positive semidefinite matrix A{ *this }.
	Only the lower triangular part of A is used.
	The inverse is returned in B. B may coincide with A in which case A is overwritten by B.

	If A has not full rank a pseudo inverse of A is computed. More specific: If rank(A) = r < n,
	then A can be represented (apart from permutaions of the rows and columns) as
		A = | A_11 A_12 |, where A_11 is an r x r matrix of full rank.
			| A_21 A_22 |
	In that case the returned pseudo inverse is given by
		A^-= | A_11^-1  0 | .
			 | 0		0 |

	Parameters:
		A{ *this } symmetric positive (semi)definite n x n matrix.
		B		(pseudo)inverse of A
		pRank	if pRank != nullptr, the rank of A is returned in *pRank.

************************************************************************************************/

	void cMatrix::InversePosDef(cMatrix& B, int* pRank) const {
		B = *this;
		B.Cholesky(B, pRank);
		for (int j = 0; j < n; j++) {
			if (B(j, j) > 0) {
				B(j, j) = 1 / B(j, j);
				for (int i = j + 1; i < n; i++) {
					if (B(i, i) > 0) {
						double Sum{ };
						for (int k = j; k < i; k++) Sum -= B(i, k) * B(k, j);
						B(i, j) = Sum / B(i, i);
					}
					else B(i, j) = 0;
				}
			}
		}
		for (int i = 0; i < n; i++)
			for (int j = 0; j <= i; j++) {
				double Sum{ };
				for (int k = i; k < n; k++) Sum += B(k, i) * B(k, j);
				B(i, j) = Sum;
				// Inverse oberhalb der Diagonalen auffüllen
				if (i != j) B(j, i) = B(i, j);
			}
	}

/************************************************************************************************

	Computes the (pseudo)inverse of a positive semidefinite matrix A{ *this }.
	For a detailed description please refer to cMatrix::InversePosDef above.

	Parameters:
		A		symmetric positive (semi)definite n x n matrix.
		pRank	if pRank != nullptr, the rank of A is returned in *pRank.
	Returns:
		B		(pseudo)inverse of A

************************************************************************************************/

	cMatrix InversePosDef(const cMatrix& A, int* pRank) {
		cMatrix B;
		A.InversePosDef(B, pRank);
		return B;
	}

/************************************************************************************************

	Computes the inverse of a given matrix A{ *this } with the Gauß-Jordan algorithm.
	A is overwritten by the inverse of A.

	Parameters:
		A{ *this } matrix to be inverted.

************************************************************************************************/

	cMatrix& cMatrix::Invert() {
		cMatrix& A{ *this };
		const double eps = 1e-15;
		assert(m == n);
		unique_ptr<int[]> p{ new int[n] };
		for (int j = 0; j < n; j++) p[j] = j;
		for (int j = 0; j < n; j++) {
			// Find pivot
			int r{ j };
			double max{ abs(A(j, j)) };
			for (int i = j + 1; i < n; i++) {
				if (abs(A(i, j)) > max) {
					max = abs(A(i, j));
					r = i;
				}
			}
			if (max < eps) throw runtime_error("Matrix singular!");
			// Change rows
			if (r > j) {
				A.XChangeRows(r, j);
				swap(p[j], p[r]);
			}
			// Transformation
			double hr = 1.0 / A(j, j);
			for (int i = 0; i < n; i++) A(i, j) *= -hr;
			A(j, j) = hr;
			for (int k = 0; k < n; k++)
				if (k != j) {
					for (int i = 0; i < n; i++)
						if (i != j) A(i, k) += A(i, j) * A(j, k);
					A(j, k) *= hr;
				}
		}
		// Change columns
		unique_ptr<double[]> tmp{ new double[n] };
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < n; k++) tmp[p[k]] = A(i, k);
			for (int k = 0; k < n; k++) A(i, k) = tmp[k];
		}
		return A;
	}

/************************************************************************************************

	Computes the inverse of a given matrix A with the Gauß-Jordan algorithm.

	Parameters:
		A	matrix to be inverted.
	Returns:
		B	inverse of A.

************************************************************************************************/

	cMatrix Inverse(cMatrix A) {
		cMatrix B{ A };
		return B.Invert();
	}

/************************************************************************************************

	Computes the reduced row ecehelon form of a given matrix A{ *this } with the Gauß-Jordan
	algorithm.
	A is overwritten by the reduced row ecehelon form of A.

	A matrix is in reduced row echelon form if
		*	all nonzero rows are above any zero rows,
		*	the leading coefficient of a nonzero row is always strictly to the right of the
			leading coefficient of the row above it,
		*	Every leading coefficient is 1 and is the only nonzero entry in its column.
	Note, that the reduced row ecehelon form of a matrix is unique.

	Parameters:
		A{ *this } matrix of which the reduced row echelon form has to be computed.

************************************************************************************************/

	cMatrix& cMatrix::Rref() {
		const double eps = 1e-14;
		cMatrix& A{ *this };
		// Zeilenskalierung
		unique_ptr<double[]> scale{ new double[m]{} };
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) scale[i] += abs(A(i, j));
			scale[i] = 1 / scale[i];
		}
		int j{-1};
		for (int k = 0; k < m; k++) {
			// Find pivot
			double max;
			int r;
			do {
				j++;
				r = k;
				max = abs(A(k, j)) * scale[k];
				for (int i = k + 1; i < m; i++) {
					double tmp{ abs(A(i, j)) * scale[i] };
					if (tmp > max) {
						max = tmp;
						r = i;
					}
				}
			} while ((max < eps) && (j < n - 1));
			if (max > eps) {
				// Change rows
				if (r != k) {
					A.XChangeRows(r, j);
					swap(scale[r], scale[j]);
				}
				// Calculate pivot row
				for (int l = j + 1; l < n; l++) A(k, l) /= A(k, j);
				A(k, j) = 1;
				// Calculate remaining rows
				for (int i = 0; i < k; i++) {
					for (int l = j + 1; l < n; l++) A(i, l) -= A(i, j) * A(k, l);
					A(i, j) = 0;
				}
				for (int i = k + 1; i < m; i++) {
					for (int l = j + 1; l < n; l++) A(i, l) -= A(i, j) * A(k, l);
					A(i, j) = 0;
				}
			}

		}
		return A;
	}

/************************************************************************************************

	Computes all eigenvalues and if necessary all eigenvetors of a symmetric matrix A{ *this }.
	as well as the decomposition  A = U diag(EW) U^t.
	This is done by applying the implicit QR algorithm with Wilkinson shift.

	Parameters:
		A{ *this } matrix whose eigenvalues have to be computed. A has to be symjmetric.
				A is not changed by this routine.
		EW		The eigenvalues of A are stored in the array EW. The eigenvalues are not
				sorted.
		U		pointer to matrix U. If U != nullptr matrix U is returned in *U.
				*U may coincide with A in which case A is overwritten by U.
				The columns of U are the eigenvectors of A.
	Returns:
		k		error code. k == 0 if no error occurred.
				If k != 0, then the (k-1)-th eigenvalue could not be determined within 30
				iterations. In that case the eigenvalues at the index position k,...,n-1
				as well as the columns of U with indices k,...,n-1 should be correct.

************************************************************************************************/

	int cMatrix::EigenSymm(double EW[], cMatrix* U) {
		cMatrix& A{ *this };
		unique_ptr<double[]> ND{ new double[n] };
		A.TriDiag(EW, ND.get(), U);
		double tst1{ 0 };
		for (int i = 0; i < n; i++) tst1 = max(tst1, abs(EW[i]) + abs(ND[i]));
		for (int k = n - 1; k >= 0; k--) {
			int its{};
			int l{};
			do {
				// test for splitting.
				l = k;
				while (tst1 + abs(ND[l]) != tst1) l--;
				// test for convergence
				if ((l != k) && (its < 30)) {
					// shift from bottom 2 by 2 minor
					its++;
					double x{ (EW[k - 1] - EW[k]) / (2 * ND[k]) };
					x = EW[l] - EW[k] + ND[k] / (x + copysign(hypot(x, 1), x));
					double s{ 1 };
					double c{ 1 };
					double p{ 0 };
					for (int i = l; i < k; i++) {
						double z{ s * ND[i + 1] };
						double b{ c * ND[i + 1] };
						double r;
						GivensRotation(c, s, x, z, r);
						if (i > l) ND[i] = r;
						x = EW[i] - p;
						r = (EW[i + 1] - x) * s + 2 * c * b;
						p = s * r;
						EW[i] = x + p;
						x = c * r - b;
						if (U) U->ApplyGivens(c, s, i, i + 1, matrix_direction::mRight);
					}
					EW[k] = EW[k] - p;
					ND[k] = x;
				}
			} while ((l != k) && (its < 30));
			if (its >= 30)	return k+1;					// set error -- no convergence to an eigenvalue after 30 iterations
		}
		return 0;
	}

/************************************************************************************************

	Computes the real Schur decomposition of a square matrix A{ *this }, i.e., A = U H U^T,
	as well as all eigenvalues of A. Here, A is an (n,n) matrix, U an orthogonal (n,n) matrix
	and H an upper quasi-triangular matrix. The eigenvalues of A are computed from the (1,1)-
	and (2,2)-blocks on the diagonal of H:

	Parameters:
		A{ *this } matrix whose eigenvalues have to be computed. A has to be a squre matrix.
				A is not changed by this routine.
		EW		The eigenvalues of A are stored in the array EW. The eigenvalues are not
				sorted.
		H		pointer to the upper quasi-triangular matrix H. If H != nullptr matrix H is
				returned in *H.
		U		pointer to the orthogonal matrix U. If U != nullptr matrix U is returned in *U.
				*U may coincide with A in which case A is overwritten by U.
	Returns:
		k		error code. k == 0 if no error occurred.
				If k != 0, then the (k-1)-th eigenvalue could not be determined within 30
				iterations. In that case the eigenvalues at the index position k,...,n-1
				as well as the columns of U with indices k,...,n-1 should be correct.

************************************************************************************************/

	int cMatrix::EigenQR(complex<double> EW[], cMatrix* _H, cMatrix* U) {
		const complex<double> _i{ 0, 1 };
		cMatrix& A{ *this };
		cMatrix& H{ (_H != nullptr) ? *_H : *(new cMatrix) };
		A.Hessenberg(H, U);
		double x{ abs(H(0, 0)) };
		for (int k = 1; k < n; k++) x = max(x, abs(H(k, k)) + abs(H(k, (k - 1))));
		double error{};
		int m{ n - 1 };
		int k{};
		while (m >= 0) {
			int its{};
			int l{};
			do {
				// test for (int splitting.
				l = m;
				while ((l >= 1) && ((abs(H(l - 1, l - 1)) + abs(H(l, l))) + abs(H(l, l - 1)) !=
					(abs(H(l - 1, l - 1)) + abs(H(l, l))))) l--;
				// test for (int convergence
				if (l == m)  {
					EW[m] = H(m, m);
					m--;
				}
				else if (l == m - 1)  {
					x = H(m, m);
					double y{ H(l, l) };
					double w{ H(l, m) * H(m, l) };
					double p{ 0.5 * (y - x) };
					double q{ p*p + w };
					double z{ sqrt(abs(q)) };
					if (q >= 0)  {
						z = p + copysign(z, p);
						EW[l] = x + z;
						EW[m] = (z != 0) ? x - w / z : x + z;
						if ((_H != nullptr) || (U != nullptr))  {
							double s, t;
							GivensRotation(s, t, z, H(m, l), x);
							if (_H) H.ApplyGivens(s, t, l, m, matrix_direction::mLeftRight);
							if (U) U->ApplyGivens(s, t, l, m, matrix_direction::mRight);
						}
					}
					else {
						EW[l] = complex<double>(x + p, z);//(x + p) + z * _i;
						EW[m] = complex<double>(x + p, -z);//(x + p) - z * _i;
					}
					m = m - 2;
				}
				else if ((l < m - 1) && (its < 30))  {
					its++;
					// Francis QR - Schritt(bearbeitet die Teilmatrix
					// H(i, j), i = l, ..., nn, j = l, ..., nn
					double s{ H(m - 1, m - 1) + H(m, m) };
					double t{ H(m - 1, m - 1) * H(m, m) - H(m - 1, m) * H(m, m - 1) };
					double xx[3];
					xx[0] = H(l, l) * H(l, l) + H(l, l + 1) * H(l + 1, l) - s * H(l, l) + t;
					xx[1] = H(l + 1, l) * (H(l, l) + H(l + 1, l + 1) - s);
					xx[2] = H(l + 1, l) * H(l + 2, l + 1);
					for (k = l; k < m; k++) {
						double ww[3], yy[3], New;
						ModifiedHouseholder(xx, min(3, m - (k - 1)), ww, yy, New);
						if (New != 0)  {
							H.ApplyModifiedHouseholder(ww, yy, k, min(k + 2, m), max(k - 1, l), matrix_direction::mLeftRight);
							if (U) U->ApplyModifiedHouseholder(ww, yy, k, min(k + 2, m), 0, matrix_direction::mRight);
						}
						xx[0] = H(k + 1, k);
						if (k < m - 1)  xx[1] = H(k + 2, k);
						if (k < m - 2)  xx[2] = H(k + 3, k);
					}
				}
			} while ((l < m - 1) && (its < 30));
			if (its >= 30) return k;
			// set error -- no convergence to an
			//   eigenvalue after 30 iterations
		}
		if (!_H) delete &H;
		return 0;
	}

/************************************************************************************************

	Computes the singular value decomposition of a (m,n)-matrix A{ *this }, i.e., A = U S V^T.
	Here, A is an (n,n)-matrix, U an orthogonal (m,m)-matrix, V an orthogonal (n,n)-matrix
	and S an (m,n)-matrix, whoe entries outside the main diagonal are all zero. The entries
	on the main diagonal of S are the so-called singular values. The algorithm is based on
	Householder	bidiagonalization and a variant of the QR algorithm.

	This routine is a translation of the algol procedure svd (Golub, G.H., and Reinsch, C.
	(1970). Singular Value Decomposition and Least Squares Solutions, Numerical Mathematics 14,
	403--420. resp. Handbook for Automatic Computation, Vol II: Linear Algebra, 134-151.)

	Parameters:
		A{ *this } is the matrix whose singualr value decomposition has to be computed.
				A is not changed by this routine.
		w		The (nonnegative) singular values of A are stored in the array w. The singualr
				values are not sorted.
		U		pointer to the orthogonal matrix U. If U != nullptr the matrix U is returned
				in *U. *U may coincide with A in which case A is overwritten by U.
		V		pointer to the orthogonal matrix V. If V != nullptr the matrix V is returned
				in *VU. If *U does not coincide with A, *V may coincide with A in which case
				A is overwritten by V.
	Returns:
		k		error code. k == 0 if no error occurred.
				If k != 0, then the (k-1)-th singular value could not be determined within 30
				iterations. In that case the singular values at the index position k,...,n-1
				as well as the corresponding columns of U and V should be correct.

************************************************************************************************/

	int cMatrix::SVD(double w[], cMatrix* _U, cMatrix* _V) {
		cMatrix& A{ *this };
		cMatrix*& U{ (m >= n) ? _U : _V };
		cMatrix*& V{ (m >= n) ? _V : _U };

		int nn{ min(m, n) };
		unique_ptr<double[]> ND{ new double[nn] };
		int ierr{};
		// Bidiagonalisierung durch Householder - Transformationen
		A.BiDiag(w, ND.get(), U, V);
		// Maximum von abs(w[i]) + abs(ND[i]) ermitteln
		double x{};
		for (int i = 0; i < nn; i++) x = max(x, abs(w[i]) + abs(ND[i]));
		// diagonalization of the bidiagonal form
		double tst1{ x };
		int its{};
		int k;
		for (k = nn - 1; k >= 0; k--) {
			its = 0;
			bool Convergence{ false };
			do {
				// test for splitting.
				bool Skip{ false };
				int l;
				double tst2;
				for (l = k; l >= 0; l--) {
					tst2 = tst1 + abs(ND[l]);
					if (tst2 == tst1) {
						Skip = true;
						break;
					}
					// ND(1) is always zero, so there is no exit
					// through the bottom of the loop
					tst2 = tst1 + abs(w[l - 1]);
					if (tst2 == tst1) break;
				}
				// cancellation of ND(l) if l greater than 1
				if (!Skip) {
					double c{ 0 };
					double s{ -1 };
					double f;
					for (int i = l; i <= k; i++) {
						f = -s * ND[i];
						ND[i] = c * ND[i];
						tst2 = tst1 + abs(f);
						if (tst2 == tst1) break;
						GivensRotation(c, s, w[i], f, w[i]);
						if (U) U->ApplyGivens(c, -s, l - 1, i, matrix_direction::mRight);
					}
				}
				// test for convergence
				double z{ w[k] };
				Convergence = (l == k);
				if (!Convergence && (its < 30)) {
					//shift from bottom 2 by 2 minor
					its++;
					double x{ w[l] };
					double y{ w[k - 1] };
					double g{ ND[k - 1] };
					double h{ ND[k] };
					double f{ 0.5e0 * (((g + z) / h) * ((g - z) / y) + y / h - h / y) };
					g = hypot(f, 1);
					f = x - (z / x) * z + (h / x) * (y / (f + copysign(g, f)) - h);
					// next qr transformation
					double c{ 1 };
					double s{ 1 };
					for (int i = l + 1; i <= k; i++) {
						g = ND[i];
						y = w[i];
						h = s * g;
						g = c * g;
						GivensRotation(c, s, f, h, ND[i - 1]);
						f = c * x + s * g;
						g = -s * x + c * g;
						h = y * s;
						y = y * c;
						if (V) V->ApplyGivens(c, s, i - 1, i, matrix_direction::mRight);
						GivensRotation(c, s, f, h, w[i - 1]);
						f = c * g + s * y;
						x = -s * g + c * y;
						if (U) U->ApplyGivens(c, s, i - 1, i, matrix_direction::mRight);
					}
					ND[l] = 0;
					ND[k] = f;
					w[k] = x;
				}
				else {
					// convergence
					if (Convergence)
						if (z < 0) {
							// w(k) is made non - negative
							w[k] = -z;
							if (V)
								for (int j = 0; j < nn; j++) (*V)(j, k) = -(*V)(j, k);
						}
				}
			} while (!Convergence && (its < 30));
		}
		// set error -- no convergence to a singular value after 30 iterations
		if (its >= 30) return k+1;
		return 0;
		// if successfull A = U * Sigma * V^T
	}

}
